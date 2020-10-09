from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef50, SSN_record
from retrobiocat_web.analysis.all_by_all_blast import AllByAllBlaster
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
import mongoengine as db
import networkx as nx
import time
import json
from rq.registry import StartedJobRegistry
from pathlib import Path
import os
import pandas as pd
from bson.binary import Binary
from collections import Counter
import random
import numpy as np
import palettable


class ClusterPositioner(object):

    def __init__(self, max_width=20):

        self.scale = 1000
        self.move = 0
        self.h_move = 0
        self.center = [0, 0]
        self.max_width = max_width * self.scale
        self.node_space = 200 * self.scale

    def measure_cluster(self, cluster):
        cluster_size = len(cluster) * self.node_space
        self.move = int(np.sqrt(cluster_size))
        if self.move > self.h_move:
            self.h_move = self.move

    def move_horizontal(self):
        self.center[0] += self.move

    def move_vertical(self):
        if self.center[0] > self.max_width:
            self.center[0] = 0
            self.center[1] += self.h_move*1.25
            self.h_move = 0

    def move_pos_dict(self, pos_dict):
        new_pos_dict = {}
        for key in pos_dict:
            new_pos_dict[key] = [pos_dict[key][0] + self.center[0], pos_dict[key][1] + self.center[1]]

        return new_pos_dict

class SSN_Visualiser(object):

    def __init__(self, enzyme_type, log_level=0):
        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
        self.node_metadata = self._find_uniref_metadata()

        self.edge_colour = {'color': 'lightgrey', 'opacity': 0.8}
        self.edge_width = 0.4
        self.uniref_border_width = 1
        self.uniref_border_colour = 'black'
        self.biocatdb_border_width = 3
        self.biocatdb_border_colour = 'darkred'
        self.border_width_selected = 4
        self.opacity = 0.9
        self.luminosity = 'bright'
        self.node_colour = f'rgba(5, 5, 168, {self.opacity})'
        self.node_size = 100
        self.node_shape = 'dot'

        self.log_level = log_level
        self.cluster_positioner = ClusterPositioner()

    def visualise(self, ssn, alignment_score, include_all_edges=True):
        graph = ssn.get_graph_filtered_edges(alignment_score)
        clusters = list(nx.connected_components(graph))
        clusters.sort(key=len, reverse=True)

        graph = self._add_cluster_node_colours(graph, clusters)
        pos_dict = self._get_cluster_positions(graph, clusters)

        if include_all_edges == True:
            nodes, edges = self._get_nodes_and_edges(graph, pos_dict, full_graph=ssn.graph)
        else:
            nodes, edges = self._get_nodes_and_edges(graph, pos_dict)
        return nodes, edges

    def _get_cluster_positions(self, graph, clusters):

        pos_dict = {}
        for i, cluster in enumerate(clusters):
            self.log(f"Getting layout for cluster {i+1} of {len(clusters)}")
            self.cluster_positioner.measure_cluster(cluster)
            self.cluster_positioner.move_horizontal()
            sub_graph = graph.subgraph(cluster)
            scale = 100+(35*len(cluster))

            if len(cluster) > 100:
                cluster_positions = nx.nx_pydot.pydot_layout(sub_graph, prog="sfdp")
            elif len(cluster) > 10:
                cluster_positions = nx.nx_pydot.pydot_layout(sub_graph, prog="neato")
            else:
                cluster_positions = nx.spring_layout(sub_graph, k=2, iterations=200, weight=None)

            cluster_positions = nx.rescale_layout_dict(cluster_positions, scale=scale)
            cluster_positions = self.cluster_positioner.move_pos_dict(cluster_positions)
            pos_dict.update(cluster_positions)
            self.cluster_positioner.move_horizontal()
            self.cluster_positioner.move_vertical()

        return pos_dict

    def _add_cluster_node_colours(self, graph, clusters):
        self.log(f"Colouring clusters.. (opacity={self.opacity})")
        colours = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
        for i, cluster in enumerate(clusters):
            c = colours.pop(0)
            colours.append(c)
            cc = f"rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},{self.opacity})"
            for node in cluster:
                graph.nodes[node]['colour'] = cc

        return graph

    def _get_nodes_and_edges(self, graph, pos_dict, full_graph=None):
        nodes = []
        edges = []
        for name in graph.nodes:
            colour = graph.nodes[name].get('colour', None)
            nodes.append(self._get_vis_node(name, pos_dict=pos_dict, colour=colour))

        if full_graph is None:
            for edge in graph.edges:
                weight = graph.get_edge_data(edge[0], edge[1], default={'weight': 0})['weight']
                edges.append(self._get_vis_edge(edge[0], edge[1], weight))
        else:
            for edge in full_graph.edges:
                weight = full_graph.get_edge_data(edge[0], edge[1], default={'weight': 0})['weight']
                edges.append(self._get_vis_edge(edge[0], edge[1], weight))

        nodes = self._sort_biocatdb_nodes_to_front(nodes)

        return nodes, edges

    def _get_vis_node(self, node_name, pos_dict=None, colour=None):
        if colour is None:
            colour = self.node_colour

        if 'UniRef50' in node_name:
            border = self.uniref_border_colour
            border_width = self.uniref_border_width
            node_type = 'uniref'
        else:
            border = self.biocatdb_border_colour
            border_width = self.biocatdb_border_width
            node_type = 'biocatdb'

        metadata = self.node_metadata.get(node_name, {})
        protein_name = metadata.get('protein_name', '')
        tax = metadata.get('tax', '')
        if protein_name != '':
            title = f"{protein_name} - {tax}"
        else:
            title = node_name

        node = {'id': node_name,
                'size': self.node_size,
                'borderWidth': border_width,
                'borderWidthSelected': self.border_width_selected,
                'color': {'background': colour,
                          'border': border,
                          'highlight': {'border': border}},
                'title': title,
                'label': '',
                'shape': self.node_shape,
                'node_type': node_type,
                'metadata': metadata}

        if pos_dict is not None:
            x, y = tuple(pos_dict.get(node_name, (0, 0)))
            node['x'] = x
            node['y'] = y

        return node

    def _get_vis_edge(self, edge_one, edge_two, weight):
        #weight = self.graph.get_edge_data(edge_one, edge_two, default={'weight': 0})['weight']
        edge = {'id': f"from {edge_one} to {edge_two}",
                'from': edge_one,
                'to': edge_two,
                'weight': weight,
                'width': self.edge_width,
                'color': self.edge_colour}
        return edge

    def _sort_biocatdb_nodes_to_front(self, vis_nodes):
        """ Returns vis_nodes with any nodes marked as node_type='biocatdb' at the front """

        biocatdb_nodes = []
        other_nodes = []

        for node in vis_nodes:
            if 'biocatdb' in node.get('node_type', ''):
                biocatdb_nodes.append(node)
            else:
                other_nodes.append(node)

        return other_nodes + biocatdb_nodes

    def _find_uniref_metadata(self):
        node_metadata = {}
        unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj).exclude('id', 'enzyme_type', 'sequence', "result_of_blasts_for")

        for seq_obj in unirefs:
            node_metadata[seq_obj.enzyme_name] = json.loads(seq_obj.to_json())

        return node_metadata

    def log(self, msg, level=1):
        if level >= self.log_level:
            print(f"SSN_Visualiser: {msg}")

class SSN(object):

    def __init__(self, enzyme_type, aba_blaster=None, log_level=0):

        self.graph = nx.Graph()

        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

        if aba_blaster is None:
            self.aba_blaster = AllByAllBlaster(enzyme_type, log_level=log_level)
        else:
            self.aba_blaster = aba_blaster

        self.node_metadata = {}

        self.log_level = log_level

        self.save_path = str(Path(__file__).parents[0]) + f'/analysis_data/ssn/{self.enzyme_type}'
        if not os.path.exists(self.save_path):
            os.mkdir(self.save_path)

        self.db_object = self._get_db_object()

        self.log(f"SSN object initialised for {enzyme_type}")

    def save(self):

        t0 = time.time()

        df_graph = nx.to_pandas_edgelist(self.graph)

        att_dict = {}
        for node in list(self.graph):
            att_dict[node] = self.graph.nodes[node]

        df_graph.to_csv(f"{self.save_path}/graph.csv")

        with open(f'{self.save_path}/attributes.json', 'wb') as outfile:
            outfile.write(json.dumps(att_dict).encode("utf-8"))

        t1 = time.time()

        self.log(f"Saved SSN for {self.enzyme_type} in {round(t1 - t0, 1)} seconds")

    def load(self, include_mutants=True, only_biocatdb=False):

        t0 = time.time()
        if not os.path.exists(f"{self.save_path}/graph.csv") or not os.path.exists(f"{self.save_path}/attributes.json"):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        df_graph = pd.read_csv(f"{self.save_path}/graph.csv")
        att_dict = json.load(open(f'{self.save_path}/attributes.json'))

        self.graph = nx.from_pandas_edgelist(df_graph, edge_attr=['weight'])

        # Nodes with no edges are not in edge list..
        for node in att_dict:
            if node not in self.graph.nodes:
                self._add_protein_node(node)

        if include_mutants is False:
            self.filter_out_mutants()
        if only_biocatdb is True:
            self.filer_out_uniref()

        nx.set_node_attributes(self.graph, att_dict)

        t1 = time.time()
        self.log(f"Loaded SSN for {self.enzyme_type} in {round(t1 - t0, 1)} seconds")

    def add_protein(self, seq_obj):
        """ Add the protein to the graph, along with any proteins which have alignments """

        self.log(f"Adding node - {seq_obj.enzyme_name} and making alignments..")
        t0 = time.time()

        name = seq_obj.enzyme_name
        self._add_protein_node(name, alignments_made=True)
        alignment_names, alignment_scores = self.aba_blaster.get_alignments(seq_obj)
        #self.graph.nodes[name]['attributes']['alignments_made'] = True

        count = 0
        for i, protein_name in enumerate(alignment_names):
            count += self._add_protein_node(protein_name)
            self._add_alignment_edge(seq_obj.enzyme_name, protein_name, alignment_scores[i])

        t1 = time.time()
        self.log(f"{count} new nodes made for alignments, with {len(alignment_names)} edges added")
        self.log(f"Protein {seq_obj.enzyme_name} processed in {round(t1-t0,0)} seconds")

    def add_multiple_proteins(self, list_seq_obj):
        for seq_obj in list_seq_obj:
            self.add_protein(seq_obj)

    def nodes_need_alignments(self, max_num=None):
        """ Return nodes which needs alignments making, maximum of max_num"""

        t0 = time.time()
        need_alignments = []
        count = 0
        for node in list(self.graph.nodes):
            if self.graph.nodes[node]['alignments_made'] == False:
                seq_obj = self._get_sequence_object(node)
                need_alignments.append(seq_obj)
                count += 1
                if count == max_num:
                    break

        t1 = time.time()
        self.log(f"Identified {count} nodes which need alignments making in {round(t1-t0,1)} seconds")

        return need_alignments

    def nodes_not_present(self, only_biocatdb=False, max_num=None):
        """ Return a list of enzymes which are not in the ssn """

        # Get a list of all sequence objects of enzyme type
        t0 = time.time()
        sequences = Sequence.objects(db.Q(enzyme_type=self.enzyme_type) &
                                     db.Q(sequence__ne="") &
                                     db.Q(sequence__ne=None) &
                                     db.Q(sequence_unavailable__ne=True))
        if only_biocatdb is True:
            seq_objects = list(sequences)
        else:
            unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj)
            seq_objects = list(sequences) + list(unirefs)

        # Get sequences not in nodes
        not_in_nodes = []
        for seq_obj in seq_objects:
            if seq_obj.enzyme_name not in list(self.graph.nodes):
                if seq_obj.sequence != None:
                    if len(seq_obj.sequence) > 12:
                        not_in_nodes.append(seq_obj)

        # Return only up to the maximum number of sequences
        if max_num != None:
            if len(not_in_nodes) > max_num:
                not_in_nodes = not_in_nodes[0:max_num]

        t1 = time.time()
        self.log(f"Identified {len(not_in_nodes)} {self.enzyme_type} proteins which need adding, in {round(t1 - t0, 1)} seconds")
        return not_in_nodes

    def remove_nonexisting_seqs(self):

        t0 = time.time()
        sequences = Sequence.objects(enzyme_type=self.enzyme_type).distinct('enzyme_name')
        unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj).distinct('enzyme_name')
        protein_names = list(sequences) + list(unirefs)
        count = 0
        for node in list(self.graph.nodes):
            if node not in protein_names:
                self.log(f"Node: {node} not in the database - removing")
                self.graph.remove_node(node)
                count += 1

        t1 = time.time()
        self.log(f"Identified {count} sequences which were in SSN but not in database, in {round(t1-t0,1)} seconds")

    def get_graph_filtered_edges(self, alignment_score):
        sub_graph = nx.Graph([(u, v, d) for u, v, d in self.graph.edges(data=True) if d['weight'] >= alignment_score])
        for node in self.graph.nodes:
            if node not in sub_graph.nodes:
                sub_graph.add_node(node)

        return sub_graph

    def filter_out_mutants(self):
        t0 = time.time()
        mutants = Sequence.objects(db.Q(enzyme_type=self.enzyme_type) &
                                   (db.Q(mutant_of__ne='') & db.Q(mutant_of__ne=None))).distinct('enzyme_name')

        for mutant in list(mutants):
            if mutant in self.graph.nodes:
                self.graph.remove_node(mutant)

        t1 = time.time()
        self.log(f'Filtered mutants from graph in {round(t1-t0,1)} seconds')

    def filer_out_uniref(self):
        t0 = time.time()
        for node in list(self.graph.nodes):
            if 'UniRef50' in node:
                self.graph.remove_node(node)

        t1 = time.time()
        self.log(f'Filtered uniref50 sequences from graph in {round(t1 - t0, 1)} seconds')

    def _add_protein_node(self, node_name, alignments_made=False):
        """ If a protein is not already in the graph, then add it """
        if 'UniRef50' in node_name:
            node_type = 'uniref'
        else:
            node_type = 'biocatdb'

        if node_name not in self.graph.nodes:
            self.graph.add_node(node_name, node_type=node_type,
                                alignments_made=alignments_made)
            return 1

        if alignments_made == True:
            self.graph.nodes[node_name]['alignments_made'] = True

        return 0

    def _add_alignment_edge(self, node_name, alignment_node_name, alignment_score):
        if node_name != alignment_node_name:
            self.graph.add_edge(node_name, alignment_node_name, weight=alignment_score)

    def log(self, msg, level=10):
        if level >= self.log_level:
            print("SSN: " + msg)

    @staticmethod
    def _get_sequence_object(enzyme_name):
        if 'UniRef50' in enzyme_name:
            return UniRef50.objects(enzyme_name=enzyme_name)[0]
        else:
            return Sequence.objects(enzyme_name=enzyme_name)[0]

    def _get_db_object(self):
        """ Either finds existing db entry for ssn of enzyme type, or makes a new one """

        query = SSN_record.objects(enzyme_type=self.enzyme_type_obj)
        if len(query) == 0:
            db_ssn = SSN_record(enzyme_type=self.enzyme_type_obj)
        else:
            db_ssn = query[0]

        return db_ssn

    def set_status(self, status):
        self.db_object.status = status
        self.db_object.save()




def task_expand_ssn(enzyme_type, log_level=1, max_num=200):
    current_app.app_context().push()

    aba_blaster = AllByAllBlaster(enzyme_type, log_level=log_level)
    aba_blaster.make_blast_db()

    ssn = SSN(enzyme_type, aba_blaster=aba_blaster, log_level=log_level)
    ssn.load()
    ssn.set_status('Chekcing SSN')
    ssn.remove_nonexisting_seqs()

    biocatdb_seqs = ssn.nodes_not_present(only_biocatdb=True, max_num=max_num)
    if len(biocatdb_seqs) != 0:
        ssn.set_status('Adding and aligning BioCatDB sequences')
        ssn.add_multiple_proteins(biocatdb_seqs)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    need_alignments = ssn.nodes_need_alignments(max_num=max_num)
    if len(need_alignments) != 0:
        ssn.set_status('Aligning sequences in SSN')
        ssn.add_multiple_proteins(need_alignments)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    not_present = ssn.nodes_not_present(max_num=max_num)
    if len(not_present) != 0:
        ssn.set_status('Adding UniRef sequences which are not yet present')
        ssn.add_multiple_proteins(not_present)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)

        return

    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    ssn.set_status('Complete')
    enz_type_obj.save()
    ssn.save()

def new_expand_ssn_job(enzyme_type):
    time.sleep(3)
    active_process_jobs = list(StartedJobRegistry(queue=current_app.alignment_queue).get_job_ids())
    active_process_jobs.extend(current_app.alignment_queue.job_ids)

    job_name = f"{enzyme_type}_expand_ssn"
    if job_name not in active_process_jobs:
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type, job_id=job_name)


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    aad_ssn = SSN('AAD', log_level=1)
    aad_ssn.load()

    aad_vis = SSN_Visualiser('AAD', log_level=1)
    nodes, edges = aad_vis.visualise(aad_ssn, 75)

    # 1. Set minimum number of nodes for a cluster
    # 2. Move down alignment score.  For each score where there is a different number of scores, visualise this.



