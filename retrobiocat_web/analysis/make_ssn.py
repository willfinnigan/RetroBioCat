from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef50, SSN_record
from retrobiocat_web.analysis.all_by_all_blast import AllByAllBlaster
from flask import current_app
import mongoengine as db
import networkx as nx
import time
import json
from pathlib import Path
import os
import pandas as pd
import palettable
import statistics
import dask.dataframe
import gc

class SSN_Cluster_Precalculator(object):

    def __init__(self, ssn, log_level=1):
        self.ssn = ssn

        self.start = 40
        self.end = 250
        self.step = 5
        self.min_cluster_size = 6

        self.log_level = log_level

    def precalulate(self, num=5, current_num_clusters=0):
        precalculated_nodes = {}
        num_at_alignment_score = {}
        identity_at_alignment_score = {}
        count = 0
        for alignment_score in range(self.start, self.end, self.step):
            visualiser = SSN_Visualiser(self.ssn.enzyme_type, log_level=0)
            clusters, graph = visualiser.get_clusters_and_subgraph(self.ssn, alignment_score)
            num_clusters = self._get_num_clusters(clusters)
            if (num_clusters > current_num_clusters) and (num_clusters != 0):
                count += 1
                self.log(f"Adding new set of {num_clusters} clusters at alignment score {alignment_score}")

                current_num_clusters = num_clusters

                pos_dict = visualiser.get_cluster_positions(graph, clusters)
                nodes, edges = visualiser.get_nodes_and_edges(graph, pos_dict)
                ident = self._get_ident(graph)

                num_at_alignment_score[str(alignment_score)] = num_clusters
                precalculated_nodes[str(alignment_score)] = nodes
                identity_at_alignment_score[str(alignment_score)] = ident

            if count == num:
                return precalculated_nodes, num_at_alignment_score, identity_at_alignment_score

        return precalculated_nodes, num_at_alignment_score, identity_at_alignment_score

    def _get_ident(self, graph):
        identities = []
        for edge in list(graph.edges(data=True)):
            identities.append(edge[2]['i'])

        if len(identities) >= 2:
            ident = [round(statistics.mean(identities), 2),
                     round(statistics.stdev(identities, statistics.mean(identities)), 2)]
        else:
            ident = [0, 0]

        return ident

    def _get_num_clusters(self, clusters):
        num = 0
        for cluster in clusters:
            if len(cluster) >= self.min_cluster_size:
                num += 1
        return num

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"Cluster_Precalc({self.ssn.enzyme_type}): {msg}")

class ClusterPositioner(object):

    def __init__(self, max_width=20):

        self.scale = 2000
        self.v_move = 0
        self.v_move_factor = 1.1
        self.max_width = max_width * self.scale
        self.gutter = 4000
        self.current_location = [0, 0]

    def get_cluster_dimensions(self, pos_dict):
        min_x, max_x = 0, 0
        min_y, max_y = 0, 0

        for node, pos in pos_dict.items():
            if pos[0] < min_x:
                min_x = pos[0]
            if pos[0] > max_x:
                max_x = pos[0]
            if pos[1] < min_y:
                min_y = pos[1]
            if pos[1] > max_y:
                max_y = pos[1]

        h_dim = max_x - min_x
        v_dim = max_y - min_y
        v_dim_move = 0
        if v_dim * 0.8 > self.v_move:
            self.v_move = v_dim * 0.8
            v_dim_move = v_dim * 0.2

        return h_dim, v_dim_move

    @staticmethod
    def move(pos_dict, move):
        new_pos_dict = {}
        for key in pos_dict:
            new_pos_dict[key] = [pos_dict[key][0] + move[0], pos_dict[key][1] + move[1]]
        return new_pos_dict

    @staticmethod
    def make_coords_positive(pos_dict):
        new_pos_dict = {}
        move_x, move_y = 0, 0
        for key, value in pos_dict.items():
            x, y = value[0], value[1]
            if x > move_x:
                move_x = x
            if y > move_y:
                move_y = y

        for key, value in pos_dict.items():
            x, y = value[0], value[1]
            x += move_x
            y += move_y
            new_pos_dict[key] = [x, y]

        return new_pos_dict

    def check_max_width(self):
        if self.current_location[0] > self.max_width:
            self.current_location[0] = 0
            self.current_location[1] += (self.v_move * self.v_move_factor) + self.gutter
            self.v_move = 0
            return True
        return False

    def position(self, pos_dict):

        h_dim, v_dim = self.get_cluster_dimensions(pos_dict)

        pos_dict = self.make_coords_positive(pos_dict)
        pos_dict = self.move(pos_dict, self.current_location)

        if self.check_max_width() == False:
            self.current_location[0] += self.gutter + h_dim
            self.current_location[1] += v_dim

        return pos_dict

    def round_positions(self, pos_dict, round_to=2):
        rounded_pos_dict = {}
        for key in pos_dict:
            rounded_pos_dict[key] = [round(pos_dict[key][0], round_to),
                                     round(pos_dict[key][1], round_to)]
        return rounded_pos_dict

class SSN_Visualiser(object):

    def __init__(self, enzyme_type, hidden_edges=True, log_level=0):
        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
        self.node_metadata = self._find_uniref_metadata()


        self.edge_colour = {'color': 'black'}
        self.edge_width = 4
        self.hidden_edges = hidden_edges
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

    def visualise(self, ssn, alignment_score, ):
        clusters, graph = self.get_clusters_and_subgraph(ssn, alignment_score)
        pos_dict = self.get_cluster_positions(graph, clusters)
        nodes, edges = self.get_nodes_and_edges(graph, pos_dict)

        return nodes, edges


    def get_cluster_positions(self, graph, clusters):

        pos_dict = {}
        for i, cluster in enumerate(clusters):
            self.log(f"Getting layout for cluster {i + 1} of {len(clusters)}")
            sub_graph = graph.subgraph(cluster)
            scale = 750 + (20 * len(cluster))

            #if len(cluster) > 200:
            #    cluster_positions = nx.nx_pydot.pydot_layout(sub_graph, prog="sfdp")
            #elif len(cluster) > 10:
            #    cluster_positions = nx.nx_pydot.pydot_layout(sub_graph, prog="neato")
            #else:
            cluster_positions = nx.spring_layout(sub_graph, k=1, iterations=200, scale=scale, weight=None)

            cluster_positions = self.cluster_positioner.position(cluster_positions)
            cluster_positions = self.cluster_positioner.round_positions(cluster_positions, round_to=0)
            pos_dict.update(cluster_positions)

        return pos_dict

    def get_clusters_and_subgraph(self, ssn, alignment_score):
        graph = ssn.get_graph_filtered_edges(alignment_score)
        clusters = list(nx.connected_components(graph))
        clusters.sort(key=len, reverse=True)
        graph = self._add_cluster_node_colours(graph, clusters)

        return clusters, graph

    def _add_cluster_node_colours(self, graph, clusters):
        self.log(f"Colouring clusters.. (opacity={self.opacity})")
        colours = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
        for i, cluster in enumerate(clusters):
            c = colours.pop(0)
            colours.append(c)
            cc = f"rgba({int(c[0] * 255)},{int(c[1] * 255)},{int(c[2] * 255)},{self.opacity})"
            for node in cluster:
                graph.nodes[node]['colour'] = cc

        return graph

    def get_nodes_and_edges(self, graph, pos_dict, full_graph=None):
        nodes = []
        edges = []
        for name in graph.nodes:
            colour = graph.nodes[name].get('colour', None)
            nodes.append(self._get_vis_node(name, pos_dict=pos_dict, colour=colour))

        if full_graph is None:
            for edge in graph.edges:
                weight = graph.get_edge_data(edge[0], edge[1], default={'weight': 0})['weight']
                edges.append(self.get_vis_edge(edge[0], edge[1], weight))
        else:
            for edge in full_graph.edges:
                weight = full_graph.get_edge_data(edge[0], edge[1], default={'weight': 0})['weight']
                edges.append(self.get_vis_edge(edge[0], edge[1], weight))

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
                'node_type': node_type}

        if pos_dict is not None:
            x, y = tuple(pos_dict.get(node_name, (0, 0)))
            node['x'] = x
            node['y'] = y

        return node

    def get_vis_edge(self, edge_one, edge_two, weight):
        # weight = self.graph.get_edge_data(edge_one, edge_two, default={'weight': 0})['weight']
        edge = {'id': f"from {edge_one} to {edge_two}",
                'from': edge_one,
                'to': edge_two,
                'hidden': self.hidden_edges,
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
        unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj).exclude('id', 'enzyme_type', 'sequence',
                                                                             "result_of_blasts_for")

        for seq_obj in unirefs:
            node_metadata[seq_obj.enzyme_name] = json.loads(seq_obj.to_json())

        return node_metadata

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"SSN_Visualiser: {msg}")

class SSN(object):

    def __init__(self, enzyme_type, aba_blaster=None, log_level=0):

        self.graph = nx.Graph()

        self.min_score = 40

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


        df_graph.to_csv(f"{self.save_path}/graph.csv", index=False)

        with open(f'{self.save_path}/attributes.json', 'wb') as outfile:
            outfile.write(json.dumps(att_dict).encode("utf-8"))

        t1 = time.time()

        self.log(f"Saved SSN for {self.enzyme_type} in {round(t1 - t0, 1)} seconds")

    def load(self, include_mutants=True, only_biocatdb=False, mode='pandas'):

        t0 = time.time()
        if not os.path.exists(f"{self.save_path}/graph.csv") or not os.path.exists(f"{self.save_path}/attributes.json"):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        #df_graph = pd.read_csv(f"{self.save_path}/graph.csv")

        if mode == 'dask':
            df_graph = dask.dataframe.read_csv(f"{self.save_path}/graph.csv")
        elif mode == 'dask_pandas':
            df_graph = dask.dataframe.read_csv(f"{self.save_path}/graph.csv")
            df_graph = df_graph.compute()
        else:
            df_graph = pd.read_csv(f"{self.save_path}/graph.csv")
        att_dict = json.load(open(f'{self.save_path}/attributes.json'))
        t1 = time.time()

        self.graph = nx.from_pandas_edgelist(df_graph, edge_attr=['weight', 'i'])

        t2 = time.time()

        # Nodes with no edges are not in edge list..
        for node in list(att_dict.keys()):
            if node not in list(self.graph.nodes):
                self._add_protein_node(node)

        if include_mutants is False:
            self.filter_out_mutants()
        if only_biocatdb is True:
            self.filer_out_uniref()

        nx.set_node_attributes(self.graph, att_dict)

        t3 = time.time()
        self.log(f"Loaded SSN for {self.enzyme_type} in {round(t3 - t0, 1)} seconds (mode = {mode})")
        self.log(f"- {round(t1 - t0, 1)} seconds to load dataframes")
        self.log(f"- {round(t2 - t1, 1)} seconds to load actual ssn from edge list")
        self.log(f"- {round(t3 - t2, 1)} seconds for attributes")

    def add_protein(self, seq_obj):
        """ Add the protein to the graph, along with any proteins which have alignments """

        self.log(f"Adding node - {seq_obj.enzyme_name} and making alignments..")
        t0 = time.time()

        name = seq_obj.enzyme_name
        self._add_protein_node(name, alignments_made=True)
        alignment_names, alignment_scores, identities, coverages = self.aba_blaster.get_alignments(seq_obj)
        # self.graph.nodes[name]['attributes']['alignments_made'] = True

        count = 0
        for i, protein_name in enumerate(alignment_names):
            count += self._add_protein_node(protein_name)
            self._add_alignment_edge(seq_obj.enzyme_name, protein_name, alignment_scores[i], identities[i],
                                     coverages[i])

        t1 = time.time()
        seq_obj.alignments_made = True
        seq_obj.save()

        self.log(f"{count} new nodes made for alignments, with {len(alignment_names)} edges added")
        self.log(f"Protein {seq_obj.enzyme_name} processed in {round(t1 - t0, 0)} seconds")

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
        self.log(f"Identified {count} nodes which need alignments making in {round(t1 - t0, 1)} seconds")

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
        self.log(
            f"Identified {len(not_in_nodes)} {self.enzyme_type} proteins which need adding, in {round(t1 - t0, 1)} seconds")
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
        self.log(f"Identified {count} sequences which were in SSN but not in database, in {round(t1 - t0, 1)} seconds")

    def remove_seqs_marked_with_no_alignments(self):

        t0 = time.time()
        enz_type_q = db.Q(enzyme_type=self.enzyme_type)
        align_q = db.Q(alignments_made__ne=True)
        protein_names = Sequence.objects(enz_type_q & align_q).distinct('enzyme_name')

        count = 0
        for node in list(self.graph.nodes):
            if node in protein_names:
                self.log(f"Node: {node} marked as having no alignments made - removing")
                self.graph.remove_node(node)
                count += 1

        t1 = time.time()
        self.log(f"Identified {count} sequences which were in SSN but not are marked has not having their alignments made, in {round(t1 - t0, 1)} seconds")

    def get_graph_filtered_edges(self, alignment_score, min_ident=0):

        sub_graph = nx.Graph([(u, v, d) for u, v, d in self.graph.edges(data=True) if d['weight'] >= alignment_score and d['i'] >= min_ident])
        for node in self.graph.nodes:
            if node not in sub_graph.nodes:
                sub_graph.add_node(node)

        return sub_graph

    def delete_edges_below_alignment_score(self, alignment_score, min_ident=0.0):
        self.log(f"Removing edges below alignment score: {alignment_score}, min identity: {min_ident}")
        edge_list_to_remove = []
        for edge in self.graph.edges(data=True):
            if edge[2]['weight'] < alignment_score or edge[2]['weight'] < min_ident:
                edge_list_to_remove.append(edge)

        for edge in edge_list_to_remove:
            self.graph.remove_edge(edge[0], edge[1])

    def filter_out_mutants(self):
        t0 = time.time()
        mutants = Sequence.objects(db.Q(enzyme_type=self.enzyme_type) &
                                   (db.Q(mutant_of__ne='') & db.Q(mutant_of__ne=None))).distinct('enzyme_name')

        for mutant in list(mutants):
            if mutant in self.graph.nodes:
                self.graph.remove_node(mutant)

        t1 = time.time()
        self.log(f'Filtered mutants from graph in {round(t1 - t0, 1)} seconds')

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

    def _add_alignment_edge(self, node_name, alignment_node_name, alignment_score, i, c):
        if node_name != alignment_node_name:
            if alignment_score > self.min_score:
                self.graph.add_edge(node_name, alignment_node_name, weight=alignment_score, i=i)

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

    def clear_position_information(self):
        self.log(f"Clearing old node position information", level=1)
        if self.db_object.num_at_alignment_score != {}:
            self.db_object.num_at_alignment_score = {}
            self.db_object.pos_at_alignment_score = {}
            self.db_object.precalcuated_vis = {}
            self.db_object.identity_at_alignment_score = {}
            self.db_object.save()

class SSN_quickload(object):

    def __init__(self, enzyme_type, log_level=0):
        self.enzyme_type = enzyme_type
        self.save_path = str(Path(__file__).parents[0]) + f'/analysis_data/ssn/{self.enzyme_type}'
        #self.ssn = SSN(enzyme_type, log_level=log_level)
        self.log_level = log_level
        self.df = None
        self.vis = SSN_Visualiser(enzyme_type, hidden_edges=False, log_level=log_level)

    def load_df(self):
        t0 = time.time()
        if not os.path.exists(f"{self.save_path}/graph.csv") or not os.path.exists(f"{self.save_path}/attributes.json"):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        #self.df = dask.dataframe.read_csv(f"{self.save_path}/graph.csv")
        #self.df = self.df.compute()
        self.df = pd.read_csv(f"{self.save_path}/graph.csv")
        t1 = time.time()

        self.log(f"Loaded df in {round(t1-t0, 1)} seconds")

    def get_edges(self, selected_node, alignment_score):
        t0 = time.time()
        edges = []
        df = self.df[(self.df['source'] == selected_node) | (self.df['target'] == selected_node)]
        df = df[df['weight'] >= alignment_score]

        for index, row in df.iterrows():
            edges.append(self.vis.get_vis_edge(row['target'], row['source'], row['weight']))
        t1 = time.time()

        self.log(f"Retrieved edges for {self.enzyme_type} node at alignment score {alignment_score} in {round(t1-t0, 1)} seconds")

        return edges

    def get_connected_nodes(self, list_nodes, alignment_score):
        t0 = time.time()
        connected_nodes = set()

        df = self.df[((self.df['source'].isin(list_nodes)) | (self.df['target'].isin(list_nodes)))]
        df = df[df['weight'] >= alignment_score]

        for index, row in df.iterrows():
            if row['target'] not in list_nodes:
                connected_nodes.add(row['target'])
            if row['source'] not in list_nodes:
                connected_nodes.add(row['source'])

        connected_nodes = list(connected_nodes)
        t1 = time.time()
        self.log(f"Retrieved {len(connected_nodes)} connected nodes for {self.enzyme_type} node at alignment score {alignment_score} in {round(t1 - t0, 1)} seconds")

        return connected_nodes




    def log(self, msg, level=1):
        if self.log_level >= level:
            print("SSN_quickload: " + msg)



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    ql = SSN_quickload('IRED', log_level=1)
    ql.load_df()
    edges = ql.get_edges('UniRef50_Q2TW47', 45)
    nodes = ql.get_connected_nodes(['UniRef50_Q2TW47'], 45)


    #ssn = SSN('IRED', log_level=1)
    #ssn.load(mode='pandas')
    #ssn.load(mode='dask')
    #ssn.load(mode='dask_pandas')
    #ssn.delete_edges_below_alignment_score(40)
    #ssn.save()

    """
    print(f"Num edges default = {len(ssn.graph.edges)}")

    
    """

    #ssn.delete_edges_below_alignment_score(40, min_ident=0.35)
    #print(f"Num edges = {len(ssn.graph.edges)}")

    #ssn.delete_edges_below_alignment_score(45)
    #print(f"Num edges = {len(ssn.graph.edges)}")

    """
    filtered_graph = ssn.get_graph_filtered_edges(40)
    print(f"Num edges alignment 40 = {len(filtered_graph.edges)}")

    filtered_graph = ssn.get_graph_filtered_edges(45)
    print(f"Num edges alignment 45 = {len(filtered_graph.edges)}")

    filtered_graph = ssn.get_graph_filtered_edges(50)
    print(f"Num edges alignment 50 = {len(filtered_graph.edges)}")

    filtered_graph = ssn.get_graph_filtered_edges(55)
    print(f"Num edges alignment 55 = {len(filtered_graph.edges)}")
    
    """

    #for i in range(20,100,10):
     #   new_graph = ssn.get_graph_filtered_edges(i)
      #  num_edges = len(new_graph.edges)
       # print(f"Score = {i}, num edges = {num_edges}")

    # 1. Set minimum number of nodes for a cluster
    # 2. Move down alignment score.  For each score where there is a different number of scores, visualise this.
