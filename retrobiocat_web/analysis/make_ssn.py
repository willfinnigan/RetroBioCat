from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef90, DB_SSN
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

class SSN(object):

    def __init__(self, enzyme_type, aba_blaster=None, print_log=False):
        self.graph = nx.Graph()
        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

        if aba_blaster is None:
            self.aba_blaster = AllByAllBlaster(enzyme_type, print_log=print_log)
        else:
            self.aba_blaster = aba_blaster

        self.print_log = print_log

        self.save_path = str(Path(__file__).parents[0]) + f'/analysis_data/ssn/{self.enzyme_type}'
        if not os.path.exists(self.save_path):
            os.mkdir(self.save_path)

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

        self.log(f"Saved SSN for {self.enzyme_type} in {round(t1-t0,1)} seconds")

    def load(self):

        t0 = time.time()
        if not os.path.exists(f"{self.save_path}/graph.csv") or not os.path.exists(f"{self.save_path}/attributes.json"):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        df_graph = pd.read_csv(f"{self.save_path}/graph.csv")
        att_dict = json.load(open(f'{self.save_path}/attributes.json'))

        self.graph = nx.from_pandas_edgelist(df_graph, edge_attr=['weight'])
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
        sequences = Sequence.objects(enzyme_type=self.enzyme_type)
        if only_biocatdb is True:
            seq_objects = list(sequences)
        else:
            unirefs = UniRef90.objects(enzyme_type=self.enzyme_type_obj)
            seq_objects = list(sequences) + list(unirefs)

        # Get sequences not in nodes
        not_in_nodes = []
        for seq_obj in seq_objects:
            if seq_obj.enzyme_name not in list(self.graph.nodes):
                not_in_nodes.append(seq_obj)

        # Return only up to the maximum number of sequences
        if max_num != None:
            if len(not_in_nodes) > max_num:
                not_in_nodes = not_in_nodes[0:max_num]

        t1 = time.time()
        self.log(f"Identified {len(not_in_nodes)} {self.enzyme_type} proteins which need adding, in {round(t1 - t0, 1)} seconds")
        return not_in_nodes

    def remove_nonexisting_seqs(self):
        sequences = Sequence.objects(enzyme_type=self.enzyme_type).distinct('enzyme_name')
        unirefs = UniRef90.objects(enzyme_type=self.enzyme_type_obj).distinct('enzyme_name')
        protein_names = list(sequences) + list(unirefs)

        for node in list(self.graph.nodes):
            if node not in protein_names:
                self.log(f"Node: {node} not in the database - removing")
                self.graph.remove_node(node)

    def _add_protein_node(self, node_name, alignments_made=False):
        """ If a protein is not already in the graph, then add it """
        if 'UniRef90' in node_name:
            node_type = 'uniref'
        else:
            node_type = 'biocatdb'

        if node_name not in list(self.graph.nodes):
            self.graph.add_node(node_name, node_type=node_type, alignments_made=alignments_made)
            return 1

        if alignments_made == True:
            self.graph.nodes[node_name]['alignments_made'] = True

        return 0

    def _add_alignment_edge(self, node_name, alignment_node_name, alignment_score):
        self.graph.add_edge(node_name, alignment_node_name, weight=alignment_score)

    def log(self, msg):
        if self.print_log == True:
            print("SSN: " + msg)

    @staticmethod
    def _get_sequence_object(enzyme_name):
        if 'UniRef90' in enzyme_name:
            return UniRef90.objects(enzyme_name=enzyme_name)[0]
        else:
            return Sequence.objects(enzyme_name=enzyme_name)[0]

def task_expand_ssn(enzyme_type, print_log=True):
    current_app.app_context().push()

    ssn = SSN(enzyme_type, print_log=print_log)
    ssn.load()
    ssn.remove_nonexisting_seqs()

    max_num = 200

    biocatdb_seqs = ssn.nodes_not_present(only_biocatdb=True, max_num=max_num)
    if len(biocatdb_seqs) != 0:
        ssn.add_multiple_proteins(biocatdb_seqs)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    need_alignments = ssn.nodes_need_alignments(max_num=max_num)
    if len(need_alignments) != 0:
        ssn.add_multiple_proteins(need_alignments)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    not_present = ssn.nodes_not_present(max_num=max_num)
    if len(not_present) != 0:
        ssn.add_multiple_proteins(need_alignments)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)

        return

    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    enz_type_obj.bioinformatics_status = 'Idle'
    enz_type_obj.save()

def new_expand_ssn_job(enzyme_type):

    active_process_jobs = list(StartedJobRegistry(queue=current_app.alignment_queue).get_job_ids())
    active_process_jobs.extend(current_app.alignment_queue.job_ids)

    job_name = f"{enzyme_type}_expand_ssn"
    if job_name not in active_process_jobs:
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type, job_id=job_name)



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    aad_ssn = SSN('AAD', print_log=True)
    aad_ssn.load()

    biocatdb_seqs = aad_ssn.nodes_not_present(only_biocatdb=True, max_num=10)
    aad_ssn.add_multiple_proteins(biocatdb_seqs)

    need_alignments = aad_ssn.nodes_need_alignments(max_num=10)
    aad_ssn.add_multiple_proteins(need_alignments)

    aad_ssn.save()


