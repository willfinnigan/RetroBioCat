import os
from pathlib import Path
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef50
import mongoengine as db
import subprocess as sp
import shutil
from pathlib import Path
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
import numpy as np
import shutil
import time

from decimal import Decimal

class AllByAllBlaster(object):
    """ For performing all by all blasts for a given enzyme type """

    def __init__(self, enzyme_type, print_log=False, num_threads=2):
        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

        self.all_by_all_blast_folder = str(Path(__file__).parents[0]) + '/analysis_data/all_by_all_blast'
        self.directory = f"{self.all_by_all_blast_folder}/{enzyme_type}"
        self.database = f"{self.directory}/{enzyme_type}.fasta"
        self.cdhit_output = f"{self.directory}/cd_hit"
        self.num_threads = num_threads
        self.max_alignments = 10000

        self.max_e = 5
        self.min_coverage = 0.8
        self.min_identity = 0.3

        self.print_log = print_log

    def make_blast_db(self):
        """ Create a blast database for a single enzyme type """
        self.log(f"Making blast database for {self.enzyme_type}")

        if os.path.exists(self.directory):
            shutil.rmtree(self.directory)
        os.mkdir(self.directory)

        self._make_db_fasta()

        command = f"makeblastdb -in {self.directory}/{self.enzyme_type}.fasta -dbtype prot"
        sp.run(command, shell=True)

    def get_alignments(self, seq_obj):
        """ Get the alignments for a given sequence object (either Sequence or UniRef90)"""

        self.log(f" - Getting alignments for sequence: {seq_obj.enzyme_name}..")
        if (seq_obj.sequence != None):
            if len(seq_obj.sequence) > 12:
                blast_record = self._blast_seq(seq_obj)
                alignment_names, alignment_scores = self._process_blast_record(seq_obj, blast_record)
                self.log(f"{len(alignment_names)} made.")
                return alignment_names, alignment_scores

        return [], []

    def get_clusters(self, identity):

        t0 = time.time()
        if not os.path.exists(self.cdhit_output):
            os.mkdir(self.cdhit_output)

        out_file = f"{self.cdhit_output}/{identity}_identity"
        cmd = f"cd-hit -i {self.database} -o {out_file} -c {identity}"

        sp.run(cmd, shell=True)

        t1 = time.time()
        self.log(f"Clustered sequences for {self.enzyme_type} at {identity} identity in {round(t1-t0,1)} seconds")

    def _blast_seq(self, seq_obj):
        """ Run blast and parse output into a biopython blast record """
        t0 = time.time()
        protein_seq = seq_obj.sequence.replace('\n', '')
        protein_name = seq_obj.enzyme_name
        fasta_path = f"{self.directory}/{protein_name}.fasta"
        with open(fasta_path, 'w') as file:
            file.write(f'>{protein_name}\n')
            file.write(f"{protein_seq}\n")

        output = NcbiblastpCommandline(query=fasta_path, db=self.database, evalue=self.max_e,
                                       num_threads=self.num_threads, outfmt=5, num_alignments=self.max_alignments)()[0]
        blast_record = NCBIXML.read(StringIO(output))

        os.remove(fasta_path)
        t1 = time.time()
        self.log(f"sequence blasted in {round(t1-t0,0)} seconds")

        return blast_record

    def _process_blast_record(self, query_object, blast_record):
        """ Process the blast record generating a list of alignments"""

        t0 = time.time()
        query_length = len(query_object.sequence)

        alignment_names = []
        alignment_scores = []
        for alignment in blast_record.alignments:
            if len(alignment.hsps) == 1:
                hsp = alignment.hsps[0]

                subject_name = alignment.title.replace(f"{alignment.hit_id} ", "")

                coverage = hsp.align_length / query_length
                identity = hsp.identities / hsp.align_length

                if (subject_name != query_object.enzyme_name) and (identity >= self.min_identity) and (coverage >= self.min_coverage):
                    score = self._calc_alignment_score(hsp.bits, query_length, hsp.align_length)
                    alignment_names.append(subject_name)
                    alignment_scores.append(score)
        t1 = time.time()
        self.log(f"processed alignments in {round(t1-t0,0)} seconds")

        return alignment_names, alignment_scores

    def _make_db_fasta(self):
        """ Create a fasta file containing all the sequences of an enzyme type """

        seqs = Sequence.objects(db.Q(enzyme_type=self.enzyme_type) &
                                db.Q(sequence__ne="") &
                                db.Q(sequence__ne=None) &
                                db.Q(sequence_unavailable__ne=True))
        bioinf_seqs = UniRef50.objects(db.Q(enzyme_type=self.enzyme_type_obj))

        with open(f"{self.directory}/{self.enzyme_type}.fasta", 'w') as file:
            for seq in list(seqs) + list(bioinf_seqs):
                name = seq.enzyme_name
                seq = seq.sequence.replace('\n', '')

                file.write(f'>{name}\n')
                file.write(f"{seq}\n")

    def log(self, msg):
        if self.print_log == True:
            print("ABA_Blaster: " + msg)

    @staticmethod
    def _calc_alignment_score(bitscore, query_length, subject_length):
        two = Decimal(2)
        bitscore = Decimal(bitscore)
        x = np.power(two, -bitscore) * query_length * subject_length
        alignment_score = int(-np.log10(x))
        return alignment_score

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()


    seq_obj = Sequence.objects(enzyme_type='AAD')[0]
    etb = AllByAllBlaster('AAD', print_log=True)
    etb.make_blast_db()
    alignment_names, alignment_scores = etb.get_alignments(seq_obj)

