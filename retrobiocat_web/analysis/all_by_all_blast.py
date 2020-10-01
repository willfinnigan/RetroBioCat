import os
from pathlib import Path
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef90, Alignment
import mongoengine as db
import subprocess as sp
import shutil
from pathlib import Path
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
import numpy as np
import shutil
import datetime
import time
import networkx as nx

ALLBYALL_BLAST_FOLDER = str(Path(__file__).parents[0]) + '/all_by_all_blast'

def make_fasta(enzyme_type, directory):
    """ Create a fasta file containing all the sequences of an enzyme type """

    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(sequence__ne=""))
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    bioinf_seqs = UniRef90.objects(db.Q(enzyme_type=enzyme_type_obj))

    with open(f"{directory}/{enzyme_type}.fasta", 'w') as file:
        for seq in list(seqs) + list(bioinf_seqs):
            name = seq.enzyme_name
            seq = seq.sequence.replace('\n', '')

            file.write(f'>{name}\n')
            file.write(f"{seq}\n")

def make_blast_db_for_enzyme_type(enzyme_type):
    """ Create a blast database for a single enzyme type """

    directory = f"{ALLBYALL_BLAST_FOLDER}/{enzyme_type}"
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.mkdir(directory)

    make_fasta(enzyme_type, directory)

    command = f"makeblastdb -in {directory}/{enzyme_type}.fasta -dbtype prot"
    sp.run(command, shell=True)

def run_blast(fasta_to_blast, database):
    """ Run blast and parse output into a biopython blast record """

    output = NcbiblastpCommandline(query=fasta_to_blast, db=database, outfmt=5)()[0]

    return output

def make_single_seq_fasta(name, sequence, enzyme_type):
    """ Make a fasta file for a single sequence in the directory """

    directory = f"{ALLBYALL_BLAST_FOLDER}/{enzyme_type}"
    fasta_path = f"{directory}/{name}.fasta"
    with open(fasta_path, 'w') as file:
        sequence = sequence.replace('\n', '')
        file.write(f'>{name}\n')
        file.write(f"{sequence}\n")

    return fasta_path

def get_sequence_object(name):
    seq_query = Sequence.objects(enzyme_name=name)
    if len(seq_query) != 0:
        return seq_query[0]

    seq_query = UniRef90.objects(enzyme_name=name)
    if len(seq_query) != 0:
        return seq_query[0]

    return None

def create_alignments_from_blast_record(seq_obj, blast_record, enzyme_type,
                                        max_e=0.0005, min_identity=0.4, min_coverage=0.7):
    """ Create alignment database entries from the parsed blast record """

    alignment_count = 0
    query_length = len(seq_obj.sequence)
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    for alignment in blast_record.alignments:
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]
            sbjct_name = alignment.title.replace(f"{alignment.hit_id} ", "")
            sbjct_obj = get_sequence_object(sbjct_name)

            coverage = hsp.align_length / query_length
            identity = hsp.identities / hsp.align_length

            """alignment_query = Alignment.objects(db.Q(proteins__in=[seq_obj])
                                                & db.Q(proteins__in=[sbjct_obj])
                                                & db.Q(enzyme_type=enzyme_type_obj))"""
            if (hsp.expect <= max_e) and (identity >= min_identity) and (coverage >= min_coverage):

                sbjct_length = len(sbjct_obj.sequence)
                bitscore = hsp.bits
                x = 2 - bitscore * (query_length / sbjct_length)
                alignment_score = np.log10(-x)
                e_value = hsp.expect

                """
                if seq_obj.enzyme_name != sbjct_obj.enzyme_name:
                    alignment = Alignment(enzyme_type=enzyme_type_obj,
                                          proteins=[seq_obj, sbjct_obj],
                                          coverage=coverage,
                                          bitscore=bitscore,
                                          identity=identity,
                                          alignment_score=alignment_score,
                                          e_value=e_value)
                    alignment.save()
                    alignment_count += 1
                """

def do_all_by_all_blast(enzyme_type, make_alignment=True, time_func=False):
    """ Loop over sequences of an enzyme type creating alignments between all sequences of that type """

    directory = f"{ALLBYALL_BLAST_FOLDER}/{enzyme_type}"
    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(sequence__ne=""))
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    bioinf_seqs = UniRef90.objects(db.Q(enzyme_type=enzyme_type_obj))
    database = f"{directory}/{enzyme_type}.fasta"
    graph = nx.Graph()

    seqs_to_iterate = list(seqs)+list(bioinf_seqs)
    t0 = time.time()
    for i, seq in enumerate(seqs_to_iterate):
        if seq.alignments_made is None:
            fasta_path = make_single_seq_fasta(seq.enzyme_name, seq.sequence, enzyme_type)
            output = run_blast(fasta_path, database)
            os.remove(fasta_path)
            blast_record = NCBIXML.read(StringIO(output))
            if make_alignment == True:
                create_alignments_from_blast_record(seq, blast_record, enzyme_type)
                seq.alignments_made = datetime.datetime.now()
                seq.save()

        if (i % 50 == 0) and time_func is True:
            t1 = time.time()
            print(f"Time per 50 blasts is {round(t1 - t0, 1)} seconds")
            t0 = time.time()
        elif time_func is False:
            print(f"{seq.enzyme_name} blasted against all {enzyme_type}s - ({i} of {len(seqs_to_iterate)})")

def single_enz_all_by_all_blast(seq_name, seq_type, enzyme_type, make_alignment=True):

    directory = f"{ALLBYALL_BLAST_FOLDER}/{enzyme_type}"
    database = f"{directory}/{enzyme_type}.fasta"

    if seq_type == 'biocatdb':
        seq_obj = Sequence.objects(enzyme_name=seq_name)[0]
    else:
        seq_obj = UniRef90.objects(enzyme_name=seq_name)[0]

    if seq_obj.alignments_made is None:
        fasta_path = make_single_seq_fasta(seq_obj.enzyme_name, seq_obj.sequence, enzyme_type)
        os.remove(fasta_path)

        output = run_blast(fasta_path, database)
        os.remove(fasta_path)
        if make_alignment == True:
            blast_record = NCBIXML.read(StringIO(output))
            create_alignments_from_blast_record(seq_obj, blast_record, enzyme_type)
            seq_obj.alignments_made = datetime.datetime.now()
            seq_obj.save()

        print(f"{seq_obj.enzyme_name} blasted against all {enzyme_type}s")


if __name__ == "__main__":

    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    #Alignment.drop_collection()

    enzyme_type = 'AADH'
    #make_blast_db_for_enzyme_type(enzyme_type)
    do_all_by_all_blast(enzyme_type, make_alignment=False, time_func=True)


