import os
from pathlib import Path
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef90, Alignments
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
from decimal import Decimal

ALLBYALL_BLAST_FOLDER = str(Path(__file__).parents[0]) + '/all_by_all_blast'

def calc_alignment_score(bitscore, query_length, subject_length):
    two = Decimal(2)
    bitscore = Decimal(bitscore)
    x = np.power(two, -bitscore) * query_length * subject_length
    alignment_score = -np.log10(x)
    return alignment_score

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

def get_alignments_obj(enzyme_type):
    query = Alignments.objects(enzyme_type=enzyme_type)
    if len(query) == 0:
        enzyme_type_object = EnzymeType(enzyme_type=enzyme_type)[0]
        alignments_obj = Alignments(enzyme_type=enzyme_type_object)

def blast_seq_against_local_enzyme_db(seq_obj, enzyme_type,
                                      max_e=0.0005, min_identity=0.3, min_coverage=0.7):

    directory = f"{ALLBYALL_BLAST_FOLDER}/{enzyme_type}"
    database = f"{directory}/{enzyme_type}.fasta"
    query_length = len(seq_obj.sequence)

    fasta_path = make_single_seq_fasta(seq_obj.enzyme_name, seq_obj.sequence, enzyme_type)
    output = run_blast(fasta_path, database)
    blast_record = NCBIXML.read(StringIO(output))
    os.remove(fasta_path)

    alignments = []
    for alignment in blast_record.alignments:
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]

            sbjct_name = alignment.title.replace(f"{alignment.hit_id} ", "")
            sbjct_obj = get_sequence_object(sbjct_name)

            coverage = hsp.align_length / query_length
            identity = hsp.identities / hsp.align_length

            if (hsp.expect <= max_e) and (identity >= min_identity) and (coverage >= min_coverage):
                pass

