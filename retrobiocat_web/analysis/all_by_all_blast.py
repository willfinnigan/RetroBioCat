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
    blast_record = NCBIXML.read(StringIO(output))
    return blast_record

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

def create_alignments_from_blast_record(seq_obj, blast_record, enzyme_type):
    """ Create alignment database entries from the parsed blast record """

    query_length = len(seq_obj.sequence)
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    for alignment in blast_record.alignments:
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]
            sbjct_name = alignment.title.replace(f"{alignment.hit_id} ", "")
            sbjct_obj = get_sequence_object(sbjct_name)

            alignment_query = Alignment.objects(db.Q(proteins__in=[seq_obj])
                                                & db.Q(proteins__in=[sbjct_obj])
                                                & db.Q(enzyme_type=enzyme_type_obj))
            if len(alignment_query) == 0:
                coverage = hsp.align_length / query_length
                identity = hsp.identities / hsp.align_length
                sbjct_length = len(sbjct_obj.sequence)
                bitscore = hsp.bits
                x = 2 - bitscore * (query_length * sbjct_length)
                alignment_score = np.log10(-x)
                e_value = hsp.expect

                alignment = Alignment(enzyme_type=enzyme_type_obj,
                                      proteins=[seq_obj,sbjct_obj],
                                      coverage=coverage,
                                      bitscore=bitscore,
                                      identity=identity,
                                      alignment_score=alignment_score,
                                      e_value=e_value)
                alignment.save()

def do_all_by_all_blast(enzyme_type):
    """ Loop over sequences of an enzyme type creating alignments between all sequences of that type """

    directory = f"{ALLBYALL_BLAST_FOLDER}/{enzyme_type}"
    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(sequence__ne=""))
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    bioinf_seqs = UniRef90.objects(db.Q(enzyme_type=enzyme_type_obj))
    database = f"{directory}/{enzyme_type}.fasta"

    for seq in list(seqs)+list(bioinf_seqs):
        fasta_path = make_single_seq_fasta(seq.enzyme_name, seq.sequence, enzyme_type)
        blast_record = run_blast(fasta_path, database)
        os.remove(fasta_path)
        create_alignments_from_blast_record(seq, blast_record, enzyme_type)
        print(f"{seq.enzyme_name} blasted against all {enzyme_type}s")



if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    Alignment.drop_collection()

    enzyme_type = 'AAD'
    make_blast_db_for_enzyme_type(enzyme_type)
    do_all_by_all_blast(enzyme_type)
