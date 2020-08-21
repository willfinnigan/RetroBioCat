import requests
from Bio.Blast import NCBIXML
from io import StringIO
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, UniRef90, EnzymeType
import time
import mongoengine as db

def start_blast_job(sequence, exp='1e-4', alignments='1000', db='uniref90'):

    url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
    data = {'email': 'wjafinnigan@gmail.com',
               'program': 'blastp',
               'stype': 'protein',
               'sequence': sequence,
               'database': db,
               'exp': exp,
               'alignments': alignments}

    response = requests.post(url, data=data)

    if response.status_code == 200:
        job_id = response.content.decode('UTF-8')
        print(f'BLAST job started - {job_id}')
        return job_id

    else:
        print(f"Error starting job - response code {response.status_code}")
        return False

def poll_till_complete(job_id, poll_time=10):
    status = False
    print(f"Waiting for job completion - {job_id}")

    while status != "FINISHED":
        time.sleep(poll_time)
        url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}'
        response = requests.get(url)

        if response.status_code == 200:
            status = response.content.decode('UTF-8')
        else:
            print(f"Error polling job status - response code {response.status_code}")

    print(f"BLAST job complete - {job_id}")
    return True

def get_blast_record(job_id):
    url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/out?format=5"

    response = requests.get(url)

    print(f"Parsing job - {job_id}")

    if response.status_code == 200:
        xml = response.content.decode('UTF-8')
        blast_record = NCBIXML.read(StringIO(xml))
        return blast_record

    else:
        print(response.status_code)
        return False

def get_sequence(identifier):
    url = f"https://www.uniprot.org/uniref/{identifier}.fasta"
    req = requests.get(url)

    if req.status_code in [200]:
        fasta = req.text
        seq = fasta[fasta.find('\n'):]
    else:
        seq = None

    return seq

def get_name_from_header(header):
    name_start = header.find(' ') + 1
    name_end = header.find(' n=')
    name = header[name_start:name_end]
    return name

def get_tax_from_header(header):
    tax_start = header.find('Tax=') + 4
    tax_end = header[tax_start:].find(' TaxID') + tax_start
    tax = header[tax_start:tax_end]
    return tax

def get_tax_id_from_header(header):
    tax_start = header.find('TaxID=') + 6
    tax_end = header[tax_start:].find(' ') + tax_start
    tax = header[tax_start:tax_end]
    return tax


def parse_blast_record(blast_record, seq_seed,
                       min_coverage=0.8, min_identity=0.4,
                       min_size_frac=0.8, max_size_frac=1.2,
                       identifier_head='UR90:',
                       blast_round=1):
    print(f"Adding sequences to database from blast of {seq_seed.enzyme_name}")

    query_length = len(seq_seed.sequence)
    min_size = query_length*min_size_frac
    max_size = query_length*max_size_frac
    enzyme_type_obj = EnzymeType.objects(enzyme_type=seq_seed.enzyme_type)[0]

    for alignment in blast_record.alignments:
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]
            identifier = alignment.hit_id.replace(identifier_head,'')
            coverage = hsp.align_length / query_length
            identity = hsp.identities / hsp.align_length
            if (coverage >= min_coverage) and (identity >= min_identity):
                db_query = UniRef90.objects(db.Q(enzyme_name=identifier) & db.Q(enzyme_type=enzyme_type_obj))
                if len(db_query) == 0:
                    protein_sequence = get_sequence(identifier)
                    protein_sequence = protein_sequence.replace('\n', '')
                    if protein_sequence is not None:
                        if (len(protein_sequence) >= min_size) and (len(protein_sequence) <= max_size):

                            seq_query = Sequence.objects(sequence=protein_sequence)
                            if len(seq_query) != 0:
                                biocatdb_seq = seq_query[0]
                                biocatdb_seq.other_identifiers.append(identifier)
                                biocatdb_seq.save()

                            else:
                                print(f"Adding sequence for {alignment.title}")
                                print(f"Coverage = {coverage}, Identity = {identity}")
                                uniref_seq = UniRef90(enzyme_name=identifier,
                                                      protein_name=get_name_from_header(alignment.title),
                                                      tax=get_tax_from_header(alignment.title),
                                                      tax_id=get_tax_id_from_header(alignment.title),
                                                      sequence=protein_sequence,
                                                      enzyme_type=enzyme_type_obj,
                                                      result_of_blasts_for=[seq_seed],
                                                      blast_round=blast_round)
                                uniref_seq.save()

                else:
                    uniref_seq = db_query[0]
                    print(f"Already present - {uniref_seq.enzyme_name}")
                    if seq_seed not in uniref_seq.result_of_blasts_for:
                        uniref_seq.result_of_blasts_for.append(seq_seed)
                        uniref_seq.save()



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    UniRef90.drop_collection()

    seq = Sequence.objects(enzyme_name='mpCAR')[0]
    protein_seq = seq.sequence

    job_id = start_blast_job(protein_seq)
    poll_till_complete(job_id)
    blast_record = get_blast_record(job_id)
    parse_blast_record(blast_record, seq)


