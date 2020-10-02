import requests
from Bio.Blast import NCBIXML
from io import StringIO
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, UniRef90, EnzymeType, Alignment
import time
import mongoengine as db
import datetime
from flask import current_app

def get_sequence(identifier):
    url = f"https://www.uniprot.org/uniref/{identifier}.fasta"
    req = requests.get(url)

    if req.status_code in [200]:
        fasta = req.text
        seq = fasta[fasta.find('\n'):]
        seq = seq.replace('\n', '')
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

class BlastRunner(object):

    def __init__(self):
        self.polltime = 30
        self.run_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        self.status_url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status'
        self.result_url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result'

        self.exp = '1e-4'
        self.alignments = '1000'
        self.db = 'uniref90'
        self.email = 'wjafinnigan@gmail.com'

    def run(self, sequence):
        job_id = self._start_job(sequence)
        self._poll_till_complete(job_id)
        xml = self._get_blast_record(job_id)

        return xml

    def _start_job(self, sequence):
        data = {'email': self.email,
                 'program': 'blastp',
                 'stype': 'protein',
                 'sequence': sequence,
                 'database': self.db,
                 'exp': self.exp,
                 'alignments': self.alignments}

        response = requests.post(self.run_url, data=data)

        if response.status_code == 200:
            job_id = response.content.decode('UTF-8')
            print(f'BLAST job started - {job_id}')
            return job_id

        else:
            print(f"Error starting job - response code {response.status_code}")
            return False

    def _poll_till_complete(self, job_id):
        if job_id is False:
            return False

        status = False
        print(f"Waiting for job completion - {job_id}")

        while status != "FINISHED":
            time.sleep(self.polltime)
            url = f'{self.status_url}/{job_id}'
            response = requests.get(url)

            if response.status_code == 200:
                status = response.content.decode('UTF-8')
            else:
                print(f"Error polling job status - response code {response.status_code}")

        print(f"BLAST job complete - {job_id}")
        return True

    def _get_blast_record(self, job_id):
        url = f"{self.result_url}/{job_id}/out?format=5"
        response = requests.get(url)

        print(f"Parsing job - {job_id}")

        if response.status_code == 200:
            xml = response.content.decode('UTF-8')
            return xml

        else:
            print(response.status_code)
            return False

class BlastParser(object):

    def __init__(self):
        self.min_coverage = 0.8
        self.min_identity = 0.3

        self.min_size_frac = 0.8
        self.max_size_frac = 1.2
        self.max_e = 0.0005
        self.identifier_head = 'UR90:'
        self.blast_round = 1

    def parse(self, xml, seq_obj):
        blast_record = NCBIXML.read(StringIO(xml))

        query_length = len(seq_obj.sequence)
        enzyme_type_obj = EnzymeType.objects(enzyme_type=seq_obj.enzyme_type)[0]

        for alignment in blast_record.alignments:
            identifier = alignment.hit_id.replace(self.identifier_head, '')

            if self._alignment_filters(alignment, query_length):
                db_query = UniRef90.objects(db.Q(enzyme_name=identifier) & db.Q(enzyme_type=enzyme_type_obj))
                if len(db_query) == 0:
                    protein_sequence = get_sequence(identifier)
                    if self._sequence_filters(protein_sequence, query_length) :
                        if self._check_biocatdb_match(protein_sequence, identifier) is False:
                            print(f"Adding sequence for {identifier}")
                            self._add_uniref(alignment, identifier, protein_sequence, enzyme_type_obj, seq_obj)

    def _add_uniref(self, alignment, identifier, sequence, enzyme_type_obj, seq_seed):

        try:
            uniref_seq = UniRef90(enzyme_name=identifier,
                                  protein_name=get_name_from_header(alignment.title),
                                  tax=get_tax_from_header(alignment.title),
                                  tax_id=get_tax_id_from_header(alignment.title),
                                  sequence=sequence,
                                  enzyme_type=enzyme_type_obj,
                                  result_of_blasts_for=[seq_seed],
                                  blast_round=self.blast_round)
            uniref_seq.save()
        except Exception as e:
            print(e)

    def _alignment_filters(self, alignment, query_length):
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]
            coverage = hsp.align_length / query_length
            identity = hsp.identities / hsp.align_length
            if (coverage >= self.min_coverage) and (identity >= self.min_identity) and (hsp.expect <= self.max_e):
                return True
        return False

    def _sequence_filters(self, sequence, query_length):
        if sequence is not None:
            min_size = query_length * self.min_size_frac
            max_size = query_length * self.max_size_frac
            if (len(sequence) >= min_size) and (len(sequence) <= max_size):
                return True
        return False

    def _check_biocatdb_match(self, sequence, identifier):
        seq_query = Sequence.objects(sequence=sequence)
        if len(seq_query) != 0:
            biocatdb_seq = seq_query[0]
            biocatdb_seq.other_identifiers.append(identifier)
            biocatdb_seq.save()
            return True
        return False


def set_up_blast_job(enzyme_name, enzyme_type):
    current_app.app_context().push()

    seq = Sequence.objects(db.Q(enzyme_name=enzyme_name))[0]
    if seq.blast is None:
        print(f'Starting blast for sequence: {seq.enzyme_name}')
        xml = BlastRunner().run(seq.sequence)
        current_app.process_blasts_queue.enqueue(parse_blast_results, enzyme_name, xml,
                                                job_id=f'{enzyme_type}_{enzyme_name}_process_blast')

def parse_blast_results(enzyme_name, xml):
    seq = Sequence.objects(db.Q(enzyme_name=enzyme_name))[0]

    if xml is not None:
        BlastParser().parse(xml, seq)

        seq.blast = datetime.datetime.now()
        seq.save()
        print(f'Finished blast of sequence {seq.enzyme_name}')


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    seq = Sequence.objects(enzyme_name='CgDAADH')[0]
    protein_seq = seq.sequence

    xml = BlastRunner().run(protein_seq)
    BlastParser().parse(xml, seq)


