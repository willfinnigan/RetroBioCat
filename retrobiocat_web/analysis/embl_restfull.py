import requests
from Bio.Blast import NCBIXML
from io import StringIO
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, UniRef50, EnzymeType, SSN_record
import time
import mongoengine as db
import datetime
from flask import current_app


class BlastRunner(object):

    def __init__(self):
        self.polltime = 30
        self.run_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        self.status_url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status'
        self.result_url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result'

        self.alignments = '1000'
        self.db = 'uniref50'
        self.email = 'wjafinnigan@gmail.com'

    def run(self, sequence):
        job_id = self._start_job(sequence)
        self._poll_till_complete(job_id)
        xml = self._get_blast_record(job_id)
        if xml is None:
            return None

        blast_record = NCBIXML.read(StringIO(xml))

        return blast_record

    def _start_job(self, sequence):
        data = {'email': self.email,
                 'program': 'blastp',
                 'stype': 'protein',
                 'sequence': sequence,
                 'database': self.db,
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

    def __init__(self, log_level=0):
        self.min_coverage = 0.8
        self.min_identity = 0.3

        self.min_size_frac = 0.8
        self.max_size_frac = 1.2
        self.max_e = 5
        self.identifier_head = 'UR50:'
        self.blast_round = 1
        self.log_level = log_level

    def parse(self, output, seq_obj):
        blast_record = output
        query_length = len(seq_obj.sequence)
        enzyme_type_obj = EnzymeType.objects(enzyme_type=seq_obj.enzyme_type)[0]

        for alignment in blast_record.alignments:
            identifier = alignment.hit_id.replace(self.identifier_head, '')

            if self._alignment_filters(alignment, query_length):
                db_query = UniRef50.objects(db.Q(enzyme_name=identifier) & db.Q(enzyme_type=enzyme_type_obj)).select_related()
                if len(db_query) == 0:
                    protein_sequence = self._get_sequence(identifier)
                    if self._sequence_filters(protein_sequence, query_length):
                        self.log(f"Adding sequence for {identifier}")
                        self._add_uniref(alignment, identifier, protein_sequence, enzyme_type_obj, seq_obj)
                else:
                    uniref_obj = db_query[0]
                    self._add_result_of_blasts_for(seq_obj, uniref_obj)

    def _add_result_of_blasts_for(self, blasted_seq, uniref_obj):
        if blasted_seq not in uniref_obj.result_of_blasts_for:
            uniref_obj.result_of_blasts_for.append(blasted_seq)
            uniref_obj.save()
            self.log(f"Adding {blasted_seq.enzyme_name} as a blast source for for {uniref_obj.enzyme_name}")
        else:
            self.log(f"{blasted_seq.enzyme_name} is -already- a blast source for for {uniref_obj.enzyme_name}")

    def _add_uniref(self, alignment, identifier, sequence, enzyme_type_obj, seq_seed):

        try:
            uniref_seq = UniRef50(enzyme_name=identifier,
                                  protein_name=self._get_name_from_header(alignment.title),
                                  tax=self._get_tax_from_header(alignment.title),
                                  tax_id=self._get_tax_id_from_header(alignment.title),
                                  sequence=sequence,
                                  enzyme_type=enzyme_type_obj,
                                  result_of_blasts_for=[seq_seed],
                                  blast_round=self.blast_round)
            uniref_seq.save()
        except Exception as e:
            self.log(e)

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

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"BlastParser({level}): {msg}")

    @staticmethod
    def _get_sequence(identifier):
        url = f"https://www.uniprot.org/uniref/{identifier}.fasta"
        req = requests.get(url)

        if req.status_code in [200]:
            fasta = req.text
            seq = fasta[fasta.find('\n'):]
            seq = seq.replace('\n', '')
        else:
            seq = None

        return seq

    @staticmethod
    def _get_name_from_header(header):
        name_start = header.find(' ') + 1
        name_end = header.find(' n=')
        name = header[name_start:name_end]
        return name

    @staticmethod
    def _get_tax_from_header(header):
        tax_start = header.find('Tax=') + 4
        tax_end = header[tax_start:].find(' TaxID') + tax_start
        tax = header[tax_start:tax_end]
        return tax

    @staticmethod
    def _get_tax_id_from_header(header):
        tax_start = header.find('TaxID=') + 6
        tax_end = header[tax_start:].find(' ') + tax_start
        tax = header[tax_start:tax_end]
        return tax

def set_up_blast_job(enzyme_name):
    current_app.app_context().push()

    seq = Sequence.objects(db.Q(enzyme_name=enzyme_name))[0]
    if seq.blast is None:
        print(f'Starting blast for sequence: {seq.enzyme_name}')
        output = BlastRunner().run(seq.sequence)
        current_app.process_blasts_queue.enqueue(parse_blast_results, enzyme_name, output)

def parse_blast_results(enzyme_name, output):
    current_app.app_context().push()
    seq = Sequence.objects(db.Q(enzyme_name=enzyme_name))[0]

    if output is not None:
        BlastParser().parse(output, seq)

        seq.blast = datetime.datetime.now()
        seq.save()
        print(f'Finished blast of sequence {seq.enzyme_name}')

    current_app.task_queue.enqueue(check_blast_status, seq.enzyme_type)

def check_blast_status(enzyme_type):
    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(bioinformatics_ignore__ne=True) & db.Q(reviewed=True))
    all_complete = True
    for seq in seqs:
        if seq.blast is None:
            all_complete = False

    if all_complete == True:
        enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
        enz_type_obj.bioinformatics_status = 'Complete'
        enz_type_obj.save()

        ssn_q = SSN_record.objects(enzyme_type=enz_type_obj)
        if len(ssn_q) == 1:
            ssn_record = SSN_record.objects(enzyme_type=enz_type_obj)[0]
            ssn_record.status = 'Queued for update'
            ssn_record.save()
        else:
            print(f'Warning - multiple SSN records for {enz_type_obj.enzyme_type}, extras need deleting')

    else:
        set_blast_jobs(enzyme_type)

def set_bioinformatics_status(enzyme_type, status):
    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    enz_type_obj.bioinformatics_status = status
    enz_type_obj.save()

def set_blast_jobs(enzyme_type):
    set_bioinformatics_status(enzyme_type, 'Blasts Queued')
    current_app.blast_queue.enqueue(set_bioinformatics_status, enzyme_type, 'Running Blasts')

    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(bioinformatics_ignore__ne=True) & db.Q(reviewed=True))
    for seq in seqs:
        if seq.sequence != '' and seq.sequence is not None and seq.blast is None:
            if len(seq.sequence) > 50:
                name = str(seq.enzyme_name)
                current_app.blast_queue.enqueue(set_up_blast_job, name)
                print(f'Queued blast for {seq.enzyme_name}')
            else:
                print(f'Not blasting {seq.enzyme_name}')
                seq.blast = datetime.datetime.now()
        else:
            seq.blast = datetime.datetime.now()
        seq.save()

    current_app.task_queue.enqueue(check_blast_status, enzyme_type)



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    seq = Sequence.objects(enzyme_type='AAD')[1]
    protein_seq = seq.sequence

    xml = BlastRunner().run(protein_seq)
    BlastParser().parse(xml, seq)


