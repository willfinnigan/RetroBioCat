from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, SSN_record, UniRef50, Sequence
from retrobiocat_web.analysis import embl_restfull
from flask import current_app
from retrobiocat_web.analysis import ssn_tasks
from redis import Redis
from rq import Queue
from rq.registry import ScheduledJobRegistry
from datetime import datetime
from datetime import timedelta
from retrobiocat_web.analysis.retrieve_uniref_info import UniRef_Parser
import mongoengine as db
import random
import time

"""
Every hour, check all blasts and all ssn's to see what needs updating
If bioinformatics status is not 'Complete', run blasts
If ssn status is not 'Complete', run ssn
"""

# 1. Check blast status - run blast if queued
# 2. Check ssn status, run ssn if queues
# 3. Every 30 minutes, check a random uniref sequence for updates, if updated... check all.

def schedual_jobs(repeat_in=30):
    registry = ScheduledJobRegistry(queue=current_app.auto_jobs)
    if 'schedule_job' not in list(registry.get_job_ids()):
        for job_id in registry.get_job_ids():
            registry.remove(job_id)

        print('Setting repeat jobs..')
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in), task_check_uniref_has_blast_source)
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in+4), check_random_uniref)
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in+8), task_check_blast_status)
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in+12), task_check_ssn_status)
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in+16), schedual_jobs, job_id='schedule_job')

def task_check_blast_status():
    if len(current_app.blast_queue.jobs) + len(current_app.process_blasts_queue.jobs) + len(current_app.alignment_queue.jobs) == 0:
        print('Checking blast status')
        enzyme_types = EnzymeType.objects()
        for enz_type in enzyme_types:
            embl_restfull.check_blast_status(enz_type.enzyme_type)

    else:
        print(f"Length blast queue = {len(current_app.blast_queue.jobs)}")
        print(f"Length process blast queue = {len(current_app.process_blasts_queue.jobs)}")
        print(f"Length alignment queue = {len(current_app.alignment_queue.jobs)}")

def task_check_ssn_status():
    for enzyme_type in EnzymeType.objects():
        ssn_query = list(SSN_record.objects(enzyme_type=enzyme_type))
        if len(ssn_query) > 1:
            print(f'Warning - multiple ssn records for {enzyme_type} - deleting extras')
            for i in range(1, len(ssn_query)):
                ssn_query[i].delete()


    if len(current_app.blast_queue.jobs) + len(current_app.process_blasts_queue.jobs) + len(current_app.alignment_queue.jobs) == 0:
        print('Checking ssn status')
        ssn_records = SSN_record.objects().select_related()

        for ssn_r in ssn_records:
            if ssn_r.status != 'Complete' and ssn_r.enzyme_type.bioinformatics_status == 'Complete':
                if len(UniRef50.objects(enzyme_type=ssn_r.enzyme_type)) != 0:
                    enzyme_type = ssn_r.enzyme_type.enzyme_type
                    job_name = f"{enzyme_type}_expand_ssn"
                    current_app.alignment_queue.enqueue(ssn_tasks.task_expand_ssn, enzyme_type, job_id=job_name)
                    print(f'Queued SSN job for {enzyme_type}')

        for enz_type_obj in EnzymeType.objects():
            if enz_type_obj.bioinformatics_status == 'Complete':
                if enz_type_obj not in SSN_record.objects().distinct('enzyme_type'):
                    unirefs = UniRef50.objects(enzyme_type=enz_type_obj)
                    biocatdb_seqs = list(Sequence.objects(db.Q(enzyme_type=enz_type_obj.enzyme_type) & db.Q(bioinformatics_ignore__ne=True)))
                    biocatdb_seqs = [seq for seq in biocatdb_seqs if seq.sequence != '' and seq.sequence is not None]

                    if len(unirefs) + len(biocatdb_seqs) != 0:
                        print(f"No SSN for {enz_type_obj.enzyme_type}, but blasts are complete and sequences present..  creating SSN.")
                        job_name = f"{enz_type_obj.enzyme_type}_expand_ssn"
                        current_app.alignment_queue.enqueue(ssn_tasks.task_expand_ssn, enz_type_obj.enzyme_type, job_id=job_name)

    else:
        print(f"Length blast queue = {len(current_app.blast_queue.jobs)}")
        print(f"Length process blast queue = {len(current_app.process_blasts_queue.jobs)}")
        print(f"Length alignment queue = {len(current_app.alignment_queue.jobs)}")

def task_check_uniref_has_blast_source():
    print('Checking for unirefs with no blast source..')
    uniref_query = UniRef50.objects(result_of_blasts_for__size=0)
    for uniref in uniref_query:
        print(f"Deleting {uniref.enzyme_name}")
        uniref.delete()

def check_random_uniref(num_to_check=25):
    for enzyme_type in EnzymeType.objects():

        unirefs = UniRef50.objects(enzyme_type=enzyme_type)

        all_match = True
        if len(unirefs) != 0:
            for i in range(num_to_check):
                rand_uniref = random.choice(unirefs)
                name = rand_uniref.enzyme_name
                ref_parser = UniRef_Parser()
                ref_parser.load_xml(name)
                time.sleep(0.2)

                if ref_parser.check_id_match(name) == False:
                    all_match = False

            if all_match != True:
                print(f'Identified mismatches with online uniref entries..  full uniref check for {enzyme_type.enzyme_type}')
                full_uniref_check(enzyme_type)

    print(f'Uniref checks complete ')

def full_uniref_check(enzyme_type_obj):
    unirefs = UniRef50.objects(enzyme_type=enzyme_type_obj).select_related()
    if len(unirefs) != 0:
        for ur in unirefs:
            print(f'Checking {ur.enzyme_name}..')
            ref_parser = UniRef_Parser()
            ref_parser.load_xml(ur.enzyme_name)
            time.sleep(0.2)

            if ref_parser.check_id_match(ur.enzyme_name) == False:
                print(f"{ur.enzyme_name} doesnt match cluster id online, deleting..")
                for seq in ur.result_of_blasts_for:
                    seq.blast = None
                    seq.save()
                ur.delete()

    ssn_query = SSN_record.objects(enzyme_type=enzyme_type_obj)
    if len(ssn_query) != 0:
        ssn_record = SSN_record.objects(enzyme_type=enzyme_type_obj)[0]
        ssn_record.status = 'Queued for update'
        ssn_record.save()

    enzyme_type_obj.bioinformatics_status = 'Queued for update'
    enzyme_type_obj.save()

    print(f"Full UniRef50 update complete for {enzyme_type_obj.enzyme_type}")



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    check_random_uniref(num_to_check=20)
