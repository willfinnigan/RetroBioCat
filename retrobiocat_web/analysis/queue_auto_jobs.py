from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, SSN_record, UniRef50
from retrobiocat_web.app.db_analysis.routes.bioinformatics import set_blast_jobs
from flask import current_app
from retrobiocat_web.analysis import ssn_tasks
from redis import Redis
from rq import Queue
from rq.registry import ScheduledJobRegistry
from datetime import datetime
from datetime import timedelta

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
    if registry.count < 3:
        for job_id in registry.get_job_ids():
            registry.remove(job_id)

        print('Setting repeat jobs..')
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in), task_check_blast_status)
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in+2), task_check_ssn_status)
        current_app.auto_jobs.enqueue_in(timedelta(minutes=repeat_in+4), schedual_jobs)

    print(registry.count)
    print(registry.get_job_ids())

def task_check_blast_status():
    if len(current_app.blast_queue.jobs) + len(current_app.process_blasts_queue.jobs) + len(current_app.alignment_queue.jobs) == 0:
        print('Checking blast status')
        enzyme_types = EnzymeType.objects()
        for enz_type in enzyme_types:
            if enz_type.bioinformatics_status != 'Complete':
                set_blast_jobs(enz_type.enzyme_type)
    else:
        print(f"Length blast queue = {len(current_app.blast_queue.jobs)}")
        print(f"Length process blast queue = {len(current_app.process_blasts_queue.jobs)}")
        print(f"Length alignment queue = {len(current_app.alignment_queue.jobs)}")

def task_check_ssn_status():
    if len(current_app.blast_queue.jobs) + len(current_app.process_blasts_queue.jobs) + len(current_app.alignment_queue.jobs) == 0:
        print('Checking ssn status')
        ssn_records = SSN_record.objects().select_related()

        for ssn_r in ssn_records:
            if ssn_r.status != 'Complete':
                if len(UniRef50.objects(enzyme_type=ssn_r.enzyme_type)) != 0:
                    enzyme_type = ssn_r.enzyme_type.enzyme_type
                    job_name = f"{enzyme_type}_expand_ssn"
                    current_app.alignment_queue.enqueue(ssn_tasks.task_expand_ssn, enzyme_type, job_id=job_name)
                    print(f'Queued SSN job for {enzyme_type}')

        for enz_type_obj in EnzymeType.objects():
            if enz_type_obj.bioinformatics_status == 'Complete':
                if enz_type_obj not in SSN_record.objects().distinct('enzyme_type'):
                    enzyme_type = enz_type_obj.enzyme_type
                    print(f"No SSN for {enzyme_type}, but blasts are complete..  creating SSN.")
                    job_name = f"{enzyme_type}_expand_ssn"
                    current_app.alignment_queue.enqueue(ssn_tasks.task_expand_ssn, enzyme_type, job_id=job_name)

    else:
        print(f"Length blast queue = {len(current_app.blast_queue.jobs)}")
        print(f"Length process blast queue = {len(current_app.process_blasts_queue.jobs)}")
        print(f"Length alignment queue = {len(current_app.alignment_queue.jobs)}")

def task_check_uniref_has_blast_source():
    uniref_query = UniRef50.objects(result_of_blasts_for__size=0)
    for uniref in uniref_query:
        uniref.delete()

def check_random_uniref():
    pass

def create_check_all_uniref_jobs():
    pass