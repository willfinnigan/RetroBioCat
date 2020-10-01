from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Molecule, Tag, EnzymeType, UniRef90, Alignment
from retrobiocat_web.analysis import embl_restfull, all_by_all_blast
from rq.registry import StartedJobRegistry
import datetime
import mongoengine as db
from rq.job import Job
import time

def set_bioinformatics_status(enzyme_type, status):
    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    enz_type_obj.bioinformatics_status = status
    enz_type_obj.save()

def set_blast_jobs(enzyme_type):
    blast_jobs = []
    last_blast_job = None
    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(bioinformatics_ignore__ne=True))
    for seq in seqs:
        print(seq.enzyme_name)
        if seq.sequence != '' and seq.sequence is not None and seq.blast is None:
            if len(seq.sequence) > 50:
                name = str(seq.enzyme_name)
                last_blast_job = current_app.blast_queue.enqueue(embl_restfull.set_up_blast_job, name, enzyme_type,
                                                                 job_id=f'{enzyme_type}_{seq.enzyme_name}_run_blast')
                blast_jobs.append(last_blast_job.id)
                print(f'Queued blast for {seq.enzyme_name}')
            else:
                print(f'Not blasting {seq.enzyme_name}')
                seq.blast = datetime.datetime.now()
        else:
            seq.blast = datetime.datetime.now()
        seq.save()

    return blast_jobs, last_blast_job

@bp.route('/_find_homologs', methods=['GET', 'POST'])
@roles_required('admin')
def find_homologs():
    enzyme_type = request.form['enzyme_type']
    set_bioinformatics_status(enzyme_type, 'Blasts Queued')
    current_app.blast_queue.enqueue(set_bioinformatics_status, enzyme_type, 'Running Blasts')

    blast_jobs, last_blast_job = set_blast_jobs(enzyme_type)

    current_app.task_queue.enqueue(check_if_all_blasts_complete, enzyme_type, blast_jobs, depends_on=last_blast_job)

    result = {'status': 'success',
              'msg': f"Started job to blast all {enzyme_type}'s",
              'issues': []}

    flash(f"Started job to blast all {enzyme_type}'s", "success")

    return jsonify(result=result)


@bp.route('/_find_all_homologs', methods=['GET', 'POST'])
@roles_required('admin')
def find_allhomologs():

    enzyme_types = EnzymeType.objects().distinct('enzyme_type')

    for enzyme_type in enzyme_types:
        set_bioinformatics_status(enzyme_type, 'Blasts Queued')
        current_app.blast_queue.enqueue(set_bioinformatics_status, enzyme_type, 'Running Blasts')
        blast_jobs, last_blast_job = set_blast_jobs(enzyme_type)
        current_app.task_queue.enqueue(check_if_all_blasts_complete, enzyme_type, blast_jobs, depends_on=last_blast_job)

    result = {'status': 'success',
              'msg': f"Started job to blast all enzyme_type's",
              'issues': []}

    flash(f"Started job to blast all enzyme_type's", "success")

    return jsonify(result=result)


def check_if_all_blasts_complete(enzyme_type, job_list):
    current_app.app_context().push()

    time.sleep(60)
    active_process_jobs = list(StartedJobRegistry(queue=current_app.process_blasts_queue).get_job_ids())
    active_process_jobs.extend(current_app.process_blasts_queue.job_ids)
    for job_id in active_process_jobs:
        if f"{enzyme_type}_" in job_id:
            job_list.append(job_id)

    jobs = Job.fetch_many(job_list, connection=current_app.redis)

    jobs_still_pending = []
    for job in jobs:
        status = job.get_status()
        if status == 'queued' or status == 'started':
            if job.id not in jobs_still_pending:
                jobs_still_pending.append(job.id)

    if len(jobs_still_pending) == 0:
        set_bioinformatics_status(enzyme_type, 'Blasts Complete')

    else:
        current_app.task_queue.enqueue(check_if_all_blasts_complete, enzyme_type, jobs_still_pending)

@bp.route('/_reset_blast_status', methods=['GET', 'POST'])
@roles_required('admin')
def reset_blast_status():
    enzyme_type = request.form['enzyme_type']
    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(bioinformatics_ignore__ne=True))
    for seq in seqs:
        seq.blast = None
        seq.save()

    result = {'status': 'success',
              'msg': f"Reset blast status to none for {enzyme_type}'s",
              'issues': []}

    flash(f"Reset blast status to none for {enzyme_type}'s", "success")
    return jsonify(result=result)

@bp.route('/bioinformatics_admin_page', methods=['GET', 'POST'])
@roles_required('admin')
def bioinformatics_admin_page():
    enzyme_types = EnzymeType.objects().order_by('enzyme_type')

    enzyme_bioinformatics_status = {}
    for enz_type_obj in enzyme_types:
        enz_type = enz_type_obj.enzyme_type
        enzyme_bioinformatics_status[enz_type] = enz_type_obj.bioinformatics_status

    enzyme_numbers = {}
    for enz_type_obj in enzyme_types:
        enz_type = enz_type_obj.enzyme_type
        enzyme_numbers[enz_type] = {}
        enzyme_numbers[enz_type]['biocatdb'] = len(Sequence.objects(enzyme_type=enz_type))
        enzyme_numbers[enz_type]['uniref'] = len(UniRef90.objects(enzyme_type=enz_type_obj))
        enzyme_numbers[enz_type]['alignments'] = len(Alignment.objects(enzyme_type=enz_type_obj))

    enz_type_dict = {}
    for enz_type_obj in enzyme_types:
        enz_type = enz_type_obj.enzyme_type
        enz_type_dict[enz_type] = 0
        seqs = Sequence.objects(enzyme_type=enz_type)
        if len(seqs) != 0:
            for seq in seqs:
                if seq.blast is not None:
                    enz_type_dict[enz_type] += 1
            if enz_type_dict[enz_type] != 0:
                enz_type_dict[enz_type] = round((enz_type_dict[enz_type]/len(seqs))*100, 0)

    aligned_enz_types = {}
    for enz_type_obj in enzyme_types:
        enz_type = enz_type_obj.enzyme_type
        aligned_enz_types[enz_type] = 0
        seq_objs = list(Sequence.objects(enzyme_type=enz_type)) + list(UniRef90.objects(enzyme_type=enz_type_obj))
        if len(seq_objs) != 0:
            for sq in seq_objs:
                if sq.alignments_made is not None:
                    aligned_enz_types[enz_type] += 1
            if aligned_enz_types[enz_type] != 0:
                aligned_enz_types[enz_type] = round((aligned_enz_types[enz_type]/len(seq_objs))*100, 0)

    registry = StartedJobRegistry(queue=current_app.blast_queue)
    num_jobs = registry.count

    return render_template('bioinformatics/bioinformatics_admin.html',
                           blasted_enz_types=enz_type_dict,
                           aligned_enz_types=aligned_enz_types,
                           enzyme_bioinformatics_status=enzyme_bioinformatics_status,
                           num_jobs=num_jobs,
                           enzyme_numbers=enzyme_numbers)


@bp.route('/_clear_all_bioinformatics_data', methods=['GET', 'POST'])
@roles_required('admin')
def clear_all_bioinformatics_data():
    enzyme_types = EnzymeType.objects()
    seqs = Sequence.objects()

    for enz in enzyme_types:
        enz.bioinformatics_status = 'Idle'
        enz.save()

    for seq in seqs:
        seq.blast = None
        seq.alignments_made = None
        seq.save()

    Alignment.drop_collection()
    UniRef90.drop_collection()

    result = {'status': 'success',
              'msg': f"Done",
              'issues': []}

    flash(f"Done", "success")
    return jsonify(result=result)