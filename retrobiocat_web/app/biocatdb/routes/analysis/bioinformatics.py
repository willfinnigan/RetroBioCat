from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Molecule, Tag, EnzymeType
from retrobiocat_web.analysis import embl_restfull
from rq.registry import StartedJobRegistry
import datetime
import mongoengine as db


def task_find_homologs(enzyme_name):
    seq = Sequence.objects(enzyme_name=enzyme_name)[0]
    print(f'Starting blast for sequence: {seq.enzyme_name}')

    try:
        job_id = embl_restfull.start_blast_job(seq.sequence)
        embl_restfull.poll_till_complete(job_id)
        blast_record = embl_restfull.get_blast_record(job_id)
        embl_restfull.parse_blast_record(blast_record, seq, blast_round=1)

        seq.blast = datetime.datetime.now()
        seq.save()
        print(f'Finished blast of sequence {seq.enzyme_name}')
    except Exception as e:
        print(e)

    return 'Done'


@bp.route('/_find_homologs', methods=['GET', 'POST'])
@roles_required('admin')
def find_homologs():
    enzyme_type = request.form['enzyme_type']

    seqs = Sequence.objects(db.Q(enzyme_type=enzyme_type) & db.Q(bioinformatics_ignore__ne=True))
    for seq in seqs:
        if seq.sequence != '' and seq.sequence is not None:
            if len(seq.sequence) > 50:
                seq.blast = None
                name = str(seq.enzyme_name)
                print(name)
                current_app.blast_queue.enqueue(task_find_homologs, name)
            else:
                seq.blast = datetime.datetime.now()
        else:
            seq.blast = datetime.datetime.now()
        seq.save()

    result = {'status': 'success',
              'msg': f"Started job to blast all {enzyme_type}'s",
              'issues': []}

    flash(f"Started job to blast all {enzyme_type}'s", "success")

    return jsonify(result=result)

@bp.route('/bioinformatics_admin_page', methods=['GET', 'POST'])
@roles_required('admin')
def bioinformatics_admin_page():
    enzyme_types = EnzymeType.objects.distinct('enzyme_type')

    enz_type_dict = {}
    for enz_type in enzyme_types:
        enz_type_dict[enz_type] = 0
        seqs = Sequence.objects(enzyme_type=enz_type)
        if len(seqs) != 0:
            for seq in seqs:
                if seq.blast is not None:
                    enz_type_dict[enz_type] += 1
            if enz_type_dict[enz_type] != 0:
                enz_type_dict[enz_type] = round((enz_type_dict[enz_type]/len(seqs))*100, 0)

    registry = StartedJobRegistry(queue=current_app.blast_queue)
    num_jobs = registry.count

    return render_template('bioinformatics/bioinformatics_admin.html',
                           blasted_enz_types=enz_type_dict,
                           num_jobs=num_jobs)
