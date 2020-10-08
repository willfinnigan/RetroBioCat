from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Molecule, Tag, EnzymeType
from retrobiocat_web.analysis import embl_restfull, all_by_all_blast, make_ssn
from rq.registry import StartedJobRegistry
import datetime
import mongoengine as db
from rq.job import Job
from rq import get_current_job
from retrobiocat_web.analysis.make_ssn import SSN
from retrobiocat_web.app.biocatdb.forms import SSN_Form
import json


@bp.route('/ssn_page/<task_id>/', methods=['GET'])
def ssn_page(task_id):
    task = current_app.network_queue.fetch_job(task_id)
    result = task.result

    return render_template('ssn/ssn.html',
                           nodes=result['nodes'],
                           edges=result['edges'],
                           alignment_score=result['alignment_score'],
                           max_nodes=result['max_nodes'],
                           cluster_nodes=result['cluster_nodes'])

def task_get_ssn(enzyme_type, min_score, combine_mutants, only_biocatdb, max_nodes):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    ssn = SSN(enzyme_type, print_log=True)
    ssn.load(include_mutants=not combine_mutants, only_biocatdb=only_biocatdb)
    cluster_nodes = ssn.get_nodes_to_cluster_on(starting_score=300, step=-2, min_edges=6)

    nodes, edges = ssn.visualise(min_score=min_score)

    result = {'nodes': nodes,
              'edges': edges,
              'max_nodes': max_nodes,
              'cluster_nodes': cluster_nodes,
              'alignment_score': min_score}
    return result

@bp.route("/ssn_status/<task_id>", methods=["GET"])
def ssn_status(task_id):
    task = current_app.network_queue.fetch_job(task_id)
    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task.get_id(),
                "task_status": task.get_status(),
                "task_progress": progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route('/ssn_form', methods=['GET', 'POST'])
@roles_required('experimental')
def ssn_form():
    form = SSN_Form()
    form.set_choices()

    task_id = ''

    if 'ssn_task_id' in session:
        old_task_id = session['ssn_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        enzyme_type = form.data['enzyme_type']
        min_score = form.data['alignment_score']
        combine_mutants = form.data['combine_mutants']
        only_biocatdb = form.data['only_biocatdb']
        max_nodes = form.data['max_nodes']
        task = current_app.network_queue.enqueue(task_get_ssn,
                                                 enzyme_type, min_score, combine_mutants, only_biocatdb, max_nodes)

        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
            except:
                pass

        task_id = task.get_id()
        session['ssn_task_id'] = task_id

    return render_template('ssn/ssn_form.html', form=form, task_id=task_id)


