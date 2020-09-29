from flask import render_template, jsonify, session, request, redirect, url_for
from retrobiocat_web.app.biocatdb.forms import SubstrateForm
from flask import current_app
from rq.job import Job
from retrobiocat_web.retro.enzyme_identification import query_mongodb
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.retro.enzyme_identification import molecular_similarity
from retrobiocat_web.app.biocatdb import bp
from rq import get_current_job
from retrobiocat_web.mongo.models.biocatdb_models import Activity, Paper
import numpy as np
from retrobiocat_web.app.biocatdb.functions.substrate_specificity import process_activity_data
import mongoengine as db
from retrobiocat_web.app.biocatdb.forms import SubstrateScopeForm

@bp.route('/substrate_specificity_form',  methods=['GET', 'POST'])
def substrate_specificity_form():
    print('Spec form')
    form = SubstrateForm()
    form.set_choices()
    task_id = ''

    if 'specificity_task_id' in session:
        old_task_id = session['specificity_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        form_data = form.data
        task = current_app.task_queue.enqueue(task_get_spec_data, form_data)
        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
                print('Deleted task - ' + str(old_task_id))
            except:
                print('No task to delete - ' + str(old_task_id))

        task_id = task.get_id()
        session['specificity_task_id'] = task_id

    return render_template('substrate_specificity/substrate_specificity_form_page.html', form=form, substrate_specificity_task_id=task_id)


@bp.route("/substrate_specificity_form_status/<task_id>", methods=["GET"])
def substrate_specificity_status(task_id):
    task = current_app.task_queue.fetch_job(task_id)
    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task.get_id(),
                "task_status": task.get_status(),
                "task_progress" : progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route("/substrate_specificity/<task_id>/", methods=["GET"])
def substrate_specificity(task_id):
    task = current_app.task_queue.fetch_job(task_id)
    activity_data = task.result

    return render_template('substrate_specificity/table_result_specificity.html', substrate_specificity_data=activity_data, title="Substrate specificity query")


def task_get_spec_data(form_data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()
    print('Started')

    enzyme_names = list(form_data['enzymes'].split(", "))
    reactions = list(form_data['reactions'].split(", "))

    if form_data['target_smiles'] != '':
        product = rdkit_smile(form_data['target_smiles'])
    else:
        product = form_data['target_smiles']

    similarity_cutoff = form_data['similarity']
    num_choices = form_data['num_choices']
    data_level = form_data['data_level']
    max_hits = form_data['max_hits']
    include_auto_data = bool(form_data['auto_data'])

    scorer = molecular_similarity.SubstrateSpecificityScorer(print_log=False)

    activity_df = scorer.querySpecificityDf(product, reactions, enzyme_names,
                                   dataLevel=data_level,
                                   numEnzymes=num_choices,
                                   simCutoff=similarity_cutoff,
                                   numHits=max_hits,
                                   include_auto_generated=include_auto_data)

    if activity_df is None:
        print('Activity df is none')
        return []

    if len(activity_df.index) == 0:
        print('Len activity df index is 0')
        return []

    activity_df = activity_df[process_activity_data.COLUMNS]
    activity_df = activity_df.round(2)
    activity_df.replace(np.nan, '', inplace=True)
    activity_df.replace(True, 'True', inplace=True)
    activity_df.replace(False, 'False', inplace=True)

    activity_data = activity_df.to_dict(orient='records')
    activity_data = process_activity_data.process_activity_data(activity_data)
    activity_data = process_activity_data.smiles_to_svg(activity_data)

    return activity_data

@bp.route("/paper_substrates/<paper_id>", methods=["GET"])
def paper_substrate_specificity(paper_id):

    paper = Paper.objects(id=paper_id)[0]

    activity_data = list(Activity.objects(paper=paper).only(*process_activity_data.mongo_cols).as_pymongo())
    activity_data = process_activity_data.process_activity_data(activity_data)
    activity_data = process_activity_data.smiles_to_svg(activity_data)

    return render_template('substrate_specificity/table_result_specificity.html', substrate_specificity_data=activity_data, title=f"{paper.short_citation} activity data")

@bp.route("/enzyme_substrates/<enzyme_name>", methods=["GET"])
def enzyme_substrate_specificity(enzyme_name):

    activity_data = list(Activity.objects(enzyme_name=enzyme_name).only(*process_activity_data.mongo_cols).as_pymongo())
    activity_data = process_activity_data.process_activity_data(activity_data)
    activity_data = process_activity_data.smiles_to_svg(activity_data)

    return render_template('substrate_specificity/table_result_specificity.html', substrate_specificity_data=activity_data, title=f"{enzyme_name} substrate specificity")

@bp.route("/enzyme_substrates_type/<enzyme_type>", methods=["GET"])
def enzyme_substrate_specificity_type(enzyme_type):
    if enzyme_type == 'All':
        enz_type_q = db.Q()
    else:
        enz_type_q = db.Q(enzyme_type=enzyme_type)

    activity_data = list(Activity.objects(enz_type_q).only(*process_activity_data.mongo_cols).as_pymongo())
    activity_data = process_activity_data.process_activity_data(activity_data)
    activity_data = process_activity_data.smiles_to_svg(activity_data)

    return render_template('substrate_specificity/table_result_specificity.html', substrate_specificity_data=activity_data,
                           title=f"Substrate scope for all {enzyme_type} enzymes")

@bp.route("/substrate_scope_search", methods=["GET", "POST"])
def substrate_scope_search():
    form = SubstrateScopeForm()
    form.set_choices()

    if form.validate_on_submit() == True:
        form_data = form.data

        if form_data['enzyme_name'] != 'All':
            return redirect(url_for("biocatdb.enzyme_substrate_specificity", enzyme_name=form_data['enzyme_name']))
        elif form_data['enzyme_type'] != 'All':
            return redirect(url_for("biocatdb.enzyme_substrate_specificity_type", enzyme_type=form_data['enzyme_type']))

    return render_template('substrate_specificity/substrate_scope_form.html', form=form)



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    scorer = molecular_similarity.SubstrateSpecificityScorer(print_log=False)

    test_product = ''
    reactions = []
    enzymes = ['CAR']

    activity_df = scorer.querySpecificityDf(test_product, reactions, enzymes)


    if len(activity_df.index) == 0:
        activity_data = []
    else:
        activity_df = activity_df[process_activity_data.COLUMNS]
        activity_df = activity_df.round(2)
        activity_df.replace(np.nan, '', inplace=True)
        activity_df.replace(True, 'True', inplace=True)
        activity_df.replace(False, 'False', inplace=True)

        activity_data = activity_df.to_dict(orient='records')


