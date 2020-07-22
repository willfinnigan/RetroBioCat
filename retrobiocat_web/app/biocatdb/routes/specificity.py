from flask import render_template, jsonify, session
from retrobiocat_web.app.biocatdb.forms import SubstrateForm
from retrobiocat_web.app.biocatdb.routes.specificity_query_task import task_run_query
from flask import current_app
from rq.job import Job
from retrobiocat_web.retro.enzyme_identification import query_mongodb
from retrobiocat_web.app.biocatdb import bp


@bp.route('/substrate_specificity_form',  methods=['GET', 'POST'])
def substrate_specificity_form():
    print('Spec form')
    form = SubstrateForm()
    enzymes = query_mongodb.get_enzymes_in_db() + ['All']
    reactions = query_mongodb.get_reactions_in_db()
    task_id = ''

    if 'specificity_task_id' in session:
        old_task_id = session['specificity_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        form_data = form.data
        task = current_app.task_queue.enqueue(task_run_query, form_data)
        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
                print('Deleted task - ' + str(old_task_id))
            except:
                print('No task to delete - ' + str(old_task_id))

        task_id = task.get_id()
        session['specificity_task_id'] = task_id
        print(task_id)

    return render_template('substrate_specificity/substrate_specificity_form.html', form=form, enzymes=enzymes, reactions=reactions, task_id=task_id)


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
    result = task.result

    return render_template('substrate_specificity/substrate_specificity.html',
                           headings=result['headings'],
                           rows=result['rows'])