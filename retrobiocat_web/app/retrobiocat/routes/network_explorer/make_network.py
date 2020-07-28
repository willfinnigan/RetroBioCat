from retrobiocat_web.app.retrobiocat import bp, forms
from flask import render_template, jsonify, session
from rq.job import Job
from rq import get_current_job
import networkx as nx
import json
from flask import current_app
import uuid
from retrobiocat_web.retro.generation.network_generation.network import Network

@bp.route('/network_explorer_form', methods=['GET', 'POST'])
def network_explorer_form():
    form = forms.NetworkExploreForm()
    task_id = ''

    if 'network_task_id' in session:
        old_task_id = session['network_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        form_data = form.data
        task = current_app.network_queue.enqueue(task_make_network, form_data)
        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
            except:
                pass

        task_id = task.get_id()
        session['network_task_id'] = task_id

    return render_template('network_explorer_form/network_explorer_form.html', form=form, task_id=task_id)

def task_make_network(form_data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    network = Network(include_experimental=bool(form_data['include_experimental']),
                      include_two_step=bool(form_data['include_two_step']),
                      print_log=not current_app.config['PRODUCTION'])

    network.update_settings({"allow_backwards_steps": bool(form_data['allow_backwards']),
                             "remove_simple": bool(form_data['remove_small']),
                             "similarity_score_threshold": float(form_data['sub_thres']),
                             "combine_enantiomers" : bool(form_data['combine_enantiomers']),
                             "num_enzymes": 1,
                             "calculate_complexities": bool(form_data['calc_complexity']),
                             "calculate_substrate_specificity": bool(form_data['sub_sim']),
                             "max_nodes": int(form_data['max_initial_nodes'],),
                             "colour_substrates" : form_data['colour_substrates'],
                             "colour_reactions" : form_data['colour_reactions'],
                             "colour_arrows": form_data['colour_edges'],
                             "show_negative_enzymes" : form_data['show_neg_enz'],
                             "only_postitive_enzyme_data" : not form_data['show_neg_enz'],
                             "max_reactions" : form_data["max_reactions"]})

    if form_data["specificity_scoring_mode"] == 'Product + substrates (slower)':
        network.update_settings({'specificity_score_substrates' : True})

    network.generate(form_data['target_smiles'], form_data['number_steps'], calculate_scores=False)

    job.meta['progress'] = 'network_generated'
    job.save_meta()

    network.calculate_scores()

    job.meta['progress'] = 'scores_calculated'
    job.save_meta()

    nodes, edges = network.get_visjs_nodes_and_edges()

    #options = {'interaction': {'multiselect': 'true',}}
    options = {}
    default_network_name = 'Network for ' + str(network.target_smiles)

    result = {'save_id':str(uuid.uuid4()),
              'save_links' : [],
              'save_name' : default_network_name,
              'nodes':nodes,
              'edges':edges,
              'options':json.dumps(options),
              'graph_dict':json.dumps(nx.to_dict_of_lists(network.graph)),
              'target_smiles':str(network.target_smiles),
              'network_options':json.dumps(network.settings),
              'attr_dict':json.dumps(network.attributes_dict()),
              'max_reactions' : int(network.settings['max_reactions'])}

    current_app.redis.mset({job.id: json.dumps(result)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(job.id, time_to_expire)

    return result

@bp.route("/network_explorer_status/<task_id>", methods=["GET"])
def network_explorer_status(task_id):
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
                "task_progress" : progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

