from flask import render_template, jsonify, request, current_app, session
from retrobiocat_web.app.retrobiocat import bp, forms
from rq import get_current_job
import networkx as nx
import json
from flask import current_app
from rq.job import Job
from retrobiocat_web.retro.generation.network_generation.network import Network
from retrobiocat_web.retro.generation.pathway_generation.pathway import Pathway
from retrobiocat_web.retro.generation.pathway_generation.pathway_scoring import PathwayEvaluator
from retrobiocat_web.retro.generation.pathway_generation.group_pathways import group_pathways
from retrobiocat_web.app.retrobiocat.routes.pathway_explorer.pathway import evaluate_pathways, package_evaluated_pathways, package_visjs_pathways
import datetime

def load_pathways(all_pathways_nodes, all_pathway_scores, network):
    pathways = []
    for i, list_nodes in enumerate(all_pathways_nodes):
        pathway = Pathway(list_nodes, network, calc_scores=False)
        pathway.scores.scores_from_dict(all_pathway_scores[i])
        pathways.append(pathway)
    return pathways

def task_reorder_pathways(weights, pathways_id):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    pathway_settings = json.loads(current_app.redis.get(pathways_id + '__pathway_settings'))
    pathway_settings.update({'weight_num_enzymes': weights[0],
                             'weight_complexity': weights[1],
                             'weight_starting': weights[2],
                             'weight_known_enzymes': weights[3],
                             'weight_diversity': weights[4]})
    current_app.redis.mset({f"{pathways_id}__pathway_settings": json.dumps(pathway_settings)})
    current_app.redis.expire(pathways_id, 60 * 60)

    network_data = json.loads(current_app.redis.get(pathways_id + '__network'))
    graph = nx.from_dict_of_lists(json.loads(network_data['graph_dict']), create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=network_data['target_smiles'],
                      print_log=not current_app.config['PRODUCTION'])
    network.update_settings(json.loads(network_data['network_options']))
    network.add_attributes(json.loads(network_data['attr_dict']))

    all_pathways_nodes, all_scores = json.loads(current_app.redis.get(f"{pathways_id}__all_pathways"))
    pathways = load_pathways(all_pathways_nodes, all_scores, network)

    pathway_evaluator = evaluate_pathways(pathways, weights)
    print(weights)
    package_evaluated_pathways(pathway_evaluator.df, pathways_id)
    package_visjs_pathways(pathways_id)

    job = get_current_job()
    job.meta['progress'] = 'complete'
    job.save_meta()

    result = {}

    return result


@bp.route("/_reorder_pathways_status/<task_id>", methods=["GET"])
def reorder_pathways_status(task_id):
    task = current_app.pathway_queue.fetch_job(task_id)
    task_id = task.get_id()
    task_status = task.get_status()
    seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
    print(seconds_since_active)
    if seconds_since_active > 180:
        print('Job no longer active')
        task_status = 'failed'

    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        if task.get_status() == 'finished':

            response_object = {
                "status": "success",
                "data": {
                    "task_id": task_id,
                    "task_status": task_status,
                    "task_progress": progress
                }
            }
        else:
            response_object = {
                "status": "success",
                "data": {
                    "task_id": task_id,
                    "task_status": task_status,
                    "task_progress": progress
                }
            }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

#ajax call used by pathway_explorer
@bp.route('/_reorder_pathways', methods=['GET', 'POST'])
def reorder_pathways():

    weights = [json.loads(request.form['weight_num_enzymes']),
               json.loads(request.form['weight_complexity']),
               json.loads(request.form['weight_starting']),
               json.loads(request.form['weight_known_enzymes']),
               json.loads(request.form['weight_diversity'])]
    task = current_app.pathway_queue.enqueue(task_reorder_pathways, weights, request.form['id'])
    reorder_task_id = task.get_id()

    if 'pathway_reorder_task_id' in session:
        old_task_id = session['pathway_reorder_task_id']
        try:
            old_job = Job.fetch(old_task_id, connection=current_app.redis)
            old_job.delete()
        except:
            pass

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": reorder_task_id,
                "task_status": task.get_status(),
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202


@bp.route("/_change_pathway_options", methods=["POST"])
def change_pathway_options():
    reaction_colours = request.form['reaction_colours']
    edge_colours = request.form['edge_colours']
    pathways_id = request.form['task_id']
    pathway_num = int(request.form['pathway_num'])
    varient_num = int(request.form['varient_num'])

    task = current_app.pathway_queue.enqueue(change_options_job, reaction_colours, edge_colours, pathways_id, pathway_num, varient_num)
    options_task_id = task.get_id()

    if 'pathway_options_task_id' in session:
        old_task_id = session['pathway_options_task_id']
        try:
            old_job = Job.fetch(old_task_id, connection=current_app.redis)
            old_job.delete()
        except:
            pass

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": options_task_id,
                "task_status": task.get_status(),
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route("/_check_options_status/<task_id>", methods=["GET"])
def check_options_status(task_id):
    task = current_app.pathway_queue.fetch_job(task_id)

    if task:
        if task.get_status() == 'finished':

            response_object = {
                "status": "success",
                "data": {
                    "task_id": task.get_id(),
                    "task_status": task.get_status(),
                },
                'result': task.result,
            }
        else:
            response_object = {
                "status": "success",
                "data": {
                    "task_id": task.get_id(),
                    "task_status": task.get_status()
                }
            }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202


def change_options_job(reaction_colours, edge_colours, pathways_id, pathway_num, varient_num):
    network_data = json.loads(current_app.redis.get(pathways_id + '__network'))
    print(network_data)

    network_options = json.loads(network_data['network_options'])
    network_options.update({'colour_reactions': reaction_colours,
                            "colour_arrows": edge_colours})
    network_data['network_options'] = json.dumps(network_options)
    current_app.redis.mset({f"{pathways_id}__network": json.dumps(network_data)})
    current_app.redis.expire(f"{pathways_id}__network", 60 * 60)

    package_visjs_pathways(pathways_id, max_vis=100)

    pathway_data = json.loads(current_app.redis.get(f"{pathways_id}__{pathway_num}"))
    nodes, edges, max_varient = pathway_data[varient_num - 1]

    return {'nodes': nodes,
            'edges': edges}