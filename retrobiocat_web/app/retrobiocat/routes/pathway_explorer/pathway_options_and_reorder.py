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


def task_reorder_pathways(form, data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    graph_dict = json.loads(data['graph_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])
    attr_dict = json.loads(data['attr_dict'])
    all_pathways_data = json.loads(data['pathways_data'])

    pathways_id = form['pathways_id']

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.update_settings(network_options)
    network.add_attributes(attr_dict)

    pathways = []
    for list_varients in all_pathways_data:
        for item in list_varients:
            pathway = Pathway(item['nodes'], network, calc_scores=False)
            pathway.scores.scores_from_dict(item['scores_dict'])
            pathways.append(pathway)
    pathways = group_pathways(pathways)

    pathway_evaluator = PathwayEvaluator(pathways)
    pathway_evaluator.weights = {'Normalised num enzyme steps': form['weight_num_enzymes'],
                                 'Normalised change in complexity': form['weight_complexity'],
                                 'starting_material': form['weight_starting'],
                                 'postitive_enzyme_steps_score': form['weight_known_enzymes'],
                                 'Normalised Cofactor Stoichiometry': 0}

    pathway_evaluator.diversity_weight = form['weight_diversity']

    pathway_evaluator.evaluate()

    pathways_data = []
    for index, row in pathway_evaluator.df.iterrows():
        varients_data = []
        pathway_varients = [row['Pathway']]
        for pathway in row['Pathway'].other_varients:
            pathway_varients.append(pathway)

        for pathway_var in pathway_varients:
            varients_data.append({'nodes': pathway_var.list_nodes,
                                  'scores_dict': pathway_var.scores.scores_dict()
                                  })
        pathways_data.append(varients_data)

    job = get_current_job()
    job.meta['progress'] = 'complete'
    job.save_meta()

    result = {'id' : pathways_id,
              'pathways_data': json.dumps(pathways_data)}

    return result

@bp.route("/_change_pathway_options", methods=["POST"])
def change_pathway_options():
    substrate_colours = request.form['substrate_colours']
    reaction_colours = request.form['reaction_colours']
    edge_colours = request.form['edge_colours']
    pathway_num = int(request.form['pathway_num'])
    varient_num = int(request.form['varient_num'])

    task_id = request.form['task_id']
    data = json.loads(current_app.redis.get(task_id))

    graph_dict = json.loads(data['graph_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])
    attr_dict = json.loads(data['attr_dict'])

    all_pathways_data = json.loads(data['pathways_data'])
    pathway_data = all_pathways_data[pathway_num][varient_num]

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.update_settings(network_options)
    network.update_settings({'colour_substrates':substrate_colours,
                             'colour_reactions':reaction_colours,
                             "colour_arrows":edge_colours})

    network.add_attributes(attr_dict)

    pathway = Pathway(pathway_data['nodes'], network, calc_scores=False)

    nodes, edges = pathway.get_visjs_nodes_and_edges()

    data['network_options'] = json.dumps(network.settings)
    current_app.redis.mset({task_id: json.dumps(data)})

    result = {'nodes': nodes,
              'edges': edges,
              'network_options' : json.dumps(network.settings)
              }

    return jsonify(result=result)

@bp.route("/_reorder_pathways_status/<task_id>", methods=["GET"])
def reorder_pathways_status(task_id):
    task = current_app.pathway_queue.fetch_job(task_id)

    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        if task.get_status() == 'finished':

            result = task.result
            id = result.pop('id')
            data = json.loads(current_app.redis.get(id))
            data.update(result)
            current_app.redis.mset({id: json.dumps(data)})

            response_object = {
                "status": "success",
                "data": {
                    "task_id": task.get_id(),
                    "task_status": task.get_status(),
                    "task_progress": progress
                }
            }
        else:
            response_object = {
                "status": "success",
                "data": {
                    "task_id": task.get_id(),
                    "task_status": task.get_status(),
                    "task_progress": progress
                }
            }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

#ajax call used by pathway_explorer
@bp.route('/_reorder_pathways', methods=['GET', 'POST'])
def reorder_pathways():

    form = {'weight_complexity' : json.loads(request.form['weight_complexity']),
            'weight_num_enzymes' : json.loads(request.form['weight_num_enzymes']),
            'weight_starting' : json.loads(request.form['weight_starting']),
            'weight_known_enzymes' : json.loads(request.form['weight_known_enzymes']),
            'weight_diversity' : json.loads(request.form['weight_diversity']),
            'pathways_id' : request.form['id']}

    data = json.loads(current_app.redis.get(request.form['id']))
    task = current_app.pathway_queue.enqueue(task_reorder_pathways, form, data)
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
