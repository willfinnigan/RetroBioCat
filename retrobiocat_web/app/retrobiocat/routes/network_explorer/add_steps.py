from retrobiocat_web.app.retrobiocat import bp, forms
from flask import render_template, jsonify, session, request, make_response
from rq.job import Job, get_current_job
import networkx as nx
import json
from flask import current_app
from retrobiocat_web.app.retrobiocat.routes.network_explorer.functions import add_new, delete_nodes_and_edges
from retrobiocat_web.retro.generation.network_generation.network import Network

#ajax call used by network_explorer
@bp.route('/_delete_step', methods=['GET', 'POST'])
def delete_step():
    reaction = request.form['reaction']
    task_id = request.form['task_id']

    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles, print_log=not current_app.config['PRODUCTION'])
    network.update_settings(network_options)
    network.add_attributes(attr_dict)

    to_delete = network.delete_reaction_node(reaction)
    nodes = []
    edges = []

    data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(data['nodes'], nodes)
    edges = add_new(data['edges'], edges)
    nodes, edges = delete_nodes_and_edges(to_delete, nodes, edges)
    data['nodes'] = nodes
    data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(data)})
    time_to_expire = 15 * 60  # 15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    result = {'to_delete': to_delete}
    return jsonify(result=result)


#ajax call used by network_explorer
@bp.route('/_step', methods=['GET', 'POST'])
def step():
    clicked_node = request.form['smiles']
    x = request.form['x']
    y = request.form['y']
    task_id = request.form['task_id']
    max_reactions = request.form['max_reactions']

    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles, print_log=not current_app.config['PRODUCTION'])
    network.update_settings(network_options)
    network.add_attributes(attr_dict)
    network.update_settings({'max_reactions': int(max_reactions)})

    new_substrate_nodes, new_reaction_nodes = network.add_step(clicked_node)

    all_new_nodes = [clicked_node] + new_substrate_nodes + new_reaction_nodes
    subgraph = network.graph.subgraph(all_new_nodes)

    nodes, edges = network.get_visjs_nodes_and_edges(graph=subgraph)

    for i, node in enumerate(nodes):
        nodes[i].update({'x': x,
                         'y': y})

    result = {'nodes': nodes,
              'edges': edges}

    data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(data['nodes'], nodes)
    edges = add_new(data['edges'], edges)
    nodes, edges = delete_nodes_and_edges([], nodes, edges)
    data['nodes'] = nodes
    data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(data)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    return jsonify(result=result)

#ajax call used by network_explorer
@bp.route('/_aizynth_step', methods=['GET', 'POST'])
def aizynth_step():
    clicked_node = request.form['smiles']
    x = request.form['x']
    y = request.form['y']
    task_id = request.form['task_id']
    max_reactions = request.form['max_reactions']

    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles, print_log=not current_app.config['PRODUCTION'])
    network.update_settings(network_options)
    network.add_attributes(attr_dict)
    network.update_settings({'max_reactions': int(max_reactions)})

    new_substrate_nodes, new_reaction_nodes = network.add_chemical_step(clicked_node)

    all_new_nodes = [clicked_node] + new_substrate_nodes + new_reaction_nodes
    subgraph = network.graph.subgraph(all_new_nodes)

    nodes, edges = network.get_visjs_nodes_and_edges(graph=subgraph)

    for i, node in enumerate(nodes):
        nodes[i].update({'x': x,
                         'y': y})

    result = {'nodes': nodes,
              'edges': edges}

    data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(data['nodes'], nodes)
    edges = add_new(data['edges'], edges)
    nodes, edges = delete_nodes_and_edges([], nodes, edges)
    data['nodes'] = nodes
    data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(data)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    return jsonify(result=result)

def task_add_retrorule_step(form_data, network_id):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    clicked_node = form_data['smiles']
    x = form_data['x']
    y = form_data['y']

    data = json.loads(current_app.redis.get(network_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])


    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.update_settings(network_options)
    network.add_attributes(attr_dict)


    network.retrorules.retrorules_rxns = current_app.retrorules_rxns
    network.retrorules.retrorule_db = current_app.retrorules_db

    new_substrate_nodes, new_reaction_nodes = network.retrorules.add_step(clicked_node)

    all_new_nodes = [clicked_node] + new_substrate_nodes + new_reaction_nodes
    subgraph = network.graph.subgraph(all_new_nodes)

    nodes, edges = network.get_visjs_nodes_and_edges(graph=subgraph)

    for i, node in enumerate(nodes):
        nodes[i].update({'x': x,
                         'y': y})

    result = {'nodes': nodes,
              'edges': edges,
              'to_delete': [],
              }

    data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    data['attr_dict'] = json.dumps(network.attributes_dict())
    data['nodes'] = add_new(data['nodes'], nodes)
    data['edges'] = add_new(data['edges'], edges)

    current_app.redis.mset({network_id: json.dumps(data)})
    current_app.redis.expire(network_id, 5*60)

    return result

@bp.route("/_retrorules_step_status/<task_id>", methods=["GET"])
def retrorules_step_status(task_id):
    task = current_app.retrorules_queue.fetch_job(task_id)
    progress = 'queuing'

    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        if task.get_status() == 'finished':

            response_object = {
                "status": "success",
                "data": {
                    "task_id": task.get_id(),
                    "task_status": task.get_status(),
                    "task_result" : task.result,
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

#ajax call used by main_site to add custom reaction
@bp.route('/_custom_reaction', methods=['GET', 'POST'])
def custom_reaction():
    product_smiles = str(request.form['product'])
    substrate_smiles = str(request.form['substrate'])
    reaction_name = str(request.form['name'])
    task_id = request.form['task_id']

    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.update_settings(network_options)
    network.add_attributes(attr_dict)

    new_substrate_nodes, new_reaction_nodes = network.custom_reaction(product_smiles, substrate_smiles, reaction_name)

    all_new_nodes = new_substrate_nodes + new_reaction_nodes
    subgraph = network.graph.subgraph(all_new_nodes)
    nodes, edges = network.get_visjs_nodes_and_edges(graph=subgraph)

    result = {'nodes' : nodes,
              'edges' : edges,
              }

    data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(data['nodes'], nodes)
    edges = add_new(data['edges'], edges)
    data['nodes'] = nodes
    data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(data)})
    current_app.redis.expire(task_id, 5*60)

    return jsonify(result=result)

