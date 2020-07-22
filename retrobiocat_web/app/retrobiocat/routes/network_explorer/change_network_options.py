from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request
import networkx as nx
import json
from flask import current_app
from retrobiocat_web.app.retrobiocat.routes.network_explorer.functions import add_new
from retrobiocat_web.retro.generation.network_generation.network import Network

@bp.route("/_change_network_options", methods=["POST"])
def change_network_options():
    substrate_colours = request.form['substrate_colours']
    reaction_colours = request.form['reaction_colours']
    edge_colours = request.form['edge_colours']
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

    network.settings['colour_substrates'] = substrate_colours
    network.settings['colour_reactions'] = reaction_colours
    network.settings["colour_arrows"] = edge_colours

    nodes, edges = network.get_visjs_nodes_and_edges()

    data['nodes'] = add_new(data['nodes'], nodes)
    data['edges'] = add_new(data['edges'], edges)
    data['network_options'] = json.dumps(network.settings)

    current_app.redis.mset({task_id: json.dumps(data)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    result = {'nodes': nodes,
              'edges': edges,
              }

    return jsonify(result=result)

@bp.route("/_change_reaction_options", methods=["POST"])
def change_reaction_options():
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
    network.update_settings({'max_reactions' : int(request.form['max_reactions'])})

    """
    if len(retrobiocat.retrorules_diameters) != 0:
        network.settings['rr_min_diameter'] = int(request.form['rr_min_diameter'])
        network.settings['rr_min_products'] = int(request.form['rr_min_products'])
        network.settings['rr_max_reactions'] =  int(request.form['rr_max_reactions'])
    """

    data['network_options'] = json.dumps(network.settings)
    current_app.redis.mset({task_id: json.dumps(data)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    result = {'network_options': json.dumps(network.settings)}

    return jsonify(result=result)