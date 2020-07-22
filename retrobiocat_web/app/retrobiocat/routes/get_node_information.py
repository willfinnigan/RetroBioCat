import networkx as nx
import json
from retrobiocat_web.app.retrobiocat import bp
from flask import render_template, jsonify, request, current_app
from retrobiocat_web.app.retrobiocat.routes.network_explorer.functions import add_new
from retrobiocat_web.retro.evaluation import pubchem_funcs
from retrobiocat_web.retro.generation.network_generation.network import Network

from rdkit import Chem

def get_substrate_info(node):
    img = str(Chem.MolFromSmiles(node))
    img = img[:-2] + 'width="150" height="150" />'

    compound = pubchem_funcs.get_pubchem_compound_from_smiles(node)
    if compound != None:
        name = compound.iupac_name
        cid = compound.cid
        logp = compound.xlogp
    else:
        name = 'Not found'
        cid = ''
        logp = ''

    html = render_template('node_info_templates/substrate_node_info.html',
                           substrate_image=img,
                           smiles=node,
                           name=name,
                           cid=cid,
                           logp=logp)

    return {'name': node,
              'type': 'Substrate Node',
              'data': {},
              'html': html}

def get_reaction_info(node, network):
    def _get_rule_info(node, network):
        rule_info = {}
        rule_info['selected_enzyme'] = network.graph.nodes[node]['attributes']['selected_enzyme']
        rule_info['possible_enzymes'] = list(network.graph.nodes[node]['attributes']['possible_enzymes'])
        rule_info['name'] = network.graph.nodes[node]['attributes']['name']
        rule_info['img'] = ''
        return rule_info

    def _get_substrate_product_info(node, network):
        info = {}
        info['products'] = list(network.graph.successors(node))
        info['substrates'] = list(network.graph.predecessors(node))

        info['product_imgs'] = []
        for smi in info['products']:
            info['product_imgs'].append(str(Chem.MolFromSmiles(smi))[:-2] + 'width="150" height="150" />')
        info['substrate_imgs'] = []
        for smi in info['substrates']:
            info['substrate_imgs'].append(str(Chem.MolFromSmiles(smi))[:-2] + 'width="150" height="150" />')

        return info

    def _get_enzyme_index(selected_enzyme, possible_enzymes):
        if selected_enzyme not in possible_enzymes:
            index = 0
        else:
            index = possible_enzymes.index(selected_enzyme)
        return index

    html = render_template('node_info_templates/reaction_node_info.html',
                           rule_info=_get_rule_info(node, network),
                           substrate_product_info=_get_substrate_product_info(node, network))

    selected_enzyme = network.graph.nodes[node]['attributes']['selected_enzyme']
    possible_enzymes = list(network.graph.nodes[node]['attributes']['possible_enzymes'])
    enzyme_index = _get_enzyme_index(selected_enzyme, possible_enzymes)

    return {'name': node,
            'type': 'Reaction Node',
            'data': {'selected_enzyme': selected_enzyme,
                     'enzyme_select_index': enzyme_index,
                     'enzyme_options': possible_enzymes},
            'html': html}

def get_retrorule_info(node, network):
    def get_name(node, network):
        name = network.graph.nodes[node]['attributes']['name']
        return name

    def get_mnxr_name(node, network):
        name = network.graph.nodes[node]['attributes']['name']
        if name[0] == '1':
            mnxr = name[3:]
        else:
            mnxr = name[2:]
        return mnxr

    link = 'https://www.metanetx.org/equa_info/' + str(get_mnxr_name(node, network))

    html = render_template('node_info_templates/retrorule_reaction_node_info.html',
                           name=get_name(node, network),
                           mnxr_link=link)

    return {'name': node,
            'type': 'Reaction Node',
            'html': html}

@bp.route('/_node_info_reaction', methods=['GET', 'POST'])
def node_info_reaction():
    if 'node' in request.form:
        node = str(request.form['node'])
    else:
        result = {'name': 'error',
                  'type': 'error',
                  'data': {},
                  'html': ''}
        return jsonify(result=result)

    task_id = request.form['task_id']
    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)

    network = Network(graph=graph, target_smiles=target_smiles)
    network.settings = network_options
    network.add_attributes(attr_dict)
    network.get_node_types()

    if node in network.reaction_nodes:
        if 'retrorule' not in network.graph.nodes[node]['attributes']:
            result = get_reaction_info(node, network)
        else:
            result = get_retrorule_info(node, network)

    else:
        print('Error node not in reactions')
        result = {'name': 'error',
                  'type': 'error',
                  'data': {},
                  'html': ''}
    return jsonify(result=result)

@bp.route('/_node_info', methods=['GET', 'POST'])
def node_info():
    if 'node' in request.form:
        node = str(request.form['node'])
    else:
        result = {'name': 'error',
                  'type': 'error',
                  'data': {},
                  'html': ''}
        return jsonify(result=result)

    task_id = request.form['task_id']
    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)

    network = Network(graph=graph, target_smiles=target_smiles)
    network.settings = network_options
    network.add_attributes(attr_dict)
    network.get_node_types()

    if node in network.substrate_nodes:
        result = get_substrate_info(node)

    elif node in network.reaction_nodes:
        if 'retrorule' not in network.graph.nodes[node]['attributes']:
            result = get_reaction_info(node, network)
        else:
            result = get_retrorule_info(node, network)

    else:
        print('Error node not in substrates or reactions')
        result = {'name': 'error',
                  'type': 'error',
                  'data' : {},
                  'html' : ''}
    return jsonify(result=result)

@bp.route('/_change_enzyme', methods=['GET', 'POST'])
def change_enzyme():
    selected_node = request.form['selected_node']
    selected_enzyme = request.form['selected_enzyme']

    task_id = request.form['task_id']
    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.settings = network_options
    network.add_attributes(attr_dict)

    network.calculate_scores()

    network.graph.nodes[selected_node]['attributes']['selected_enzyme'] = selected_enzyme
    successors = list(network.graph.successors(selected_node))
    predecessors = list(network.graph.predecessors(selected_node))

    subgraph = network.graph.subgraph([selected_node]+successors+predecessors)
    nodes, edges = network.get_visjs_nodes_and_edges(graph=subgraph)

    data['attr_dict'] = json.dumps(network.attributes_dict())
    data['nodes'] = add_new(data['nodes'], nodes)
    data['edges'] = add_new(data['edges'], edges)

    result = {'nodes': nodes,
              'edges': edges}

    return jsonify(result=result)

