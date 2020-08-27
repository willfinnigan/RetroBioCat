from rdkit.Chem import AllChem
from retrobiocat_web.retro.visualisation import rdkit_images

def add_substrate_info(nodes, graph):
    for i, node_dict in enumerate(nodes):
        node = nodes[i]['id']
        if graph.nodes[node]['attributes']['node_type'] == 'substrate':
            info = str(node) + '<br>'
            if 'complexity' in graph.nodes[node]['attributes']:
                info += 'Complexity: ' + str(round(graph.nodes[node]['attributes']['complexity'], 3)) + '<br>'
                info += 'Complexity relative to target: ' + str(
                    round(graph.nodes[node]['attributes']['relative_complexity'], 3)) + '<br>'
            if 'is_starting_material' in graph.nodes[node]['attributes']:
                info += 'Is Building Block: ' + str(
                    graph.nodes[node]['attributes']['is_starting_material']) + '<br>'

            nodes[i]['title'] = info

    return nodes

def add_enzyme_info(nodes, graph, img_size=(150,150)):
    def format_enzyme_info(info_dict, score):
        cols_to_ignore = ['smiles_reaction', 'paper_id', 'activity_id']

        if info_dict == False:
            return 'False'

        formatted = ''
        formatted += ('Score: ' + str(score) + '<br>')

        product_smiles = list(info_dict.keys())
        if len(product_smiles) >= 1:
            smi = product_smiles[0]
            formatted += ('<u>' + str(smi) + '</u>')
            formatted += '<br>'
            mol = AllChem.MolFromSmiles(smi)
            formatted += rdkit_images.moltosvg(mol, molSize=img_size)
            formatted += '<br>'

            for key, value in info_dict[smi].items():
                if key not in cols_to_ignore:
                    formatted += (str(key) + ': ' + str(value) + '<br>')
                #formatted += '<hr>'

        return formatted

    for i, node_dict in enumerate(nodes):
        node = node_dict['id']
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            enz = graph.nodes[node]['attributes']['selected_enzyme']
            if 'enzyme_info' in graph.nodes[node]['attributes']:
                info = graph.nodes[node]['attributes']['enzyme_info'][enz]
                score = graph.nodes[node]['attributes']['specificity_scores'][enz]
                if (info != False) and (score != 0):
                    formatted_info = format_enzyme_info(info, score)
                    if formatted_info != 'False':
                        nodes[i]['title'] = formatted_info

    return nodes

