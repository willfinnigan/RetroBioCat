import networkx as nx
from rdkit.Chem import AllChem, rdmolops

def get_nodes_with_no_successors(graph):

    no_succ = []
    for node in list(graph.nodes()):
        if len(list(graph.successors(node))) == 0:
            no_succ.append(node)
    return no_succ

def get_reaction_nodes_of_list_substrates(graph, list_substrates):
    reaction_nodes = []
    for substrate in list_substrates:
        predecessors = list(graph.predecessors(substrate))
        reaction_nodes.extend(predecessors)
    return reaction_nodes

def get_substrate_nodes(graph):
    substrate_nodes = []
    for node in list(graph):
        if 'node_type' not in graph.nodes[node]['attributes']:
            print('WARNING node_type not available for node: ' + str(node))
        if graph.nodes[node]['attributes']['node_type'] == 'substrate':
            substrate_nodes.append(node)
    return substrate_nodes

def get_reaction_nodes(graph):
    reaction_nodes = []
    for node in list(graph):
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            reaction_nodes.append(node)
    return reaction_nodes

def check_substrates_nx_predecessor(graph, smiles_to_check, origin_smile):
    any_are_predecessor = False

    # Get all the substrate nodes which are predecessors of origin_smile
    substrate_predecessors = []
    reaction_predecessors = list(graph.predecessors(origin_smile))
    for reaction_node in reaction_predecessors:
        substrate_predecessors.extend(list(graph.predecessors(reaction_node)))

    # do any of these match the smiles we want to check?
    for smiles in smiles_to_check:
        if smiles in substrate_predecessors:
            any_are_predecessor = True

    return any_are_predecessor

def substrate_precursors(graph, smile):
    substrate_precursors = []
    reaction_precursors = list(graph.successors(smile))

    for reaction in reaction_precursors:
        substrate_precursors.extend(list(graph.successors(reaction)))

    return substrate_precursors

def get_num_reactions_from_target(graph, smile, target):
    """ Returns the number of reaction nodes in shortest path from target to source"""
    num_in_path = nx.shortest_path_length(graph, source=target, target=smile)
    num_reacions_in_path = num_in_path/2
    return num_reacions_in_path

def get_current_pathway_length(graph, endNodes, target):
    # if the length of the pathway is max size
    currentLength = 0
    for endNode in endNodes:
        length = get_num_reactions_from_target(graph, endNode, target)
        if length > currentLength:
            currentLength = length
    return currentLength

def sort_by_score(list_nodes, list_scores, reverse=False):
    scores_nodes = list(zip(list_scores, list_nodes))
    scores_nodes.sort(reverse=reverse)
    nodes_sorted = [node for score, node in scores_nodes]
    return nodes_sorted

def get_reaction_names(reactions_list, graph):
    """ Returns a list of reaction names (minus uuid) """
    reaction_names = []
    for reaction in reactions_list:
        name = graph.nodes[reaction]['attributes']['name']
        reaction_names.append(name)
    return reaction_names

def rdkit_smile(smile, warning=False, inchi=False):
    try:
        mol = AllChem.MolFromSmiles(smile)
        if inchi == True:
            inchi_mol = AllChem.inchi.MolToInchi(mol)
            mol = AllChem.inchi.MolFromInchi(inchi_mol)
        smile = AllChem.MolToSmiles(mol)
    except:
        if warning == True:
            print('Warning couldnt convert smiles to rdkit - ' + str(smile))
        return None
    return smile

def get_inchl(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    inchi = AllChem.inchi.MolToInchi(mol)
    return inchi

def split_smi(smi):
    mol = AllChem.MolFromSmiles(smi)
    splitMols = rdmolops.GetMolFrags(mol, asMols=True)
    split_list = []
    for mol in splitMols:
        p_smile = AllChem.MolToSmiles(mol)
        split_list.append(p_smile)
    return split_list

def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')