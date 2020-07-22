from retrobiocat_web.retro.generation.node_analysis import rdkit_smile


def add_scscore(graph, substrate_nodes, sc_score_model):
    for node in list(graph.nodes):
        if 'complexity' not in graph.nodes[node]['attributes']:
            if node in substrate_nodes:
                (smile, score) = sc_score_model.get_score_from_smi(node)
                graph.nodes[node]['attributes']['complexity'] = round(score, 4)
            else:
                graph.nodes[node]['attributes']['complexity'] = 0

    return graph

def max_complexity(graph, list_nodes):
    complexity = 0
    for node in list_nodes:
        if graph.nodes[node]['attributes']['complexity'] > complexity:
            complexity = graph.nodes[node]['attributes']['complexity']
    return complexity

def add_change_in_complexity(graph):
    for node in list(graph):
        if 'change_in_complexity' not in graph.nodes[node]['attributes']:
            successors = list(graph.successors(node))
            predecessors = list(graph.predecessors(node))
            difference = max_complexity(graph, predecessors) - max_complexity(graph, successors)
            graph.nodes[node]['attributes']['change_in_complexity'] = difference

    return graph

def add_relative_complexity(graph, target_smile):
    target_smile = rdkit_smile(target_smile)
    tm_complexity = graph.nodes[target_smile]['attributes']['complexity']

    for node in list(graph):
        if graph.nodes[node]['attributes']['node_type'] == 'substrate':
            if 'relative_complexity' not in graph.nodes[node]['attributes']:
                complexity = graph.nodes[node]['attributes']['complexity']
                relative_complexity = complexity - tm_complexity
                graph.nodes[node]['attributes']['relative_complexity'] = relative_complexity

    return graph

def add_reaction_relative_complexity(graph, target_smile):
    tm_complexity = graph.nodes[target_smile]['attributes']['complexity']

    for node in list(graph):
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            if 'reaction_avg_relative_complexity' not in graph.nodes[node]['attributes']:

                # calculate scores for reaction successor nodes
                successors = list(graph.successors(node))
                avg_suc_complexity_score = 0
                max_suc_complexity_score = 0
                min_suc_complexity_score = 5
                for suc_node in successors:
                    avg_suc_complexity_score += graph.nodes[suc_node]['attributes']['complexity']
                    if graph.nodes[suc_node]['attributes']['complexity'] > max_suc_complexity_score:
                        max_suc_complexity_score = graph.nodes[suc_node]['attributes']['complexity']
                    if graph.nodes[suc_node]['attributes']['complexity'] < min_suc_complexity_score:
                        max_suc_complexity_score = graph.nodes[suc_node]['attributes']['complexity']

                avg_suc_complexity_score = avg_suc_complexity_score / len(successors)

                graph.nodes[node]['attributes']['reaction_avg_relative_complexity'] = avg_suc_complexity_score - tm_complexity
                graph.nodes[node]['attributes']['reaction_max_relative_complexity'] = max_suc_complexity_score - tm_complexity
                graph.nodes[node]['attributes']['reaction_min_relative_complexity'] = min_suc_complexity_score - tm_complexity

def calc_complexities(target_smile, graph, substrate_nodes, sc_score_model):

    add_scscore(graph, substrate_nodes, sc_score_model)
    add_change_in_complexity(graph)
    add_relative_complexity(graph, target_smile)
    add_reaction_relative_complexity(graph, target_smile)

    return graph