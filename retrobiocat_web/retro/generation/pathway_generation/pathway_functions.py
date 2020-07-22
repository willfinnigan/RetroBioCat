import networkx as nx


def get_pathways(start, end, graph, cutoff=None):
    pathways = list(nx.all_simple_paths(graph, start, end, cutoff=cutoff))
    return pathways

def shortest_pathways(start, end, graph):
    pathways = list(nx.all_shortest_paths(graph, start, end))
    return pathways

def get_all_pathways(start, end_nodes, graph):
    all_pathways = {}
    for i, node in enumerate(end_nodes):
        for j, pathway in enumerate(get_pathways(start, node, graph)):
            name = str(i) + '.' + str(j)
            all_pathways[name] = pathway
    return all_pathways

def get_substrates_in_pathway(pathway, substrate_nodes):
    substrates = []
    for node in pathway:
        if node in substrate_nodes:
            substrates.append(node)

    return substrates

def get_shortest_pathways(start, end_nodes, graph):
    all_pathways = {}
    for i, node in enumerate(end_nodes):
        for j, pathway in enumerate(shortest_pathways(start, node, graph)):
            name = str(i) + '.' + str(j)
            all_pathways[name] = pathway
    return all_pathways

def flatten_list(list_to_flatten):
    list_to_flatten = list(list_to_flatten)
    def flatten(l):
        return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]
    flattened = flatten(list_to_flatten)
    if type(flattened[0]) == list:
        flattened = flattened[0]

    return flattened

def get_pathway_sub_graph(graph, pathway):
    pathway = list(pathway)
    flat_pathway = flatten_list(pathway)
    if len(pathway) == 0:
        flat_pathway = [pathway[0]]

    sub_graph = graph.subgraph(flat_pathway)

    return sub_graph



