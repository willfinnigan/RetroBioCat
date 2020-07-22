

def colour_edges_by_change_in_complexity(graph, edges):
    def change_to_colour(change, min_change=-0.5, max_change=0.5):
        normalised_change = 2*((change-min_change)/(max_change-min_change))-1

        if normalised_change > 1:
            normalised_change = 1
        elif normalised_change < -1:
            normalised_change = -1

        if normalised_change > 0:
            colour = (255-(255*normalised_change), 255, 0)
        elif normalised_change < 0:
            colour = (255, 255-(255*-normalised_change), 0)
        else:
            colour = (255,255,0)

        colour = 'rgb' + str(colour)
        return colour

    for i, edge in enumerate(edges):

        if graph.nodes[edge['from']]['attributes']['node_type'] == 'reaction':
            node = edge['from']
        else:
            node = edge['to']

        change = graph.nodes[node]['attributes']['change_in_complexity']
        colour = change_to_colour(change)
        edges[i]['color'] = colour

    return edges