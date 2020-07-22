
def mark_node(list_nodes, target_node, colour, borderWidth=2):
    for i, node in enumerate(list_nodes):
        if list_nodes[i]['id'] == target_node:
            if colour != '':
                list_nodes[i]['color'] = colour
            list_nodes[i]['borderWidth'] = borderWidth

    return list_nodes

def colour_substrate_nodes_by_building_block(nodes, graph):

    for node_dict in nodes:
        node = node_dict['id']
        if graph.nodes[node]['attributes']['node_type'] == 'substrate':
            if graph.nodes[node]['attributes']['is_starting_material']:
                nodes = mark_node(nodes, node, 'purple')
            else:
                nodes = mark_node(nodes, node, 'blue')

    return nodes

def colour_substrates_by_relative_complexity(nodes, graph, substrate_nodes, reaction_nodes, borderwidth=2):
    def min_max_relative_complexity(graph):
        min_complexity = 0
        max_complexity = 0
        for node in list(graph):
            if graph.nodes[node]['attributes']['node_type'] == 'substrate':
                relative_complexity = graph.nodes[node]['attributes']['relative_complexity']
                if relative_complexity > max_complexity:
                    max_complexity = relative_complexity
                if relative_complexity < min_complexity:
                    min_complexity = relative_complexity

        return (min_complexity, max_complexity)

    # This code is really hacky and should be refactored ideally.

    all_nodes = []
    old_nodes, edges = get_nodes_and_edges(graph)
    for node in substrate_nodes:
        added = False
        for node_dict in nodes:
            if node_dict['id'] == node:
                all_nodes.append(node_dict)
                added = True
        if added == False:
            for node_dict in old_nodes:
                if node_dict['id'] == node:
                    all_nodes.append(node_dict)

    # add reactions nodes recently added
    for node in reaction_nodes:
        for node_dict in nodes:
            if node_dict['id'] == node:
                all_nodes.append(node_dict)

    min_comp, max_comp = min_max_relative_complexity(graph)

    for i, node_dict in enumerate(all_nodes):
        colour = ''
        node = node_dict['id']
        if graph.nodes[node]['attributes']['node_type'] == 'substrate':
            relative_complexity = graph.nodes[node]['attributes']['relative_complexity']
            if float(relative_complexity) > 0.0:
                normalised = relative_complexity / max_comp
                colour = 'rgb' + str((255, 255 - (255 * normalised), 0))
            elif float(relative_complexity) < 0.0:
                normalised = relative_complexity / min_comp
                colour = 'rgb' + str((255 - (255 * normalised), 255, 0))
            elif float(relative_complexity) == 0.0:
                colour = 'rgb' + str((255, 255, 0))

            if colour != '':
                nodes[i]['color'] = colour
                nodes[i]['borderWidth'] = borderwidth

    return all_nodes

def colour_reactions_by_activity(nodes, graph, show_negative=False):
    for i, node_dict in enumerate(nodes):
        node = node_dict['id']
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            enz = graph.nodes[node]['attributes']['selected_enzyme']
            if 'specificity_scores' in graph.nodes[node]['attributes']:
                activity = float(graph.nodes[node]['attributes']['specificity_scores'][enz])

                if (activity == 0):
                    colour = (200, 200, 200)

                else:
                    if activity > 0:
                        colour = (255 - (255 * activity), 255, 0)

                    elif (activity < 0) and (show_negative != False):
                        colour = (255, 255 - (255 * -activity), 0)
                    else:
                        colour = (200, 200, 200)

                colour = 'rgb' + str(colour)
                nodes[i]['color'] = colour
    return nodes

def colour_reactions_by_change_in_complexity(nodes, graph):
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

    for i, node_dict in enumerate(nodes):
        node = node_dict['id']
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            change = graph.nodes[node]['attributes']['change_in_complexity']
            colour = change_to_colour(change)
            nodes[i]['color'] = colour

    return nodes