

def add_cofactors(graph, edges):
    for i, edge in enumerate(edges):
        label = ''
        colour = ''

        if graph.nodes[edge['from']]['attributes']['node_type'] == 'reaction':
            node = edge['from']
            colour = 'red'
            selected_enzyme = graph.nodes[node]['attributes']['selected_enzyme']
            cofactors = graph.nodes[node]['attributes']['enzyme_cofactors']
            if selected_enzyme in cofactors:
                selected_cofactors = cofactors[selected_enzyme]
                for cf in selected_cofactors[0]:
                    if len(label) > 0:
                        label += '\n'
                    label += '- ' + str(cf)

        elif graph.nodes[edge['to']]['attributes']['node_type'] == 'reaction':
            # if more than one successor, only add to [0]
            cofactor_info = True
            sucs = list(graph.successors(edge['to']))
            if len(sucs) > 1:
                if edge['from'] != sucs[0]:
                    cofactor_info = False

            if cofactor_info == True:
                node = edge['to']
                colour = 'green'
                selected_enzyme = graph.nodes[node]['attributes']['selected_enzyme']
                cofactors = graph.nodes[node]['attributes']['enzyme_cofactors']
                if selected_enzyme in cofactors:
                    selected_cofactors = cofactors[selected_enzyme]
                    for cf in selected_cofactors[1]:
                        if len(label) > 0:
                            label += '\n'
                        label += '+ ' + str(cf)

        edges[i].update({'label': label,
                         'font': {'size': 6,
                                  'align': 'top',
                                  'multi': 'md',
                                  'color': colour
                                  }
                         })

    return edges