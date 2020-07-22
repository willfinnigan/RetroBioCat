
cofactors_to_count = ('NADPH, NADP+, NADH, NAD+, ATP, AMP, ADP')

def add_reaction_cofactor_count(graph, cofactorDict, include=cofactors_to_count):

    for node in list(graph):
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            if 'cofactors' not in graph.nodes[node]['attributes']:
                cofactors = {}
                name = graph.nodes[node]['attributes']['name']

                for enz in cofactorDict[name]:
                    cofactor_score = 0
                    if len(cofactorDict[name][enz]) != 0:
                        cofactor_plus = cofactorDict[name][enz][0]
                        cofactor_minus = cofactorDict[name][enz][1]

                        for cf in cofactor_plus:
                            if cf in include:
                                cofactor_score += 1

                        for cf in cofactor_minus:
                            if cf in include:
                                cofactor_score -= 1

                    cofactors[enz] = cofactor_score

                graph.nodes[node]['attributes']['cofactors'] = cofactors