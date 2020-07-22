from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.generation.pathway_generation import pathway_functions
from retrobiocat_web.retro.visualisation import visualise
from retrobiocat_web.retro.generation.pathway_generation import pathway_scoring

class Pathway(object):

    def __init__(self, list_nodes, network, calc_scores=True, attr_dict=None, reaction_enzyme_map=None):
        self.list_nodes = list_nodes
        self.network = network
        self.target_smiles = network.target_smiles
        self.sub_graph = pathway_functions.get_pathway_sub_graph(self.network.graph, list_nodes)

        # for loading pre-made pathway
        self.attributes = {}
        if attr_dict == None:
            self.get_attributes_from_network()
        else:
            self.set_attributes_from_dict(attr_dict)

        self.substrates = node_analysis.get_substrate_nodes(self.sub_graph)
        self.reactions = node_analysis.get_reaction_nodes(self.sub_graph)
        self.end_nodes = node_analysis.get_nodes_with_no_successors(self.sub_graph)
        self.reaction_names = self._get_reaction_names()

        # specifies which enzyme for each step
        if reaction_enzyme_map == None:
            self.set_reaction_enzyme_map_from_network()
        else:
            self.reaction_enzyme_map = reaction_enzyme_map

        # scores
        self.scores = pathway_scoring.PathwayScores(self, calc_scores=calc_scores)

        self.other_varients = []
        self.other_varients_as_nodes = []

    def __repr__(self):
        """ Returns the pathway nodes as a string - useful for dataframe"""
        return str(self.list_nodes)

    def get_attributes_from_network(self):
        """ Gets the attributes of each node from the network graph """
        for node in self.list_nodes:
            self.attributes[node] = self.network.graph.nodes[node]['attributes']

    def set_attributes_from_dict(self, attr_dict):
        """ Loads attributes from the passed dictionary """
        for node in attr_dict:
            self.network.graph.nodes[node]['attributes'] = attr_dict[node]

    def set_reaction_enzyme_map_from_network(self):
        """ Set dict of which enzyme for each step from default options in network """
        self.reaction_enzyme_map = {}
        for node in self.reactions:
            self.reaction_enzyme_map[node] = self.network.graph.nodes[node]['attributes']['selected_enzyme']

    def visualise(self, height='750px', width='1000px', options=None):
        return visualise.Visualiser(self.network, graph=self.sub_graph).html(height=height, width=width, options=options)

    def get_visjs_nodes_and_edges(self):
        nodes, edges = visualise.Visualiser(self.network, graph=self.sub_graph).nodes_edges()
        return nodes, edges

    def _get_reaction_names(self):
        names = []
        for reaction in self.reactions:
            names.append(self.network.graph.nodes[reaction]['attributes']['name'])
        return names


    """
    This code isnt working
    def calculate_best_cofactor(self, log=False):
        ''' Determines the best combination of enzymes for minimum cofactor stoichimetry '''
        cofactor_count.add_reaction_cofactor_count(self.sub_graph)

        def get_list_enzymes():
            all_enzymes = []
            all_indexes = []

            for node in self.reactions:
                node_enzymes = []
                node_enzyme_indexes = []
                cofactors = self.sub_graph.nodes[node]['attributes']['cofactors']
                for i, enz_dict in enumerate(cofactors):
                    enz = list(enz_dict.keys())[0]
                    node_enzymes.append(enz)
                    node_enzyme_indexes.append(i)

                all_enzymes.append(node_enzymes)
                all_indexes.append(node_enzyme_indexes)

            return all_enzymes, all_indexes

        def get_combinations(all_indexes):
            combination_indexes = list(itertools.product(*all_indexes))
            combination_names = []
            for comb in combination_indexes:
                names = []
                for i, index in enumerate(comb):
                    names.append(all_enzymes[i][index])
                combination_names.append(names)

            return combination_indexes, combination_names

        all_enzymes, all_indexes = get_list_enzymes()
        combination_indexes, combination_names = get_combinations(all_indexes)

        score = {}
        for comb_indexes, comb_names in zip(combination_indexes, combination_names):
            score[comb_indexes] = {}
            for i, node in enumerate(self.reactions):
                cofactors = self.sub_graph.nodes[node]['attributes']['cofactors']
                cofactor_dict = cofactors[comb_indexes[i]][comb_names[i]]

                for cf in cofactor_dict:
                    if cf not in score[comb_indexes]:
                        score[comb_indexes][cf] = 0
                    score[comb_indexes][cf] += cofactor_dict[cf]

            for cf in score[comb_indexes]:
                if score[comb_indexes][cf] < 0:
                    score[comb_indexes][cf] = -score[comb_indexes][cf]

        if log == True:
            print(score)

        score_totals = {}
        min_score = 9999
        comb_names = combination_names[0]
        for i, comb in enumerate(combination_indexes):
            score_totals[comb] = 0
            for cf in score[comb]:

                score_totals[comb] += score[comb][cf]
            if score_totals[comb] < min_score:
                min_score = score_totals[comb]
                comb_names = combination_names[i]

        if log==True:
            print('Best enzyme combination = ' + str(comb_names))
            print('with score =  ' + str(min_score))

        if min_score == 9999:
            min_score = 0

        self.cofactor_score = min_score
    """