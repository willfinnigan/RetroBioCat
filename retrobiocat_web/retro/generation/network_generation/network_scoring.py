from retrobiocat_web.retro.generation.load import load_rule_yamls
from retrobiocat_web.retro.evaluation import complexity
from retrobiocat_web.retro.evaluation.scscore.standalone_model_numpy import SCScorer
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.enzyme_identification.molecular_similarity import SubstrateSpecificityScorer
from retrobiocat_web.retro.evaluation.starting_material import StartingMaterialEvaluator

class NetworkEvaluator():

    def __init__(self, network, print_log=False):

        self.print_log = print_log
        self.sc_score_model = SCScorer()
        self.sc_score_model.restore()
        self.buyable_scorer = StartingMaterialEvaluator()
        self.specficity_scorer = SubstrateSpecificityScorer()

        self.enzyme_reaction_map = network.rxn_obj.reaction_enzyme_map
        self.reactionEnzymeCofactorDict = network.rxn_obj.reactionEnzymeCofactorDict

    def calculate_scores(self, network):

        self._log("Calculating scores")
        self.add_scores_substrate_availability(network)
        self.add_enzymes(network)
        self.add_cofactors(network)
        self.add_scores_is_enzyme(network)
        self.add_scores_complexity(network)
        self.add_scores_specificity(network)

        self._log("Scores calculated")

    def _log(self, msg):
        if self.print_log == True:
            print(msg)

    def add_scores_complexity(self, network):
        if network.settings['calculate_complexities'] == True:
            self._log('- molecular complexity')
            complexity.add_scscore(network.graph, network.substrate_nodes, self.sc_score_model)
            complexity.add_change_in_complexity(network.graph)
            complexity.add_relative_complexity(network.graph, network.target_smiles)
            complexity.add_reaction_relative_complexity(network.graph, network.target_smiles)

    def add_enzymes(self, network):
        enzyme_map = self.enzyme_reaction_map

        for node in network.reaction_nodes:
            if 'possible_enzymes' not in network.graph.nodes[node]['attributes']:
                reaction = network.graph.nodes[node]['attributes']['name']
                if reaction in enzyme_map:
                    possible_enzymes = enzyme_map[reaction]
                else:
                    possible_enzymes = ['']
                network.graph.nodes[node]['attributes']['possible_enzymes'] = possible_enzymes
                network.graph.nodes[node]['attributes']['selected_enzyme'] = possible_enzymes[0]

        return network

    def add_cofactors(self, network):
        """
        Add a dict of enzymes and their cofactors to graph under as {Enz : ['Cofactor +', 'Cofactor -']}
        """
        self._log('- adding cofactors')

        for node in list(network.graph):
            if network.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                if 'enzyme_cofactors' not in network.graph.nodes[node]['attributes']:
                    reaction_name = network.graph.nodes[node]['attributes']['name']
                    if reaction_name in self.reactionEnzymeCofactorDict:
                        enzyme_cofactors = self.reactionEnzymeCofactorDict[reaction_name]
                    else:
                        enzyme_cofactors = {}

                    if len(list(enzyme_cofactors.keys())) > 0:
                        network.graph.nodes[node]['attributes']['selected_enzyme'] = list(enzyme_cofactors.keys())[0]

                    network.graph.nodes[node]['attributes']['enzyme_cofactors'] = enzyme_cofactors

    def add_scores_is_enzyme(self, network):
        self._log('- is enzyme')
        for node in network.reaction_nodes:
            if 'is_enzyme' not in network.graph.nodes[node]['attributes']:
                reaction_name = network.graph.nodes[node]['attributes']['name']
                if reaction_name in self.enzyme_reaction_map:
                    enzymes = self.enzyme_reaction_map[reaction_name]
                    if enzymes[0] == 'Chemical':
                        network.graph.nodes[node]['attributes']['is_enzyme'] = 0
                    else:
                        network.graph.nodes[node]['attributes']['is_enzyme'] = 1
                else:
                    network.graph.nodes[node]['attributes']['is_enzyme'] = 0

    def add_scores_specificity(self, network):
        def determine_subOne_subTwo(listSubstrateSmiles):
            subOne = None
            subTwo = None
            for smi in listSubstrateSmiles:
                if network.graph.nodes[smi]['attributes']['substrate_num'] == 1:
                    subOne = smi
                elif network.graph.nodes[smi]['attributes']['substrate_num'] == 2:
                    subTwo = smi
            return subOne, subTwo

        if network.settings['calculate_substrate_specificity'] == True:
            self._log('- specificity scores')

            only_active = network.settings["only_postitive_enzyme_data"]
            sim_cuttoff = network.settings['similarity_score_threshold']


            for node in network.reaction_nodes:
                if 'specificity_scores' not in network.graph.nodes[node]['attributes']:
                    network.graph.nodes[node]['attributes']['specificity_scores'] = {}
                    network.graph.nodes[node]['attributes']['enzyme_info'] = {}

                    possible_enzymes = network.graph.nodes[node]['attributes']['possible_enzymes']
                    reaction_name = network.graph.nodes[node]['attributes']['name']

                    product = list(network.graph.predecessors(node))[0]
                    subOne, subTwo = determine_subOne_subTwo(list(network.graph.successors(node)))

                    for enz in possible_enzymes:
                        score, info, = self.specficity_scorer.scoreReaction(reaction_name, enz, product, subOne, subTwo,
                                                                           sim_cutoff=sim_cuttoff,
                                                                           onlyActive=only_active)

                        network.graph.nodes[node]['attributes']['specificity_scores'][enz] = score
                        network.graph.nodes[node]['attributes']['enzyme_info'][enz] = info

        self.select_best_enzyme(network)

    def add_scores_substrate_availability(self, network):
        if network.settings['get_building_blocks'] == True:

            self._log('- substrate availability')

            substrate_nodes = node_analysis.get_substrate_nodes(network.graph)

            for node in substrate_nodes:
                if 'is_starting_material' not in network.graph.nodes[node]['attributes']:
                    network.graph.nodes[node]['attributes']['is_starting_material'] = self.buyable_scorer.eval(node)

    def select_best_enzyme(self, network):
        if network.settings['calculate_substrate_specificity'] == True:
            for node in network.reaction_nodes:
                current_enz = network.graph.nodes[node]['attributes']['selected_enzyme']
                current_score = network.graph.nodes[node]['attributes']['specificity_scores'][current_enz]
                current_score_neg = True
                possible_enzymes = network.graph.nodes[node]['attributes']['possible_enzymes']

                for enz in possible_enzymes:
                    score = network.graph.nodes[node]['attributes']['specificity_scores'][enz]
                    if (score > current_score and score != 0) or (abs(score) > current_score and current_score_neg==True):
                        network.graph.nodes[node]['attributes']['selected_enzyme'] = enz
                        current_score = abs(score)
                        if score < 0:
                            current_score_neg = True


if __name__ == '__main__':
    from retrobiocat_web.retro.generation.network_generation.network import Network
    network = Network()
    network.generate('CCCCC=O', 4)  # evaluator calculate scores called during generate
