import retrobiocat_web.retro.generation.node_analysis
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.visualisation import visualise
from retrobiocat_web.retro.generation.network_generation.network_scoring import NetworkEvaluator
from retrobiocat_web.retro.generation.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
from retrobiocat_web.retro.generation.retrosynthesis_engine.retrorules import RetroRules
import networkx as nx
from retrobiocat_web.retro.generation.load.rxn_class import RetroBioCat_Reactions
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile

class Network(object):

    def __init__(self, graph=nx.DiGraph(), target_smiles='', number_steps=0, max_nodes=False, include_two_step=False, include_experimental=False, print_log=False):
        self.graph = graph
        self.number_steps = number_steps
        self.target_smiles = rdkit_smile(target_smiles, warning=True)
        self.rxn_obj = RetroBioCat_Reactions(include_experimental=include_experimental,
                                             include_two_step=include_two_step)
        self.rxns = self.rxn_obj.rxns

        self.substrate_nodes = []
        self.reaction_nodes = []
        self.end_nodes = []

        self.evaluator = NetworkEvaluator(self, print_log=print_log)
        self.retrosynthesisEngine = RetrosynthesisEngine(self)
        #self.retrorules = RetroRules(self)

        self.settings = {"allow_backwards_steps": False,
                         "add_if_precursor" : True,
                         "print_precursors" : False,
                         "combine_enantiomers" : True,
                         "clean_brackets" : True,
                         "print_log" : print_log,
                         "remove_simple": False,
                         "similarity_score_threshold": 0.6,
                         "num_enzymes": 1,
                         "complexity_score": "SCS",
                         "calculate_complexities" : True,
                         "calculate_substrate_specificity" : True,
                         "get_building_blocks" : True,
                         'building_blocks_db_mode' : 'in_stock',
                         "max_nodes": max_nodes,
                         "prune_steps" : 1,
                         "rank_pathways_by_enzyme" : True,
                         "only_postitive_enzyme_data": False,
                         "colour_substrates" : 'Starting material', #Starting material, Relative complexity, or Off
                         "colour_reactions" : 'Substrate specificity', #Substrate specificity, Complexity change, or Off
                         "colour_arrows" : "None", #None or Complexity change
                         "show_negative_enzymes": True,
                         "display_cofactors" : True,
                         "molSize": (300,300),
                         'prune_on_substrates' : False,
                         'max_reactions' : False,
                         'specificity_score_substrates' : False,
                         'include_experimental' : include_experimental,
                         'include_two_step': include_two_step,
                         'rr_min_diameter': 2,
                         'rr_min_products': 10,
                         'rr_max_reactions': 1,
                         'aizynth_reaction_mode': 'policy',
                         'retrobiocat_reaction_mode': 'complexity',
                         'only_reviewed_activity_data': False}

    def update_settings(self, settings):
        self.settings.update(settings)
        self.evaluator.specficity_scorer.score_substrates = self.settings['specificity_score_substrates']
        self.evaluator.buyable_scorer.mode = self.settings['building_blocks_db_mode']

        if ('include_experimental' in settings) or ('two_step' in settings):
            self.rxn_obj = RetroBioCat_Reactions(include_experimental=self.settings['include_experimental'],
                                                 include_two_step=self.settings['include_two_step'])
            self.rxns = self.rxn_obj.rxns


    def initialise_graph(self, target_smiles=''):

        if target_smiles != '':
            self.target_smiles = rdkit_smile(target_smiles, warning=True)

        self._log(" -initialise graph, target smiles:  " + str(self.target_smiles))

        self.graph = nx.DiGraph()
        self.graph.add_node(self.target_smiles, attributes={'name': self.target_smiles,
                                                            'node_type': 'substrate',
                                                            'node_num': 0,
                                                            'substrate_num' : 1})

        return [self.target_smiles]

    def generate(self, target_smiles, number_of_steps, calculate_scores=True):
        """
        Iteratively applies reaction rules to end nodes, adding reactions and products to the graph.

        If the graph gets bigger than max nodes, will delete reactions accoring to molecular complexity
        """
        self._log("Generating Network")

        self.initialise_graph(target_smiles=target_smiles)
        self.number_steps = number_of_steps
        self.retrosynthesisEngine.generate_network(self.target_smiles, self.number_steps, self.rxns, self.graph)

        if calculate_scores == True:
            self.evaluator.calculate_scores(self)

        self._log('Network Generated')

        return self.graph

    def add_step(self, smiles, calculate_scores=True):
        """
        Add a single retrosynthetic step to graph from a single smiles node
        """

        listSmiles, listReactions = self.retrosynthesisEngine.single_step(smiles,self.rxns,self.graph)
        self.get_node_types()

        if calculate_scores == True:
            self.evaluator.calculate_scores(self)

        return listSmiles, listReactions
    
    def add_chemical_step(self, smiles, calculate_scores=True):
        listSmiles, listReactions = self.retrosynthesisEngine.single_aizynth_step(smiles, self.graph)
        self.get_node_types()

        if calculate_scores == True:
            self.evaluator.calculate_scores(self)

        return listSmiles, listReactions

    def custom_reaction(self, product_smiles, substrate_smiles, reaction_name):
        """ Add a custom reaction to self.graph"""

        product_smiles = rdkit_smile(product_smiles, warning=True)
        new_substrates, new_reactions = self.retrosynthesisEngine.custom_reaction(self.graph, product_smiles, substrate_smiles, reaction_name)

        self.substrate_nodes.extend(new_substrates)
        if product_smiles not in self.substrate_nodes:
            self.substrate_nodes.append(product_smiles)
        self.reaction_nodes.extend(new_reactions)

        self._log('Custom reaction added: ' + str(product_smiles) + '<--' + str(new_reactions) + '<--' + str(new_substrates))

        self.calculate_scores()

        new_substrates.append(product_smiles)

        return new_substrates, new_reactions

    def calculate_scores(self):
        self.evaluator.calculate_scores(self)

    def delete_reaction_node(self, reaction_node_to_remove):
        deleted = self.retrosynthesisEngine.graphPruner.delete_reaction_node(self, reaction_node_to_remove)
        self.get_node_types()
        return deleted

    def visualise(self, height='750px', width='1000px', options=None):
        return visualise.Visualiser(self).html(height=height, width=width, options=options)

    def get_visjs_nodes_and_edges(self, graph=None):
        nodes, edges = visualise.Visualiser(self, graph=graph).nodes_edges()
        return nodes, edges

    def get_node_types(self):
        self.substrate_nodes = node_analysis.get_substrate_nodes(self.graph)
        self.reaction_nodes = node_analysis.get_reaction_nodes(self.graph)
        self.end_nodes = node_analysis.get_nodes_with_no_successors(self.graph)

    def _log(self, to_log):
        if self.settings['print_log'] == True:
            print(to_log)

    def attributes_dict(self):
        att_dict = {}
        for node in list(self.graph):
            att_dict[node] = self.graph.nodes[node]['attributes']

        self.substrate_nodes = node_analysis.get_substrate_nodes(self.graph)
        self.reaction_nodes = node_analysis.get_reaction_nodes(self.graph)
        return att_dict

    def add_attributes(self, attributes_dict):
        for node in attributes_dict:
            self.graph.nodes[node]['attributes'] = attributes_dict[node]

        self.get_node_types()

    def clear_cofactor_data(self):
        for node in list(self.graph):
            if self.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                if 'cofactors' in self.graph.nodes[node]['attributes']:
                    self.graph.nodes[node]['attributes'].pop('cofactors')

if __name__ == '__main__':
    network = Network(max_nodes=100)
    network.generate('O=C(O)C(=O)Cc1ccccc1', 3)

