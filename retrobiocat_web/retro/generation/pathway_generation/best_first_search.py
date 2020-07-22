"""
Module for generating pathways using best first search
"""
import retrobiocat_web.retro.generation.node_analysis
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.generation.pathway_generation.pathway import Pathway
from retrobiocat_web.retro.generation.network_generation.network import Network
import copy
import random

class BFS():

    def __init__(self, network=None, target=None,  max_pathways=50000, max_pathway_length=5, min_weight=1, use_random=False,
                 print_log=False, score_pathways=True):
        """
        Best First Search object, for generating pathways from a network

        After initialising, run search using the .run() method

        Args:
            network: a network object which has been generated
            min_weight: the minimum weight to assign to zero complexity change (and Stop)
            max_pathways: the maximum number of pathways to generate before stopping
            use_random: set the bfs to use weighted random selection rather than always picking the best
        """
        self.score_pathways = score_pathways
        self.print_log = print_log
        self.min_weight = min_weight
        self.choices = {}
        self.max_pathways = max_pathways
        self.max_pathway_length = max_pathway_length
        self.pathways = []
        self.use_random = use_random
        self.network = network
        self.generate_network = False
        if self.network == None:
            self.target = node_analysis.rdkit_smile(target, warning=True)
            self.generate_network = True
            self.network = Network(target_smiles=self.target, number_steps=self.max_pathway_length, print_log=False)
            self.network.generate(self.target, 0)
            self.log('BFS - will generate network')
        else:
            self.target = self.network.target_smiles

    def log(self, msg):
        if self.print_log == True:
            print(msg)

    def _get_context(self, nodes):
        """ Returns the pathway context, which is a string of node numbers"""
        list_node_numbers = []
        context = ''
        for node in nodes:
            list_node_numbers.append(self.network.graph.nodes[node]['attributes']['node_num'])

        sorted_node_numbers = sorted(list_node_numbers)
        for node_num in sorted_node_numbers:
            context += str(node_num)
            context += '-'
        return context

    def _expand_network(self, smi):
        nodes_added = []
        new_substrates, new_reactions = self.network.add_step(smi)
        nodes_added.extend(new_substrates)
        nodes_added.extend(new_reactions)
        return nodes_added

    def _get_choices(self, end_nodes):
        """ Returns a list of reaction nodes (and Stop) which are choices for the next step"""

        def get_choice_scores(choices):
            scores = [0]
            for node in choices[1:]:
                scores.append(self.network.graph.nodes[node]['attributes']['change_in_complexity'])
            return scores

        def get_weighted_scores(scores):
            # invert changes so decreases in complexity are favoured
            inverted_reaction_complexity_changes = [x * -1 for x in scores]

            min_change = min(inverted_reaction_complexity_changes)
            if min_change < 0:
                min_change = -min_change
            else:
                min_change = 0

            non_neg_changes = [x + self.min_weight + min_change for x in inverted_reaction_complexity_changes]

            return non_neg_changes

        def get_choices(end_nodes, graph):
            successor_reactions = ['Stop']
            for node in end_nodes:
                successor_reactions.extend(list(graph.successors(node)))
            return successor_reactions

        def make_choice_dict(choices, scores):
            choice_dict = {}
            for i, choice in enumerate(choices):
                choice_dict[choice] = scores[i]
            return choice_dict

        choices = get_choices(end_nodes, self.network.graph)
        scores = get_choice_scores(choices)
        weighted_scores = get_weighted_scores(scores)
        choice_dict = make_choice_dict(choices, weighted_scores)

        return choice_dict

    def _pick_choice(self, context):
        """ Given a context, picks an option to extend (or stop) that pathway """
        def pick_best(choices, scores):
            sorted_options = node_analysis.sort_by_score(choices, scores, reverse=False)
            return sorted_options[0]

        def pick_weighted_random(choices, scores):
            return random.choices(choices, scores, k=1)[0]

        def get_lists_choices_scores(choices_dict):
            list_choices = []
            list_scores = []
            for choice in choices_dict:
                list_choices.append(choice)
                list_scores.append(choices_dict[choice])

            return list_choices, list_scores

        choices, scores = get_lists_choices_scores(self.choices[context])

        if self.use_random == False:
            option = pick_best(choices, scores)
        else:
            option = pick_weighted_random(choices, scores)

        return option

    def _add_reaction(self, reaction_choice):
        new_end_nodes = list(self.network.graph.successors(reaction_choice))
        added_nodes = [reaction_choice] + new_end_nodes

        return added_nodes, new_end_nodes

    def _check_pathway_has_end(self, nodes):
        pathway_subgraph = self.network.graph.subgraph(nodes)
        end_nodes = node_analysis.get_nodes_with_no_successors(pathway_subgraph)
        if len(end_nodes) == 0:
            return False
        return True

    def _make_pathway(self, nodes):
        """ Create pathway object from list of nodes"""
        return Pathway(nodes, self.network, calc_scores=self.score_pathways)

    def _check_if_should_expand_network(self, end_nodes, pathway_nodes):
        if self.generate_network == True:
            if self._num_reactions(pathway_nodes) < self.max_pathway_length:
                for node in end_nodes:
                    if len(list(self.network.graph.successors(node))) == 0:
                        self._expand_network(node)

    def _is_node_already_in_pathway(self, current_nodes, new_nodes):
        for node in new_nodes:
            if node in current_nodes:
                return True
        return False

    def _num_reactions(self, nodes):
        count = 0
        for node in nodes:
            if self.network.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                count += 1
        return count

    def run(self):
        """
        Generate pathways using best first search

        Returns: list of pathways
        """
        self.log('Run BFS')
        self.pathways = []
        self.choices = {}

        nodes = [self.target]
        context = self._get_context(nodes)

        self._check_if_should_expand_network(nodes, nodes)
        self.choices[context] = self._get_choices(nodes)
        start_context = copy.deepcopy(context)

        while (len(self.pathways) < self.max_pathways) and (len(self.choices[start_context]) > 0):
            nodes = [self.network.target_smiles]
            context = self._get_context(nodes)

            while len(self.choices[context]) > 0:
                best_choice = self._pick_choice(context)

                if best_choice == 'Stop':
                    if self._check_pathway_has_end(nodes) == True:
                        self.pathways.append(nodes)
                    self.choices[context].pop('Stop')
                    break

                else:
                    added_nodes, new_end_nodes = self._add_reaction(best_choice)

                    if self._is_node_already_in_pathway(nodes, added_nodes) == True:
                        self.choices[context].pop(best_choice)
                        break
                    else:
                        new_context = self._get_context(nodes+added_nodes)
                        if new_context not in self.choices:
                            self._check_if_should_expand_network(new_end_nodes, nodes+added_nodes)
                            self.choices[new_context] = self._get_choices(new_end_nodes)

                        if len(self.choices[new_context]) == 0:
                            self.choices[context].pop(best_choice)
                        else:
                            nodes = nodes+added_nodes
                            context = new_context
        self.log('BFS complete')
        return self.pathways

    def get_pathways(self):
        pathway_objects = []
        for list_nodes in self.pathways:
            pathway_objects.append(self._make_pathway(list_nodes))
        return pathway_objects


if __name__ == '__main__':
    def bfs_with_network_generation():
        network = Network()
        network.generate('CCCCC=O', 5)

        bfs = BFS(network=network, max_pathways=10000)
        bfs.run()
        pathways = bfs.get_pathways()
        return pathways

    def bfs_without_network_generation():
        # currently produces less pathways than starting with the complete network, not sure why.
        target = 'CCCCC=O'
        bfs = BFS(target=target, max_pathway_length=5, max_pathways=10000)
        pathways = bfs.run()

        return pathways


    p1 = bfs_with_network_generation()
    p2 = bfs_without_network_generation()

