import numpy as np
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.generation.pathway_generation.pathway import Pathway
from retrobiocat_web.retro.generation.network_generation.network import Network



class ReactionExpander():

    def __init__(self, network, maxLength):
        self.network = network
        self.maxLength = maxLength
        self.buyable_stop_score = 1

        self.min_weight = 1

    def expand_network(self, endNodes):
            for node in endNodes:
                if node_analysis.get_num_reactions_from_target(self.network.graph, node, self.network.target_smiles) < self.maxLength:
                    if len(list(self.network.graph.successors(node))) == 0:
                        self.network.add_step(node)

    def get_reaction_choices(self, endNodes):
        """ Returns a list of reaction nodes (and Stop) which are choices for the next step"""

        def get_choice_scores(choices):

            if self._are_end_nodes_buyable == True:
                scores = [self.buyable_stop_score]
            else:
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

        choices = get_choices(endNodes, self.network.graph)
        scores = get_choice_scores(choices)
        weighted_scores = get_weighted_scores(scores)
        choice_dict = make_choice_dict(choices, weighted_scores)

        return choice_dict

    def _are_end_nodes_buyable(self, endNodes):
        for node in endNodes:
            if self.network.graph.nodes[node]['attributes']['is_starting_material'] != 1:
                return False
        return True


class Rollout():

    def __init__(self,  currentTreeNode, network, mcts):
        self.network = network
        self.currentTreeNode = currentTreeNode
        self.currentEndNodes = currentTreeNode.endNodes
        self.currentReactionChoices = {}
        self.pathwayNodes = currentTreeNode.pathwayNodes
        self.maxLength = mcts.max_length
        self.minLength = mcts.min_length

        self.reactionExpander = ReactionExpander(network, self.maxLength)

    def run(self):

        while True:
            # check if should stop
            if node_analysis.get_current_pathway_length(self.network.graph,
                                                        self.currentEndNodes,
                                                        self.network.target_smiles) >= self.maxLength:
                return self.score_pathway()

            # expand network if necessary, then get choices
            self.reactionExpander.expand_network(self.currentEndNodes)
            self.currentReactionChoices = self.reactionExpander.get_reaction_choices(self.currentEndNodes)

            # pick a node, add to pathway
            chosenReaction = self._pick_best_choice()
            if chosenReaction == 'Stop':
                return self.score_pathway()

            added_nodes, new_end_nodes = self._get_new_pathway_nodes(chosenReaction)
            if self._is_node_already_in_pathway(added_nodes):
                return self.score_pathway()
            else:
                self.pathwayNodes.extend(added_nodes)
                self.currentEndNodes = new_end_nodes
                self.currentReactionChoices = {}

    def score_pathway(self):
        pathway = Pathway(self.pathwayNodes, self.network)

        complexity_score = pathway.scores.change_in_complexity
        buyable_score = pathway.scores.starting_material*2

        score = complexity_score + buyable_score

        if len(pathway.reactions) < self.minLength:
            score -= 2

        return score

    def _pick_best_choice(self):
        """ Given a context, picks an option to extend (or stop) that pathway """
        def pick_best(choices, scores):
            sorted_options = node_analysis.sort_by_score(choices, scores, reverse=False)
            return sorted_options[0]

        def get_lists_choices_scores(choices_dict):
            list_choices = []
            list_scores = []
            for choice in choices_dict:
                list_choices.append(choice)
                list_scores.append(choices_dict[choice])

            return list_choices, list_scores

        choices, scores = get_lists_choices_scores(self.currentReactionChoices)
        option = pick_best(choices, scores)

        return option

    def _get_new_pathway_nodes(self, chosenReaction):
        new_end_nodes = list(self.network.graph.successors(chosenReaction))
        added_nodes = [chosenReaction] + new_end_nodes

        return added_nodes, new_end_nodes

    def _is_node_already_in_pathway(self, new_nodes):
        for node in new_nodes:
            if node in self.pathwayNodes:
                return True
        return False


class MCTS():

    def __init__(self, target_smiles, exploration=1):

        self.network = Network(target_smiles=target_smiles)
        self.network.initialise_graph()

        self.smilesWhichHaveBeenExpanded = []
        self.min_weight = 0.5
        self.root = MCTS_Node([self.network.target_smiles], [self.network.target_smiles], self.network, self)
        self.treeNodes = [self.root]
        self.exploration=exploration
        self.max_length = 10
        self.min_length = 2

        self.reactionExpander = ReactionExpander(self.network, self.max_length)

    def calculate_UCB(self, treeNode, exploration):
        N = treeNode.parent.visits
        ni = treeNode.visits
        C = exploration

        wi = treeNode.value

        if (ni == 0) or (N == 0):
            return np.inf

        ucb = (wi/ni) + (C*(np.sqrt(np.log(N)/ni)))

        return ucb

    def selection(self, currentTreeNode):
        print('--Selection--')

        max_score = 0
        while len(currentTreeNode.children) != 0:
            max_score = 0
            for child in currentTreeNode.children:
                ucb = self.calculate_UCB(child, self.exploration)
                print(str(child.pathwayNodes) + ' - ' + str(ucb))
                if ucb > max_score:
                    currentTreeNode = child
                    max_score = ucb
            if max_score == 0:
                currentTreeNode = currentTreeNode.children[0]

        print('Selected:')
        print(str(currentTreeNode.pathwayNodes))
        print(str(max_score))

        return currentTreeNode

    def rollout(self, currentTreeNode):
        print('--Rollout--')
        score = Rollout(currentTreeNode, self.network, self).run()
        print('score = ' + str(score))
        return score

    def backpropagate(self, currentTreeNode, score):
        print('--Backpropagage--')
        currentTreeNode.value += score
        currentTreeNode.visits += 1
        while currentTreeNode.parent is not None:
            currentTreeNode = currentTreeNode.parent
            currentTreeNode.value += score
            currentTreeNode.visits += 1

    def run(self):
        for i in range(10000):
            currentTreeNode = self.selection(self.root)
            currentTreeNode.expand()
            currentTreeNode = self.selection(currentTreeNode)
            score = self.rollout(currentTreeNode)
            self.backpropagate(currentTreeNode, score)

    def get_best_pathways(self, n):
        def non_picked_children(currentTreeNode):
            children = []
            for child in currentTreeNode.children:
                if child not in currentTreeNode.picked_children:
                    children.append(child)
            return children

        pathways = []
        for i in range(n):
            currentTreeNode = self.root
            children = non_picked_children(currentTreeNode)

            while len(children) != 0:
                max_score = 0
                for child in children:
                    if child.visits == 0:
                        score = 0
                    else:
                        score = child.value / child.visits

                    if score >= max_score:
                        currentTreeNode = child
                        max_score = score
                    else:
                        currentTreeNode = children[0]

                children = non_picked_children(currentTreeNode)

            currentTreeNode.parent.picked_children.append(currentTreeNode)
            pathways.append(currentTreeNode.pathwayNodes)

        return pathways







class MCTS_Node():

    def __init__(self, pathwayNodes, endNodes, network, mcts, parent=None):
        self.network = network
        self.pathwayNodes = pathwayNodes
        self.endNodes = endNodes
        self.mcts = mcts
        self.visits = 0
        self.value = 0
        self.terminal_visits = 0
        self.parent = parent

        self.children = []
        self.reactions = {} # reactions with complexity scores
        self.children_reactions = {} # key is child node, value is reaction to get there
        self.reaction_children = {} # key is reaction, value is child node
        self.parent=parent

        self.terminal = False
        self.picked_children = []

    def expand(self):
        print('--Expand Node--')
        if self.terminal != True:

            self.mcts.reactionExpander.expand_network(self.endNodes)
            self.reactions = self.mcts.reactionExpander.get_reaction_choices(self.endNodes)
            print(self.reactions)

            for reaction in self.reactions:
                new_mcts_node = self._create_new_mcst_node(reaction)
                self.children_reactions[new_mcts_node] = reaction
                self.reaction_children[reaction] = new_mcts_node
                self.children.append(new_mcts_node)

    def _create_new_mcst_node(self, reaction):
        additionalPathwayNodes, newEndNodes = self._select_reaction(reaction)
        new_pathway_nodes = self.pathwayNodes+additionalPathwayNodes
        new_MCTS_node = MCTS_Node(new_pathway_nodes, newEndNodes, self.network, self.mcts, parent=self)
        if newEndNodes == []:
            new_MCTS_node.terminal = True
        return new_MCTS_node

    def _select_reaction(self, reaction_choice):
        if reaction_choice == 'Stop':
            new_end_nodes = []
            added_nodes = []
        else:
            new_end_nodes = list(self.network.graph.successors(reaction_choice))
            added_nodes = [reaction_choice] + new_end_nodes

        return added_nodes, new_end_nodes




if __name__ == "__main__":
    mcts = MCTS('[C@H]1(C2=CC=CC=C2)NCCCC1')
    mcts.run()
    pathways = mcts.get_best_pathways(200)

    for p in pathways:
        print(p)
        print()

