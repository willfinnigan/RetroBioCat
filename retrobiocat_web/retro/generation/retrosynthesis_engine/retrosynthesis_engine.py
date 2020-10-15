from rdkit.Chem import AllChem, rdmolops
from rdkit import Chem
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.retro.generation import node_analysis
import uuid
import networkx as nx
from retrobiocat_web.retro.rdchiral.main import rdchiralReactants, rdchiralRun
from retrobiocat_web.retro.rdchiral.clean import combine_enantiomers_into_racemic
from retrobiocat_web.retro.generation.retrosynthesis_engine import aizynthfinder_actions



class BracketCleaner():

    def __init__(self):
        # smi = clean_brackets("CCC[C]CCC")
        # print(smi)  #--> CCCCCCC

        self.bracket_required_chars = ['@', '+', '-', 'H', 'h']
        self.bracket_remove = ['[', ']']

    def clean_brackets(self, smi):
        """
        If smi contains '[', check if brackets are needed, if not then remove.
        """

        if '[' in smi:
            bracket_locations = self._get_bracket_locations(smi)

            for loc in bracket_locations:
                if self._check_bracket_required(smi, loc) == False:
                    smi = self._remove_bracket(smi, loc)

        return smi

    def _get_bracket_locations(self, smi):
        brackets = []
        for i, l in enumerate(smi):
            if l == '[':
                for j, l2 in enumerate(smi[i:]):
                    if l2 == ']':
                        brackets.append([i, j + i + 1])
                        # print(smi[i:j+i+1])
        return brackets

    def _check_bracket_required(self, smi, loc):
        bracket = smi[loc[0]:loc[1]]

        for l in bracket:
            if l in self.bracket_required_chars:
                return True

        return False

    def _remove_bracket(self, smi, loc):
        new_brac = ''
        for l in smi[loc[0]:loc[1]]:
            if l not in self.bracket_remove:
                new_brac += l

        new_smi = smi[0:loc[0]] + new_brac + smi[loc[1]:]
        return new_smi

class ReactionSelector():

    def __init__(self, network, graphPruner):
        self.network = network
        self.graphPruner = graphPruner

    def remove_disallowed_products(self, disallowedProducts, listProducts, listReactions):
        if disallowedProducts == None:
            return listProducts, listReactions

        listProductsToKeep = listProducts
        self.network.get_node_types()
        for product in listProducts:
            if product in disallowedProducts:
                for reaction in list(self.network.graph.predecessors(product)):
                    if reaction in listReactions:
                        self.graphPruner.delete_terminal_reaction_node(self.network, reaction)
                        listReactions.remove(reaction)
                if product in listProductsToKeep:
                    listProductsToKeep.remove(product)
        return listProductsToKeep, listReactions

    def select_best_by_complexity(self, listProducts, listReactions, max_reactions):
        if (len(listReactions) < max_reactions) or (max_reactions == False):
            return listProducts, listReactions

        self.network.get_node_types()
        self.network.evaluator.add_scores_complexity(self.network)
        changes_in_complexity = []
        for reaction in listReactions:
            changes_in_complexity.append(self.network.graph.nodes[reaction]['attributes']['change_in_complexity'])

        sorted_reactions = node_analysis.sort_by_score(listReactions, changes_in_complexity, reverse=True)
        for reaction in sorted_reactions[max_reactions:]:
            self.graphPruner.delete_terminal_reaction_node(self.network, reaction)

        reactions_to_keep = sorted_reactions[0:max_reactions]
        smiles_to_keep = []
        for reaction in reactions_to_keep:
            for smi in list(self.network.graph.successors(reaction)):
                smiles_to_keep.append(smi)


        return smiles_to_keep, reactions_to_keep

    def select_best_aizynth_by_policy(self, listProducts, listReactions, max_reactions):
        if (len(listReactions) < max_reactions) or (max_reactions == False):
            return listProducts, listReactions

        policy_scores = []
        for reaction in listReactions:
            policy_scores.append(self.network.graph.nodes[reaction]['attributes']['metadata']['policy_probability'])

        sorted_reactions = node_analysis.sort_by_score(listReactions, policy_scores, reverse=True)
        for reaction in sorted_reactions[max_reactions:]:
            self.graphPruner.delete_terminal_reaction_node(self.network, reaction)

        reactions_to_keep = sorted_reactions[0:max_reactions]
        smiles_to_keep = []
        for reaction in reactions_to_keep:
            for smi in list(self.network.graph.successors(reaction)):
                smiles_to_keep.append(smi)

        return smiles_to_keep, reactions_to_keep

class RuleApplicator():

    def __init__(self, network):
        self.network = network
        self.small_precursors = ['N', 'O', 'O=O', 'H+', '[H+]']
        self.bracket_cleaner = BracketCleaner()

    def run(self, smile, rxns, graph, explicit_hydrogens=False):
        precursor_dict = self.apply_rules(smile, rxns)
        precursor_dict = self._remove_precursors_already_in_graph(graph, smile, precursor_dict)

        if self.network.settings['remove_simple'] == True:
            for name in precursor_dict:
                precursor_dict[name] = self._remove_simple_precursors(precursor_dict[name])

        return precursor_dict

    def _split_products(self, smi):
        mol = AllChem.MolFromSmiles(smi)
        splitMols = rdmolops.GetMolFrags(mol, asMols=True)
        split_list = []
        for mol in splitMols:
            p_smile = AllChem.MolToSmiles(mol)
            split_list.append(p_smile)
        return split_list

    def _rdkit_products(self, listSmi):
        new_list = []
        for smi in listSmi:
            new_list.append(rdkit_smile(smi))
        return new_list

    def apply_rules(self, smile, rxns):

        reactants = rdchiralReactants(smile)

        rxn_product_dict = {}
        for rxn_name, rxn_list in rxns.items():
            reaction_products = []
            for rxn in rxn_list:
                try:
                    reaction_products.extend(rdchiralRun(rxn, reactants, combine_enantiomers=False))
                except:
                    reaction_products = []
                    print('Error running reactants for: ' + str(smile) + ' ' + str(rxn_name))

            if self.network.settings["combine_enantiomers"] == True:
                reaction_products = combine_enantiomers_into_racemic(set(reaction_products))

            parsed_reaction_products = []
            for smi in reaction_products:
                if self.network.settings["clean_brackets"] == True:
                    smi = self.bracket_cleaner.clean_brackets(smi)

                if self._check_valid_smile(smi, rxn_name=rxn_name) == True:
                    products = self._split_products(smi)
                    products = self._rdkit_products(products)
                    parsed_reaction_products.append(products)

            if len(parsed_reaction_products) != 0:
                rxn_product_dict[rxn_name] = parsed_reaction_products

        return rxn_product_dict

    def _check_valid_smile(self, smile, rxn_name='', log=False):
        if smile == None:
            if log == True:
                print('WARNING: smile = None - ' + str(rxn_name))
            return False

        m = AllChem.MolFromSmiles(smile)
        if m == None:
            if log == True:
                print('WARNING: invalid smile = ' + str(smile) + ' - ' + str(rxn_name))
            return False
        else:
            return True

    def _remove_precursors_already_in_graph(self, graph, origin_smile, precursor_dict):
        # {'name' : list_precursors}
        # for substrate in list,
        new_precursor_dict = {}

        current_origin_reactions = []
        for reaction in list(graph.successors(origin_smile)):
            current_origin_reactions.append(graph.nodes[reaction]['attributes']['name'])

        for name in precursor_dict:
            if name not in current_origin_reactions:
                new_precursor_dict[name] = precursor_dict[name]

        return new_precursor_dict

    def _remove_simple_precursors(self, list_smiles, list_to_remove='default'):
        if list_to_remove == 'default':
            list_to_remove = self.small_precursors

        list_to_return = []
        for smiles in list_smiles:
            new_list = []
            for smile in smiles:
                if smile not in list_to_remove:
                    new_list.append(smile)
            if len(new_list) > 0:
                list_to_return.append(new_list)
        return list_to_return

class GraphManipulator():

    def __init__(self, network):
        self.network = network

    def add_custom_reaction(self, graph, product_smiles, substrate_smiles, reaction_name):
        p_mol = Chem.MolFromSmiles(product_smiles)
        p_smile = Chem.MolToSmiles(p_mol)

        if 'node_type' not in graph.nodes[p_smile]['attributes']:
            graph.nodes[p_smile]['attributes']['node_type'] = 'substrate'

        reaction_node = self._add_reaction_node_to_graph(graph, reaction_name, p_smile, 'custom', {})

        s_mol = Chem.MolFromSmiles(substrate_smiles)
        s_frags = Chem.GetMolFrags(s_mol, asMols=True)
        list_s_smiles = []
        for i, m in enumerate(s_frags):
            s_smiles = Chem.MolToSmiles(m)
            list_s_smiles.append(s_smiles)

            graph.add_node(s_smiles, attributes={'name': s_smiles,
                                                 'node_type': 'substrate',
                                                 'node_num': self._get_node_number(s_smiles, graph),
                                                 'substrate_num': i + 1})
            graph.add_edge(reaction_node, s_smiles)

        return list_s_smiles, [reaction_node]

    def add_nodes_to_graph(self, rxnsSubstratesDict, target_smi, graph, reaction_type, metadata={}):
        list_products, list_reactions = [], []
        for name in rxnsSubstratesDict:
            for precursor_list in rxnsSubstratesDict[name]:
                if self._check_if_reaction_goes_backwards(graph, precursor_list, target_smi) == False:
                    reaction_node = self._add_reaction_node_to_graph(graph, name, target_smi, reaction_type, metadata)
                    if reaction_node not in list_reactions:
                        list_reactions.append(reaction_node)

                    reaction_products = self._add_substrate_node_to_graph(graph, precursor_list, reaction_node)
                    for product in reaction_products:
                        if product not in list_products:
                            list_products.append(product)
        return list_products, list_reactions

    def _add_reaction_node_to_graph(self, graph, reaction_name, source_smile, reaction_type, metadata):
        unique_reaction_name = reaction_name + '(' + str(uuid.uuid1()) + ')'
        graph.add_node(unique_reaction_name, attributes={'name': reaction_name,
                                                         'node_type': 'reaction',
                                                         'reaction_type': reaction_type,
                                                         'metadata': metadata.get(reaction_name, {}),
                                                         'node_num': self._get_node_number(unique_reaction_name, graph)})
        graph.add_edge(source_smile, unique_reaction_name)
        return unique_reaction_name

    def _add_substrate_node_to_graph(self, graph, substrate_names, source_reaction):
        products = []
        for i, precursor in enumerate(substrate_names):
            graph.add_node(precursor, attributes={'name': precursor,
                                                  'node_type': 'substrate',
                                                  'node_num': self._get_node_number(precursor, graph),
                                                  'substrate_num': i + 1})
            graph.add_edge(source_reaction, precursor)
            products.append(precursor)
        return products

    def _get_node_number(self, node, graph):
        if node in list(graph.nodes()):
            return graph.nodes[node]['attributes']['node_num']
        else:
            return len(list(graph.nodes())) + 1

    def _check_if_reaction_goes_backwards(self, graph, reactionProducts, targetSmi):
        are_predecessor = node_analysis.check_substrates_nx_predecessor(graph, reactionProducts, targetSmi)
        if (are_predecessor == True) and (self.network.settings['allow_backwards_steps'] == False):
            return True
        return False

class GraphPruner():

    def __init__(self, print_log=False):
        self.print_log = print_log

    def prune(self, network, smiles):
        def check_if_nodes_need_to_be_removed(network):
            if (network.settings['max_nodes'] != False) and (len(network.substrate_nodes) > network.settings['max_nodes']):
                return True
            return False

        def get_smiles_not_removed(removed, smiles):
            """Return list of smiles which are not in list removed"""
            not_removed = []
            for smi in smiles:
                if smi not in removed:
                    not_removed.append(smi)
            return not_removed

        if check_if_nodes_need_to_be_removed(network):
            removed_nodes = self._run_prune(network, steps=1)
            not_removed = get_smiles_not_removed(removed_nodes, smiles)
            return not_removed

        else:
            return smiles

    def delete_reaction_node(self, network, node_to_remove):
        if node_to_remove not in network.reaction_nodes:
            print('Can not delete non-reaction node')
            return []

        to_delete = [node_to_remove]
        network.graph.remove_node(node_to_remove)
        network.reaction_nodes.remove(node_to_remove)

        # if any node is now not connected to the target, delete also
        connected_nodes = list(nx.dfs_preorder_nodes(network.graph, source=network.target_smiles))

        while len(connected_nodes) != len(list(network.graph.nodes)):
            for node in list(network.graph.nodes):
                if node not in connected_nodes:
                    to_delete.append(node)
                    network.graph.remove_node(node)
                    if node in network.reaction_nodes:
                        network.reaction_nodes.remove(node)
                    if node in network.substrate_nodes:
                        network.substrate_nodes.remove(node)
            connected_nodes = list(nx.dfs_preorder_nodes(network.graph, source=network.target_smiles))

        return to_delete

    def delete_terminal_reaction_node(self, network, node_to_remove):
        if node_to_remove not in network.reaction_nodes:
            print('Can not delete non-reaction node')
            return []

        substrate_successors = list(network.graph.successors(node_to_remove))

        for substrate in substrate_successors:
            if len(list(network.graph.successors(substrate))) != 0:
                return []

        to_delete = [node_to_remove]

        network.graph.remove_node(node_to_remove)
        network.reaction_nodes.remove(node_to_remove)

        for node in substrate_successors:
            if len(list(network.graph.predecessors(node))) == 0:
                network.graph.remove_node(node)
                network.substrate_nodes.remove(node)
                to_delete.append(node)

        return to_delete

    def _get_nodes_to_prune(self, graph, num_to_remove, on_substrates=False):
        def prune_on_reactions(end_nodes, num_to_remove):
            end_node_reactions = node_analysis.get_reaction_nodes_of_list_substrates(graph, end_nodes)
            end_node_reactions = self._check_other_substrates_of_end_node_reactions(graph, end_node_reactions)
            scores = []
            for reaction in end_node_reactions:
                scores.append(graph.nodes[reaction]['attributes']['change_in_complexity'])

            sorted_reactions = node_analysis.sort_by_score(end_node_reactions, scores)
            return sorted_reactions[0:num_to_remove]

        def prune_on_substrates(end_nodes, num_to_remove):
            scores = []
            for substrate in end_nodes:
                scores.append(graph.nodes[substrate]['attributes']['complexity'])
            sorted_substrates = node_analysis.sort_by_score(end_nodes, scores, reverse=True)
            reactions_to_prune = []
            while len(reactions_to_prune) < num_to_remove:
                reactions = node_analysis.get_reaction_nodes_of_list_substrates(graph, [sorted_substrates.pop(0)])
                reactions_to_prune.extend(reactions)
                reactions_to_prune = self._check_other_substrates_of_end_node_reactions(graph, reactions_to_prune)
            return reactions_to_prune

        end_nodes = node_analysis.get_nodes_with_no_successors(graph)
        if on_substrates == True:
            return prune_on_substrates(end_nodes, num_to_remove)
        else:
            return prune_on_reactions(end_nodes, num_to_remove)

    def _check_other_substrates_of_end_node_reactions(self, graph, end_node_reactions):
        def are_all_substrates_terminal(reaction):
            substrate_successors = list(graph.successors(reaction))
            for substrate in substrate_successors:
                if len(list(graph.successors(substrate))) != 0:
                    return False
            return True

        reactions_to_keep = []
        for reaction in end_node_reactions:
            if are_all_substrates_terminal(reaction) == True:
                reactions_to_keep.append(reaction)

        return reactions_to_keep

    def _run_prune(self, network, steps=5, on_substrates=True):
        self._log('Prune network')
        network.evaluator.add_scores_complexity(network)
        nodes_removed = []
        self._log('- delete nodes')
        while len(network.substrate_nodes) >= network.settings['max_nodes']:
            to_remove = self._get_nodes_to_prune(network.graph, steps, on_substrates=on_substrates)

            for node in to_remove:
                if node in list(network.graph.nodes()):
                    deleted_nodes = self.delete_terminal_reaction_node(network, node)
                    nodes_removed.extend(deleted_nodes)
                    if len(network.substrate_nodes) >= network.settings['max_nodes']:
                        break
                else:
                    print('Node not present - could not delete - ' + str(node))
            network.get_node_types()

        return nodes_removed

    def _log(self, msg):
        if self.print_log == True:
            print(msg)

class AIZynthfinder_RuleApplicator(RuleApplicator):

    def __init__(self, network):
        super().__init__(network)
        self.action_applier = aizynthfinder_actions.aizynth_action_applier

    def run(self, smile, graph):
        rxns, metadata = self.action_applier.get_rxns(smile)

        precursor_dict = self.apply_rules(smile, rxns)
        precursor_dict = self._remove_precursors_already_in_graph(graph, smile, precursor_dict)

        if self.network.settings['remove_simple'] == True:
            for name in precursor_dict:
                precursor_dict[name] = self._remove_simple_precursors(precursor_dict[name])

        return precursor_dict, metadata

class RetrosynthesisEngine():

    def __init__(self, network):
        self.network = network
        self.ruleApplication = RuleApplicator(network)
        self.graphManipulator = GraphManipulator(network)
        self.reactionSelector = ReactionSelector(network, GraphPruner())
        self.graphPruner = GraphPruner()
        self.aizynth_rule_application = AIZynthfinder_RuleApplicator(network)

    def single_step(self, smile, rxns, graph, disallowedProducts=None):
        rxn_type = 'retrobiocat'
        rxn_mode = self.network.settings.get('retrobiocat_reaction_mode', 'complexity')

        if self._should_rules_be_applied(smile, graph) == False:
            return [],[]

        rxnsSubstratesDict = self.ruleApplication.run(smile, rxns, graph)
        listProducts, listReactions = self.graphManipulator.add_nodes_to_graph(rxnsSubstratesDict, smile, graph, rxn_type)
        self.network.get_node_types()
        listProducts, listReactions = self.reactionSelector.remove_disallowed_products(disallowedProducts, listProducts,listReactions)
        if rxn_mode == 'complexity':
            listProducts, listReactions = self.reactionSelector.select_best_by_complexity(listProducts, listReactions, self.network.settings['max_reactions'])
        return listProducts, listReactions

    def single_aizynth_step(self, smile, graph, disallowedProducts=None):
        rxn_type = 'aizynth'
        rxn_mode = self.network.settings.get('aizynth_reaction_mode', 'complexity')

        if self._should_rules_be_applied(smile, graph) == False:
            return [],[]

        rxnsSubstratesDict, metadata = self.aizynth_rule_application.run(smile, graph)
        listProducts, listReactions = self.graphManipulator.add_nodes_to_graph(rxnsSubstratesDict, smile, graph, rxn_type, metadata=metadata)
        self.network.get_node_types()
        listProducts, listReactions = self.reactionSelector.remove_disallowed_products(disallowedProducts, listProducts, listReactions)

        if rxn_mode == 'complexity':
            listProducts, listReactions = self.reactionSelector.select_best_by_complexity(listProducts,
                                                                                          listReactions,
                                                                                          self.network.settings['max_reactions'])
        elif rxn_mode == 'policy':
            listProducts, listReactions = self.reactionSelector.select_best_aizynth_by_policy(listProducts,
                                                                                              listReactions,
                                                                                              self.network.settings['max_reactions'])
        return listProducts, listReactions

    def generate_network(self, target_smile, number_steps, rxns, graph, disallowedProducts=None):
        listSmiles = [target_smile]
        for i in range(number_steps):
            newListSmiles, newListReactions = [], []
            for smi in listSmiles:
                newSmiles, newReactions = self.single_step(smi, rxns, graph, disallowedProducts=disallowedProducts)
                newListSmiles.extend(newSmiles)
                newListReactions.extend(newReactions)

            self.network.get_node_types()
            self._log('-- Step ' + str(i + 1) + ' --')
            self._log(str(len(newListReactions)) + ' reactions added, ' + str(len(self.network.substrate_nodes)) + ' substrate nodes')

            listSmiles = self.graphPruner.prune(self.network, newListSmiles)

    def custom_reaction(self, graph, product_smiles, substrate_smiles, reaction_name):
        listSmiles, listReactions = self.graphManipulator.add_custom_reaction(graph, product_smiles, substrate_smiles, reaction_name)
        return listSmiles, listReactions

    def _should_rules_be_applied(self, smile, graph):
        def check_target_is_present_and_substrate(smi, graph):
            """ Returns true is the provided smiles is in the graph as a substrate"""
            if smi not in list(graph.nodes):
                print('WARNING SMILE NOT IN GRAPH ' + str(smi))
                return False
            if graph.nodes[smi]['attributes']['node_type'] != 'substrate':
                print('Warning - target SMILES was not a substrate')
                graph.nodes[smi]['attributes']['node_type'] = 'substrate'
            return True

        def are_substrates_already_in_graph(smi, graph):
            add_if_already_precursors = self.network.settings["add_if_precursor"]

            if (len(list(graph.successors(smi))) == 0) or (add_if_already_precursors == True):
                return True

            return False

        if check_target_is_present_and_substrate(smile, graph) == False:
            return False
        if are_substrates_already_in_graph(smile, graph) == False:
            return False

        return True

    def _log(self, msg):
        if self.network.settings['print_log'] == True:
            print(msg)


if __name__ == '__main__':
    app = RuleApplicator(None)



