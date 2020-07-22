
from retrobiocat_web.retro.generation.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine, RuleApplicator, ReactionSelector, GraphManipulator
from pathlib import Path
import pandas as pd
import pickle
import sqlite3
import uuid

from rdkit.Chem import AllChem, rdmolops
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.retro.generation import node_analysis

PATH_TO_RETRORULES_FOLDER = str(Path(__file__).parents[3]) + '/data/retrorules'

class RetroRulesRuleApplicator(RuleApplicator):

    def __init__(self, network):
        super().__init__(network)
        self.path_to_retrorules_folder = PATH_TO_RETRORULES_FOLDER
        self.cofactors_df = self.load_cofactors()

    def run(self, smile, rxns, graph, removeRetroRuleCofactors=True, explicit_hydrogens=False):
        precursor_dict = self.apply_retrorules(smile, rxns, explicit_hydrogens=explicit_hydrogens)
        precursor_dict = self._remove_precursors_already_in_graph(graph, smile, precursor_dict)

        if self.network.settings['remove_simple'] == True:
            for name in precursor_dict:
                precursor_dict[name] = self._remove_simple_precursors(precursor_dict[name])

        if removeRetroRuleCofactors == True:
            precursor_dict = self.remove_cofactors(precursor_dict)

        return precursor_dict

    def apply_retrorules(self, smile, rxns, explicit_hydrogens=False):
        '''Function takes a smile and dictionary of reactions, applys the reactions and
           returns a dictionary of rxn_names : products '''
        try:
            substrate_molecule = AllChem.MolFromSmiles(smile)
        except:
            return {}

        if explicit_hydrogens == True:
            substrate_molecule = rdmolops.AddHs(substrate_molecule)

        rxn_product_dict = {}
        for rxn_name, rxn in rxns.items():
            try:
                products = rxn.RunReactants((substrate_molecule,))
            except:
                products = []
                print('Error running reactants for: ' + str(smile))

            smiles_products = []
            for product in products:
                sub_list = []
                for mol in product:
                    mols = [mol]

                    if explicit_hydrogens == True:
                        mol = rdmolops.RemoveHs(mol)

                    try:
                        mols = rdmolops.GetMolFrags(mol, asMols=True)
                    except:
                        pass

                    for mol in mols:
                        try:
                            p_smile = AllChem.MolToSmiles(mol)
                            p_smile = rdkit_smile(p_smile)
                            if self._check_valid_smile(p_smile, rxn_name=rxn_name) == True:
                                sub_list.append(p_smile)
                        except:
                            pass

                if (sub_list not in smiles_products) and (len(sub_list) != 0):
                    smiles_products.append(sub_list)

            if len(smiles_products) != 0:
                rxn_product_dict[rxn_name] = smiles_products

        return rxn_product_dict

    def load_cofactors(self):
        path = self.path_to_retrorules_folder + '/retrorules_cofactors.tsv'
        df = pd.read_csv(path, sep='\t', names=['inchi', 'name', 'reactions'])
        return df

    def remove_cofactors(self, products_dict):
        cofactors = {}
        for id in products_dict:
            for i, products in enumerate(products_dict[id]):
                for smi in products:
                    inchi = retrobiocat.retro.generation.node_analysis.get_inchl(smi)
                    if inchi in list(self.cofactors_df['inchi']):
                        products_dict[id][i].remove(smi)
                        cofactors[id] = self.cofactors_df['name']
        return products_dict

class RetroRulesReactionSelector(ReactionSelector):

    def __init__(self, network, graphPruner):
        super().__init__(network, graphPruner)


    def sort_reactions_into_groups(self, listReactions, listProducts):
        groups = {}
        for reaction in listReactions:
            reactionProducts = list(self.network.graph.successors(reaction))
            reactionProducts.sort()
            reactionProducts = str(reactionProducts)
            if reactionProducts not in groups:
                groups[reactionProducts] = []
            groups[reactionProducts].append(reaction)
        return groups

    def pick_best_reaction_per_group(self, groups, max_reactions):
        groupBestReactions = {}

        for listProducts in groups:
            listReactions = groups[listProducts]
            best = self._select_best(listReactions, self.network.settings['rr_max_reactions'])
            groupBestReactions[listProducts] = best

        return groupBestReactions

    def select_best_reactions(self, listReactions, listProducts, max_reactions):
        groups = self.sort_reactions_into_groups(listReactions, listProducts)
        best_groups = self.pick_best_reaction_per_group(groups, max_reactions)

        to_keep = []
        for groupReactions in best_groups.values():
            to_keep.extend(groupReactions)

        for reaction in listReactions:
            if reaction not in to_keep:
                self.network.graph.remove_node(reaction)
                self.network.reaction_nodes.remove(reaction)

        return listProducts, to_keep

    def _get_product_reactions(self, product, listReactions):
        productReactions = []
        for reaction in list(self.network.graph.predecessors(product)):
            if reaction in listReactions:
                productReactions.append(reaction)
        return productReactions

    def _get_reaction_scores(self, listReactions):
        productReactionScores = []
        for reaction in listReactions:
            rule_name = self.network.graph.nodes[reaction]['attributes']['rule_id']
            score = self.get_retrorules_score(rule_name)
            productReactionScores.append(score)
        return productReactionScores

    def _select_best(self, listReactions, max_reactions):
        listReactionsToKeep = []
        reactionScores = self._get_reaction_scores(listReactions)
        sorted_reactions = node_analysis.sort_by_score(listReactions, reactionScores, reverse=False)
        listReactionsToKeep.extend(sorted_reactions[0:max_reactions])

        return listReactionsToKeep


    def get_retrorules_score(self, rule_name):
        conn = self.network.retrorules.retrorule_db

        c = conn.cursor()
        query = "SELECT [# Rule_ID], [Score_normalized] FROM RetroRules WHERE [# Rule_ID] = '"+ rule_name + "';"
        c.execute(query)
        result = c.fetchall()
        if len(result) != 0:
            score = result[0][1]
        else:
            print('Warning - no retrorules score found for reaction: ' + str(rule_name))
            score = 0
        return score

class RetroRulesGraphManipulator(GraphManipulator):

    def __init__(self, network):
        super().__init__(network)

    def get_retrorules_name(self, reaction_name, conn='default'):
        if conn == 'default':
            conn = self.network.retrorules.retrorule_db

        c = conn.cursor()

        query = "SELECT [# Rule_ID], [Diameter], [Reaction_ID] FROM RetroRules WHERE [# Rule_ID] = '"+ reaction_name + "';"
        c.execute(query)
        result = c.fetchall()
        if len(result) != 0:
            name = result[0][2]
            diameter = result[0][1]
            combined_name = str(diameter) + '_' + str(name)
        else:
            print('Warning - no retrorules name found for reaction: ' + str(reaction_name))
            combined_name = reaction_name

        return combined_name


    def _add_reaction_node_to_graph(self, graph, reaction_name, source_smile):
        unique_reaction_name = reaction_name + '(' + str(uuid.uuid1()) + ')'
        name = self.get_retrorules_name(reaction_name)

        graph.add_node(unique_reaction_name, attributes={'name': name,
                                                         'rule_id' : reaction_name,
                                                         'node_type': 'reaction',
                                                         'node_num': self._get_node_number(unique_reaction_name, graph),
                                                         'retrorule' : True})

        graph.add_edge(source_smile, unique_reaction_name)
        return unique_reaction_name

class RetroRulesRetrosynthesisEngine(RetrosynthesisEngine):

    def __init__(self, network, retrorules):
        super().__init__(network)
        self.reactionSelector = RetroRulesReactionSelector(network, self.graphPruner)
        self.ruleApplication = RetroRulesRuleApplicator(network)
        self.retrorules = retrorules
        self.graphManipulator = RetroRulesGraphManipulator(network)


    def single_step(self, smile, rxns, graph, disallowedProducts=None):
        if self._should_rules_be_applied(smile, graph) == False:
            return [],[]

        if disallowedProducts == None:
            disallowedProducts = []

        listProducts, listReactions = [], []
        diameter = 16

        while True:
            if diameter < self.network.settings['rr_min_diameter']:
                break
            if len(listProducts) > self.network.settings['rr_min_products']:
                break

            rxns_to_apply = rxns[diameter]
            diameter = diameter-2

            rxnsSubstratesDict = self.ruleApplication.run(smile, rxns_to_apply, graph, removeRetroRuleCofactors=True)

            newProducts, newReactions = self.graphManipulator.add_nodes_to_graph(rxnsSubstratesDict, smile, graph)

            newProducts = list(dict.fromkeys(newProducts))

            newProducts, newReactions = self.reactionSelector.remove_disallowed_products(disallowedProducts, newProducts,newReactions)

            newProducts, newReactions = self.reactionSelector.select_best_reactions(newReactions, newProducts,
                                                                                    max_reactions=self.network.settings['rr_max_reactions'])

            newProducts, newReactions = self.reactionSelector.select_best_by_complexity(newProducts, newReactions, self.network.settings['max_reactions'])

            disallowedProducts.extend(newProducts)
            listProducts.extend(newProducts)
            listReactions.extend(newReactions)

        return listProducts, listReactions

class RetroRules():

    def __init__(self, network):
        self.network = network

        self.retrosynthesisEngine = RetroRulesRetrosynthesisEngine(network, self)
        self.retrorules_rxns = None
        self.retrorule_db = None
        self.path_to_retrorules_folder = PATH_TO_RETRORULES_FOLDER

        self.diameters = retrobiocat.retrorules_diameters

    def load(self, pkl_file='', retrorule_db=''):
        print('Loading retrorules')
        if pkl_file == '':
            self.retrorules_rxns = self._load_retrorule_pkls()
        else:
            self.retrorules_rxns = pickle.load(open(pkl_file, "rb"))

        if retrorule_db == '':
            retrorule_db = self.path_to_retrorules_folder + '/retrorules.db'

        self.retrorule_db = sqlite3.connect(retrorule_db)

    def _load_retrorule_pkls(self):
        path = self.path_to_retrorules_folder + '/'

        rxns = {}
        all_d = [2,4,6,8,10,12,14,16]
        for d in all_d:
            if d in self.diameters:
                file = 'rules' + str(d) + '.pkl'
                rxns[d] = pickle.load(open(path + file, "rb"))
            else:
                rxns[d] = {}

        return rxns

    def add_step(self, smiles, calculate_scores=True):
        """
        Add a single retrosynthetic step to graph from a single smiles node
        """

        smiles = node_analysis.rdkit_smile(smiles)

        listSmiles, listReactions = self.retrosynthesisEngine.single_step(smiles, self.retrorules_rxns, self.network.graph)
        self.network.get_node_types()

        if calculate_scores == True:
            self.network.evaluator.calculate_scores(self.network)

        return listSmiles, listReactions


if __name__ == '__main__':
    from retrobiocat_web.retro.generation.network_generation.network import Network
    target = 'CCCCCO'
    network = Network()
    network.generate(target, 2)
    network.retrorules.diameters=[2]
    network.retrorules.load()
    network.retrorules.add_step('CCCCCC(C)=O')

    """
    file = str(Path(__file__).parents[3]) + '/data/reaction_rules/retrorules/retrorules_all.pkl'
    rxns = pickle.load(open(file, "rb"))

    for d in rxns:
        print(d)
        file_name = 'rules' + str(d) + '.pkl'
        with open(file_name, 'wb') as handle:
            pickle.dump(rxns[d], handle, protocol=pickle.HIGHEST_PROTOCOL)
    """