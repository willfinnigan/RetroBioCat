from pathlib import Path
import pandas as pd
import time
from retrobiocat_web.retro.generation.network_generation.network import Network
from retrobiocat_web.retro.generation.pathway_generation.best_first_search import BFS
from retrobiocat_web.retro.generation.pathway_generation.group_pathways import group_pathways
from retrobiocat_web.retro.generation.pathway_generation.pathway_scoring import PathwayEvaluator
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile, disable_rdkit_logging
from retrobiocat_web import __version__


class TestSetEvaluator(object):

    def __init__(self, excel=None, log=False,
                 weights=None,
                 num_steps=4, max_nodes=400, max_pathways=25000,
                 min_weight=1,
                 only_pos=False, spec_threshold=0.6):
        self.print_log = log

        if excel is not None:
            self.df = self.load_df(excel)

        if weights is not None:
            self.weights = weights
        else:
            self.weights = [1,1,1,1,1]

        self.num_steps = num_steps
        self.max_nodes = max_nodes
        self.max_pathways = max_pathways
        self.min_weight = min_weight
        self.only_positive_specificity = only_pos
        self.specificity_threshold = spec_threshold

        self.results = {}

    def load_df(self, excel):
        df = pd.read_excel(excel, index_col='Num')
        df = self._split_substrates_to_list(df)
        df = self._split_enzymes_to_list(df)
        return df

    def create_pathways_for_test(self, target_product):
        network = Network()

        network.settings.update({"combine_enantiomers": True,
                                 "remove_simple": True,
                                 "similarity_score_threshold": self.specificity_threshold,
                                 "num_enzymes": 1,
                                 "max_nodes": self.max_nodes,
                                 "prune_steps": 1,
                                 "only_postitive_enzyme_data": self.only_positive_specificity,
                                 'prune_on_substrates': False,
                                 'max_reactions': False,
                                 'include_experimental': False,
                                 'include_two_step': False})

        network.generate(target_product, self.num_steps)

        bfs = BFS(network=network, max_pathways=self.max_pathways, min_weight=self.min_weight)
        bfs.run()
        pathways = bfs.get_pathways()
        if len(pathways) == self.max_pathways:
            self.log(f"  ~max pathways reached")

        #pathways = group_pathways(pathways)
        return pathways

    def evaluate_pathways(self, pathways):
        pathway_evaluator = PathwayEvaluator(pathways)
        pathway_evaluator.weights = {'Normalised num enzyme steps': self.weights[0],
                                     'Normalised change in complexity': self.weights[1],
                                     'starting_material': self.weights[2],
                                     'postitive_enzyme_steps_score': self.weights[3],
                                     'Normalised Cofactor Stoichiometry': 0}
        pathway_evaluator.diversity_weight = self.weights[4]

        pathway_evaluator.evaluate()

        return pathway_evaluator

    def identify_locations_of_test_substrates(self, df, list_test_substrates, list_test_enzymes):
        list_indexes_abs = []
        list_rows_abs = []
        list_indexes_any = []
        list_rows_any = []


        for index, row in df.iterrows():

            # For absolute ranking
            if self._are_literature_end_nodes_present(list_test_substrates, row['End nodes']):
                enzymes = self._get_reaction_enzymes(row['Pathway'])
                if (all(elem in enzymes for elem in list_test_enzymes)
                    and len(list_test_enzymes) == row['Pathway'].scores.num_enzyme_steps) \
                        or len(list_test_enzymes) == 0:
                    self.log(f"  ~ absolute ranking {index}, and enzymes match")
                    list_indexes_abs.append(index)
                    list_rows_abs.append(row)
                else:
                    self.log(f"  ~ absolute ranking {index}, but enzymes do not match")

            # For any ranking
            if self._are_literature_end_nodes_present(list_test_substrates, row['Substrate nodes']):
                enzymes = self._get_reaction_enzymes(row['Pathway'])
                if all(elem in enzymes for elem in list_test_enzymes) or len(list_test_enzymes) == 0:
                    list_indexes_any.append(index)
                    list_rows_any.append(row)
                    self.log(f"  ~ any ranking {index}, and enzymes match")
                else:
                    self.log(f"  ~ any ranking {index}, but enzymes do not match")

        return list_indexes_any, list_rows_any, list_indexes_abs, list_rows_abs

    def run(self, id='default'):
        self.log(f"Running (ID = {id})")
        absolute_rank = []
        any_rank = []
        times = []

        for index, row in self.df.iterrows():
            self.log(f"Testing {row}")

            t0 = time.time()
            self.log(f" - create pathways for {row['Target Product']}")
            pathways = self.create_pathways_for_test(row['Target Product'])

            self.log(f" - evaluate pathways for {row['Target Product']}")
            evaluator = self.evaluate_pathways(pathways)

            self.log(f" - find ranking of literature pathway")
            lit_end_nodes = self._lit_end_nodes_as_rdkit_smiles(row['Substrates_list'])
            indexes_any, rows_any, indexes_abs, rows_abs = self.identify_locations_of_test_substrates(evaluator.df,
                                                                                                      lit_end_nodes,
                                                                                                      row['Enzymes_list'])
            t1 = time.time()
            t = round(t1 - t0, 0)
            times.append(t)
            self.log(f'Complete - runtime of {t} seconds.')

            if len(indexes_abs) != 0:
                absolute_rank.append(indexes_abs[0] + 1)
            else:
                absolute_rank.append('')

            if len(indexes_any) != 0:
                any_rank.append(indexes_any[0] + 1)
            else:
                any_rank.append('')

        self.results[id] = (absolute_rank, any_rank, times)

        return absolute_rank, any_rank, times

    def save(self):
        for id in self.results:
            self.df['v' + __version__ + '_' + id + '_any_rank'] = self.results[id][1]
            self.df['v' + __version__ + '_' + id + '_absolute_rank'] = self.results[id][0]
            self.df['v' + __version__ + '_' + id + '_time (s)'] = self.results[id][2]

        self.df.drop(columns=['Substrates_list', 'Enzymes_list'], inplace=True)

        self.df.to_excel(path_to_excel, index_label='Num')

    def log(self, msg):
        if self.print_log == True:
            print(msg)


    def _split_substrates_to_list(self, df):
        substrates_as_list = []

        for index, row in df.iterrows():
            if isinstance(row['Substrates'], str):
                split_list = row['Substrates'].split(", ")
                rdkit_split_list = []
                for smi in split_list:
                    rdkit_split_list.append(rdkit_smile(smi))
                substrates_as_list.append(rdkit_split_list)
            else:
                substrates_as_list.append([])

        df['Substrates_list'] = substrates_as_list
        return df

    def _split_enzymes_to_list(self, df):
        enzymes_as_list = []

        for index, row in df.iterrows():
            if isinstance(row['Enzymes'], str):
                split_list = row['Enzymes'].split(", ")
                enzymes_as_list.append(split_list)
            else:
                enzymes_as_list.append([])

        df['Enzymes_list'] = enzymes_as_list
        return df

    def _get_reaction_enzymes(self, pathway):
        enzymes = []
        for reaction in pathway.reactions:
            potential_enzymes = pathway.network.graph.nodes[reaction]['attributes']['possible_enzymes']
            enzymes.extend(potential_enzymes)
        return enzymes

    def _are_literature_end_nodes_present(self, lit_end_nodes, pathway_nodes):
        for node in lit_end_nodes:
            if node not in pathway_nodes:
                return False
        return True

    def _lit_end_nodes_as_rdkit_smiles(self, lit_end_nodes):
        rdkit_nodes = []
        for node in lit_end_nodes:
            rdkit_nodes.append(rdkit_smile(node))
        return rdkit_nodes


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    disable_rdkit_logging()

    path_to_excel = 'test_pathways.xlsx'
    test_set_eval = TestSetEvaluator(excel=path_to_excel,
                                     num_steps=4, max_nodes=400, max_pathways=40000,
                                     min_weight=1, only_pos=False, spec_threshold=0.6,
                                     log=True, weights=[1,1,1,0,1])
    test_set_eval.run(id='venv_no_lit')
    test_set_eval.save()
