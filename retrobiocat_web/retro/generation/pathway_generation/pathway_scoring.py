import pandas as pd
import numpy as np


class PathwayScores():

    def __init__(self, pathway, calc_scores=True, multiple_end_substrates_policy='max'):
        self.pathway = pathway
        self.multiple_end_substrates_policy = multiple_end_substrates_policy

        # scores
        self.change_in_complexity = 0
        self.starting_material = 0
        self.num_intermediates = 0
        self.num_enzyme_steps = 0
        self.complexity_per_enzyme = 0
        self.complexity_per_intermediate = 0
        self.positive_enzyme_steps_score = 0
        self.cofactor_score = 0

        self.simple_score = 0

        if calc_scores==True:
            self.score()

    def scores_dict(self):
        scores_dict = {'change_in_complexity' : self.change_in_complexity,
                       'starting_material' : self.starting_material,
                       'num_intermediates' : self.num_intermediates,
                       'num_enzyme_steps' : self.num_enzyme_steps,
                       'complexity_per_enzyme' : self.complexity_per_enzyme,
                       'complexity_per_intermediate' : self.complexity_per_intermediate,
                       'postitive_enzyme_steps_score' : self.positive_enzyme_steps_score,
                       'cofactor_score' : self.cofactor_score,
                       'simple_score' : self.simple_score}

        return scores_dict

    def scores_list(self):
        scores_list = []
        scores_dict = self.scores_dict()
        for name in scores_dict:
            scores_list.append(scores_dict[name])
        return scores_list

    def scores_from_dict(self, scores_dict):
        self.change_in_complexity = scores_dict['change_in_complexity']
        self.starting_material = scores_dict['starting_material']
        self.num_intermediates = scores_dict['num_intermediates']
        self.num_enzyme_steps = scores_dict['num_enzyme_steps']
        self.complexity_per_enzyme = scores_dict['complexity_per_enzyme']
        self.complexity_per_intermediate = scores_dict['complexity_per_intermediate']
        self.positive_enzyme_steps_score = scores_dict['postitive_enzyme_steps_score']
        self.cofactor_score = scores_dict['cofactor_score']

    def score(self):
        self.starting_material = self.score_starting_material(self.pathway)
        self.num_intermediates = len(self.pathway.substrates) - 2
        self.num_enzyme_steps = self.score_num_enzyme_steps(self.pathway)

        if self.pathway.network.settings['calculate_complexities'] == True:
            self.change_in_complexity = self.score_complexity_change(self.pathway)
            if self.num_enzyme_steps != 0:
                self.complexity_per_enzyme = self.change_in_complexity / self.num_enzyme_steps
            if self.num_intermediates != 0:
                self.complexity_per_intermediate = self.change_in_complexity / self.num_intermediates

        if self.pathway.network.settings['calculate_substrate_specificity'] == True:
            self.positive_enzyme_steps_score = self.score_positive_enzyme_steps(self.pathway, self.num_enzyme_steps)

        self.simple_score = self.calc_simple_score()

    def calc_simple_score(self):
        complexity = self.change_in_complexity
        steps = self.num_enzyme_steps
        enz_steps = self.positive_enzyme_steps_score

        if steps == 0:
            steps = 1

        simple_score = ((complexity/steps) + (complexity*enz_steps)) / 2
        return simple_score

    def score_complexity_change(self, pathway):

        graph = pathway.network.graph
        target_smiles = pathway.network.target_smiles

        end_complexity = []

        if type(pathway.end_nodes) != list:
            self.end_nodes = [self.end_nodes]
            print('Error end nodes are not a list - ' + str(self.end_nodes))

        for node in pathway.end_nodes:
            try:
                end_complexity.append(graph.nodes[node]['attributes']['complexity'])
            except:
                print('Error node not in graph - ' + str(node))
                end_complexity.append(0)

        try:
            if self.multiple_end_substrates_policy == 'max':
                end = np.max(end_complexity)
            elif self.multiple_end_substrates_policy == 'min':
                end = np.min(end_complexity)
            elif self.multiple_end_substrates_policy == 'mean':
                end = np.mean(end_complexity)
            else:
                end = np.mean(end_complexity)
        except:
            print('Error calculating avg complexity change from end nodes ' + str(self.end_nodes))
            end = 0

        start = graph.nodes[target_smiles]['attributes']['complexity']

        return start - end

    def score_starting_material(self, pathway):

        # average of scores
        scores = []
        for node in pathway.end_nodes:
            if pathway.network.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                raise Exception('End node is a reaction')
            else:
                scores.append(pathway.network.graph.nodes[node]['attributes']['is_starting_material'])

        if len(scores) != 0:
            score = sum(scores) / len(scores)
        else:
            score = 0

        return score

    def score_num_enzyme_steps(self, pathway):

        enzyme_steps = 0
        for node in pathway.reactions:
            try:
                enzyme_steps += pathway.network.graph.nodes[node]['attributes']['is_enzyme']
            except:
                print('Warning no is_enzyme attribute - ' + str(type(node)) + str(node))

        return enzyme_steps

    def score_positive_enzyme_steps(self, pathway, num_enzyme_steps):
        total_score = 0
        for node in pathway.reactions:
            enz = pathway.reaction_enzyme_map[node]
            if pathway.network.graph.nodes[node]['attributes']['is_enzyme'] == 1:
                score = pathway.network.graph.nodes[node]['attributes']['specificity_scores'][enz]
                total_score += score

        if num_enzyme_steps == 0:
            positive_enzyme_steps_score = 0
        else:
            positive_enzyme_steps_score = total_score / num_enzyme_steps

        return round(positive_enzyme_steps_score,2)

class PathwayEvaluator():

    def __init__(self, list_pathways):
        self.df = None
        self.pathways = list_pathways
        self.weights = {'Normalised num enzyme steps' : 1,
                        'Normalised change in complexity' : 1,
                        'starting_material' : 0,
                        'postitive_enzyme_steps_score' : 0,
                        'Normalised Cofactor Stoichiometry' : 0}

        self.diversity_weight = 0


    def evaluate(self):
        self.make_evalation_df()
        self.calculate_normalised_scores()
        self.calc_weighted_scores()
        self.df.sort_values('Weighted Score', ascending=False, inplace=True)
        self.df = self.df.reset_index(drop=True)

    def make_evalation_df(self):

        to_df = {'Pathway': [],
                 'Substrate nodes': [],
                 'Reaction nodes': [],
                 'End nodes': []}

        for pathway in self.pathways:
            to_df['Pathway'].append(pathway)
            to_df['Substrate nodes'].append(pathway.substrates)
            to_df['Reaction nodes'].append(pathway.reaction_names)
            to_df['End nodes'].append(pathway.end_nodes)

            scores_dict = pathway.scores.scores_dict()
            for name in scores_dict:
                if name not in to_df:
                    to_df[name] = []
                to_df[name].append(scores_dict[name])

        self.df = pd.DataFrame(to_df)
        return self.df

    def take_best_pathways(self, number):
        self.df = self.df.head(number)

    def calculate_normalised_scores(self):
        def normalise_complexity(df):
            correct_neg_change_in_complexity = 0
            while (df['change_in_complexity'].max() + correct_neg_change_in_complexity) <= 0:
                correct_neg_change_in_complexity += 0.05
            df['Normalised change in complexity'] = df['change_in_complexity'] / (df['change_in_complexity'].max() + correct_neg_change_in_complexity)

            return df

        def normalise_num_steps(df):
            df['Rescaled num enzyme steps'] = df['num_enzyme_steps'].max() - df['num_enzyme_steps']
            df['Normalised num enzyme steps'] = df['Rescaled num enzyme steps'] / df['Rescaled num enzyme steps'].max()
            return df

        def normalise_cofactor(df):
            df['Rescaled Cofactor Stoichiometry'] = df['cofactor_score'].max() - df['cofactor_score']

            if df['Rescaled Cofactor Stoichiometry'].max() != 0:
                df['Normalised Cofactor Stoichiometry'] = df['Rescaled Cofactor Stoichiometry'] / df['Rescaled Cofactor Stoichiometry'].max()
            else:
                df['Normalised Cofactor Stoichiometry'] = df['Rescaled Cofactor Stoichiometry']
            return df

        self.df = normalise_complexity(self.df)
        self.df = normalise_num_steps(self.df)
        self.df = normalise_cofactor(self.df)

    def calc_weighted_scores(self):

        scores_to_df = []

        for index, row in self.df.iterrows():
            score = 0
            for name in self.weights:
                score += row[name] * self.weights[name]
            scores_to_df.append(score)

        self.df['Weighted Score'] = scores_to_df

        if self.diversity_weight != 0:
            self.rank_diversity()

        return self.df

    def apply_diversity(self, top=25, reaction_only=False):

        diversity = []
        div_dict = {}

        for index, row in self.df.iterrows():

            if index <= top:
                pathway = row['Pathway']
                pathway_diversity = 0

                for node in pathway.list_nodes:
                    if (reaction_only == False) or pathway.network.graph.nodes[node]['attributes']['node_type'] == 'reaction':
                        node_name = pathway.network.graph.nodes[node]['attributes']['name']
                        if node_name not in div_dict:
                            div_dict[node_name] = 0
                        pathway_diversity += div_dict[node_name]
                        div_dict[node_name] += 1
                diversity.append(pathway_diversity)
            else:
                diversity.append(0)

        self.df['Diversity'] = diversity
        self.df['Diversity'] = self.df['Diversity'] / self.df['Diversity'].max()

        self.df['Diversity Weighted Score'] = (self.df['Diversity']*self.diversity_weight) + self.df['Weighted Score']


    def rank_diversity(self, n=10):

        for i in range(n):
            self.apply_diversity()
            self.df.sort_values('Diversity Weighted Score', ascending=False, inplace=True)