from retrobiocat_web.mongo.models.biocatdb_models import Activity, Paper
from retrobiocat_web.mongo.models.reaction_models import AutoProcessingRule, Reaction
import pandas as pd
from retrobiocat_web.retro.generation.network_generation.network import Network
from retrobiocat_web.retro.generation.pathway_generation.best_first_search import BFS
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.mongo.init_db import biocatdb_excel_to_db
from retrobiocat_web.mongo.init_db.data_checks import check_is_float, check_is_nan, get_smile, get_mol
import mongoengine as db


COLS = ['enzyme_type', 'enzyme_name', 'reaction', 'short_citation', 'html_doi', 'added_by_string',
                'substrate_1_smiles', 'substrate_2_smiles', 'product_1_smiles', 'temperature', 'ph',
                'solvent', 'other_conditions', 'notes', 'formulation', 'biocat_conc', 'reaction_vol', 'substrate_1_conc',
                'substrate_2_conc', 'conversion', 'conversion_time', 'categorical', 'binary', 'auto_generated',
                'reviewed', 'selectivity', 'kcat', 'km', 'mw', 'specific_activity']


class AutoProcessor():

    def __init__(self):
        self.converters = []
        for rule in AutoProcessingRule.objects().select_related():
            self.converters.append(Converter(rule))

    def get_activity(self, paper):
        acts = Activity.objects(paper=paper).as_pymongo().only(*COLS)
        spec_df = pd.DataFrame(list(acts))
        spec_df.drop(columns=['_id'], inplace=True)
        return spec_df

    def delete_old_autoprocessed(self, paper):
        pq = db.Q(paper=paper)
        ap = db.Q(auto_generated=True)

        acts = Activity.objects(pq & ap)
        for act in acts:
            act.delete()

    def get_only_new_rows(self, df, new_df):
        new_df = pd.merge(df, new_df, indicator=True, how='outer')
        new_df = new_df[new_df['_merge'] == 'right_only']
        new_df.drop(columns=['_merge'], inplace=True)
        return new_df

    def add_to_db(self, new_df, paper):
        if len(new_df) != 0:
            self.log(f"..adding {len(new_df)} new activity rows")

        for i, row in new_df.iterrows():
            if row['binary'] == 1:
                binary = True
            else:
                binary = False

            if row['auto_generated'] == 1:
                auto_gen = True
            else:
                auto_gen = False

            activity = Activity(enzyme_type=check_is_nan(row['enzyme_type']),
                                enzyme_name=check_is_nan(row['enzyme_name']),
                                reaction=check_is_nan(row['reaction']),
                                short_citation=check_is_nan(row['short_citation']),
                                html_doi=check_is_nan(row['html_doi']),
                                added_by_string=str(row['added_by_string']),
                                paper=paper,
                                substrate_1_smiles=get_smile(row['substrate_1_smiles']),
                                substrate_2_smiles=get_smile(row['substrate_2_smiles']),
                                product_1_smiles=get_smile(row['product_1_smiles']),
                                temperature=check_is_nan(row['temperature']),
                                ph=check_is_nan(row['ph']),
                                solvent=check_is_nan(row['solvent']),
                                other_conditions=check_is_nan(row['other_conditions']),
                                notes=check_is_nan(row['notes']),
                                reaction_vol=check_is_nan(row['reaction_vol']),
                                formulation=check_is_nan(row['formulation']),
                                biocat_conc=check_is_nan(row['biocat_conc']),
                                kcat=check_is_float(row['kcat']),
                                km=check_is_float(row['km']),
                                mw=check_is_float(row['mw']),
                                substrate_1_conc=check_is_nan(row['substrate_1_conc']),
                                substrate_2_conc=check_is_nan(row['substrate_2_conc']),
                                specific_activity=check_is_float(row['specific_activity']),
                                conversion=check_is_float(row['conversion']),
                                conversion_time=check_is_float(row['conversion_time']),
                                categorical=check_is_nan(row['categorical']),
                                binary=binary,
                                selectivity=check_is_nan(row['selectivity']),
                                auto_generated=auto_gen)
            activity.save()

    def auto_process(self, paper):

        self.log(f"running paper - {paper.short_citation}")

        # remove old autoprocessed data for this paper
        self.delete_old_autoprocessed(paper)

        df = self.get_activity(paper)
        list_new_dfs = []
        for converter in self.converters:
            new_df = converter.apply_rule(df)
            list_new_dfs.append(new_df)

        final_new_df = pd.concat(list_new_dfs)
        final_new_df = self.get_only_new_rows(df, final_new_df)

        self.add_to_db(final_new_df, paper)

    def log(self, msg):
        print(f"AUTOPROCESSOR: {msg}")

class Converter(object):
    def __init__(self, rule):
        self.rule = rule

    def apply_rule(self, df):
        new_df = self._reaction_conversion(df,
                                           self.rule.multi_step_reaction,
                                           self.rule.reactions,
                                           self.rule.min_steps,
                                           self.rule.max_steps,
                                           self.rule.ignore_substrate_two)
        new_df = new_df.drop_duplicates()
        return new_df

    def _reaction_conversion(self, df, reaction_name, part_reaction_names, min_steps, max_steps, ignore_substrate_two, log=True):
        products_pathways_dict = {}
        list_new_dfs = []
        for index, row in df.iterrows():
            if row['reaction'].lower() == reaction_name.lower():
                new_row = self._duplicate_row(row)

                product = row['product_1_smiles']
                if product not in products_pathways_dict:
                    pathways = self._get_retrobiocat_pathways(new_row, part_reaction_names)
                    pathways = self._get_pathways_to_substrates(new_row, pathways, ignore_substrate_two)

                    pathways = self._check_expected_steps(pathways, min_steps, max_steps)
                    if log == True:
                        print('Auto data extraction for ' + str(reaction_name) + ': ' + str(
                            len(pathways)) + ' route to ' + str(row['product_1_smiles']))
                    products_pathways_dict[product] = pathways

                new_df = self._get_new_df_from_pathways(df, new_row, products_pathways_dict[product])
                list_new_dfs.append(new_df)

        list_new_dfs.append(pd.DataFrame(columns=COLS))
        if len(list_new_dfs) != 0:
            new_df = pd.concat(list_new_dfs)
            new_df = self._select_only_positive_binary(new_df)
            new_df = self._select_only_single_positive_chemical_step(new_df)

        return new_df

    def _duplicate_row(self, row):
        def set_added_by(row):
            row['added_by_string'] = 'AUTOMATED ENTRY'
            row['notes'] = 'AUTOMATED ENTRY'
            return row

        new_row = row.copy()
        new_row = set_added_by(new_row)
        return new_row

    def _get_retrobiocat_pathways(self, row, reaction_names, combine_enantiomers=True):
        product = row['product_1_smiles']
        network = Network(print_log=False)
        network.settings.update({"calculate_complexities": True,
                                 "calculate_substrate_specificity": False,
                                 "get_building_blocks": False,
                                 "combine_enantiomers": combine_enantiomers})

        # Only want reactions in list
        rxns_to_keep = {}
        for rxn_name in network.rxns:
            if rxn_name in reaction_names:
                rxns_to_keep[rxn_name] = network.rxns[rxn_name]
        network.rxns = rxns_to_keep

        network.generate(product, 5)
        bfs = BFS(network=network, print_log=False, score_pathways=False)
        bfs.run()
        pathways = bfs.get_pathways()

        return pathways

    def _get_pathways_to_substrates(self, row, pathways, ignore_substrate_two):
        substrate_one = rdkit_smile(row['substrate_1_smiles'])
        substrate_two = rdkit_smile(row['substrate_2_smiles'])

        pathways_to_keep = []
        for pathway in pathways:
            if substrate_one in pathway.end_nodes:
                if ignore_substrate_two == True:
                    pathways_to_keep.append(pathway)
                elif substrate_two != None:
                    if substrate_two in pathway.end_nodes:
                        pathways_to_keep.append(pathway)

        return pathways_to_keep

    def _get_new_df_from_pathways(self, df, row, pathways):
        new_df_dict = {}

        for col in list(df.columns):
            new_df_dict[col] = []

        i = 0
        for pathway in pathways:
            for reaction in pathway.reactions:

                for col in list(df.columns):
                    new_df_dict[col].append(row[col])

                new_df_dict['auto_generated'][i] = True

                reaction_name = pathway.network.graph.nodes[reaction]['attributes']['name']
                new_df_dict['reaction'][i] = reaction_name

                reaction_products = list(pathway.network.graph.predecessors(reaction))
                assert len(reaction_products) == 1
                new_df_dict['product_1_smiles'][i] = reaction_products[0]

                reaction_substrates = list(pathway.network.graph.successors(reaction))
                new_df_dict['substrate_1_smiles'][i] = reaction_substrates[0]
                if len(reaction_substrates) == 2:
                    new_df_dict['substrate_2_smiles'][i] = reaction_substrates[1]
                else:
                    new_df_dict['substrate_2_smiles'][i] = ''
                i += 1

        new_df = pd.DataFrame(new_df_dict)
        new_df = new_df.drop_duplicates()

        return new_df

    def _check_expected_steps(self, pathways, min_steps, max_steps):
        pathways_to_keep = []
        for pathway in pathways:
            if len(pathway.reactions) > max_steps:
                print('Warning - pathway to target ' + str(pathway.target_smiles) + ' more than max steps ' + str(
                    max_steps))
            elif len(pathway.reactions) < min_steps:
                print('Warning - pathway to target ' + str(pathway.target_smiles) + ' less than min steps ' + str(
                    min_steps))
            else:
                pathways_to_keep.append(pathway)
        return pathways_to_keep

    def _select_only_single_positive_chemical_step(self, df):

        network = Network()
        rxn_obj = network.rxn_obj
        rxn_obj.load_additional_info()

        # Change enzyme name to just 'Chemical' for chemical steps

        enzymes = []
        for i, row in df.iterrows():
            name = row['reaction']
            if name in rxn_obj.reactions:
                if name in rxn_obj.rules_by_type['Chemical']:
                    enzymes.append('Chemical')
                else:
                    enzymes.append(row['enzyme_name'])
            else:
                enzymes.append(row['enzyme_name'])

        df['enzyme_name'] = enzymes

        # Filter out negative binary data which is chemical
        # df = df[(df['Binary'] != 0) & (df['Enzyme name'] == 'Chemical') | (df['Enzyme name'] != 'Chemical')]

        # Remove duplicate entries
        df = df.drop_duplicates(
            ['reaction', 'enzyme_name', 'product_1_smiles', 'substrate_1_smiles', 'substrate_2_smiles'])

        return df

    def _select_only_positive_binary(self, df):
        df = df[(df['binary'] == 1)]
        return df

def task_autoprocess():
    ap = AutoProcessor()
    paper_ids = list(Paper.objects(status='Complete').distinct('_id'))
    for id in paper_ids:
        paper = Paper.objects(id=id)[0]
        try:
            ap.auto_process(paper)
        except Exception as e:
            print()
            print(f'WARNING - ERROR PROCESSING PAPER - {paper.short_citation}')
            print()
            print(e)
            print()



if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    task_autoprocess()

    #AutoProcessingRule.drop_collection()
    #new_rule = AutoProcessingRule(multi_step_reaction='Reductive amination', reactions=['Aldehyde amination', 'Ketone amination'], min_steps=1, max_steps=1, ignore_substrate_two=True)
    #new_rule.save()

    #paper = Paper.objects(short_citation='Montgomery et al, 2017, Angew. Chem. Int. Ed.')[0]
    #print(paper)

    #ap = AutoProcessor()
    #ap.auto_process(paper)

    #for paper in Paper.objects():
    #    ap.auto_process(paper)


