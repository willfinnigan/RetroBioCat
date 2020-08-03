import pandas as pd
from retrobiocat_web.retro.generation.network_generation.network import Network
from retrobiocat_web.retro.generation.pathway_generation.best_first_search import BFS
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile


def duplicate_row(row):
    def set_added_by(row):
        row['Added by'] = 'AUTOMATED ENTRY FROM REDUCTIVE AMINATION - ' + str(row['Added by'])
        row['Notes'] = 'AUTOMATED ENTRY FROM REDUCTIVE AMINATION - ' + str(row['Notes'])
        return row

    new_row = {}
    for col in COLUMNS:
        new_row[col] = row[col]

    new_row = set_added_by(new_row)
    return new_row

def get_retrobiocat_pathways(row, reaction_names, combine_enantiomers=True):
    product = row['Product 1 SMILES']
    network = Network(print_log=False)
    network.settings.update({"calculate_complexities" : True,
                             "calculate_substrate_specificity" : False,
                             "get_building_blocks" : False,
                             "combine_enantiomers" : combine_enantiomers})

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

def get_pathways_to_substrates(row, pathways, ignore_substrate_two):
    substrate_one = rdkit_smile(row['Substrate 1 SMILES'])
    substrate_two = rdkit_smile(row['Substrate 2 SMILES'])

    pathways_to_keep = []
    for pathway in pathways:
        if substrate_one in pathway.end_nodes:
            if ignore_substrate_two == True:
                pathways_to_keep.append(pathway)
            elif substrate_two != None:
                if substrate_two in pathway.end_nodes:
                    pathways_to_keep.append(pathway)

    return pathways_to_keep

def get_new_df_from_pathways(row, pathways):
    new_df_dict = {}

    for col in COLUMNS:
        new_df_dict[col] = []

    i = 0
    for pathway in pathways:
        for reaction in pathway.reactions:

            for col in COLUMNS:
                new_df_dict[col].append(row[col])

            new_df_dict['Auto Generated'][i] = True

            reaction_name = pathway.network.graph.nodes[reaction]['attributes']['name']
            new_df_dict['Reaction'][i] = reaction_name

            reaction_products = list(pathway.network.graph.predecessors(reaction))
            assert len(reaction_products) == 1
            new_df_dict['Product 1 SMILES'][i] = reaction_products[0]

            reaction_substrates = list(pathway.network.graph.successors(reaction))
            new_df_dict['Substrate 1 SMILES'][i] = reaction_substrates[0]
            if len(reaction_substrates) == 2:
                new_df_dict['Substrate 2 SMILES'][i] = reaction_substrates[1]
            else:
                new_df_dict['Substrate 2 SMILES'][i] = ''
            i += 1

    new_df = pd.DataFrame(new_df_dict)
    new_df = new_df.drop_duplicates()

    return new_df

def check_expected_steps(pathways, min_steps, max_steps):
    pathways_to_keep = []
    for pathway in pathways:
        if len(pathway.reactions) > max_steps:
            print('Warning - pathway to target ' + str(pathway.target_smiles) + ' more than max steps ' + str(max_steps))
        elif len(pathway.reactions) < min_steps:
            print('Warning - pathway to target ' + str(pathway.target_smiles) + ' less than min steps ' + str(min_steps))
        else:
            pathways_to_keep.append(pathway)
    return pathways_to_keep

def select_only_single_positive_chemical_step(df):

    network = Network()
    rxn_obj = network.rxn_obj
    rxn_obj.load_additional_info()

    # Change enzyme name to just 'Chemical' for chemical steps

    enzymes = []
    for i, row in df.iterrows():
        name = row['Reaction']
        if name in rxn_obj.reactions:
            if name in rxn_obj.rules_by_type['Chemical']:
                enzymes.append('Chemical')
            else:
                enzymes.append(row['Enzyme name'])
        else:
            enzymes.append(row['Enzyme name'])

    df['Enzyme name'] = enzymes

    # Filter out negative binary data which is chemical
    #df = df[(df['Binary'] != 0) & (df['Enzyme name'] == 'Chemical') | (df['Enzyme name'] != 'Chemical')]

    # Remove duplicate entries
    df = df.drop_duplicates(['Reaction', 'Enzyme name', 'Product 1 SMILES', 'Substrate 1 SMILES', 'Substrate 2 SMILES'])

    return df

def select_only_positive_binary(df):
    df = df[(df['Binary'] == 1)]
    return df


def reaction_coversion(df, df_reaction, retrobiocat_reactions,
                       min_steps=2, max_steps=2, ignore_substrate_two=False,
                       log=True):

    products_pathways_dict = {}
    list_new_dfs = []
    for index, row in df.iterrows():
        if row['Reaction'].lower() == df_reaction.lower():
            new_row = duplicate_row(row)
            product = row['Product 1 SMILES']
            if product not in products_pathways_dict:
                pathways = get_retrobiocat_pathways(new_row, retrobiocat_reactions)
                pathways = get_pathways_to_substrates(new_row, pathways, ignore_substrate_two)

                pathways = check_expected_steps(pathways, min_steps, max_steps)
                if log == True:
                    print('Auto data extraction for ' + str(df_reaction) +': ' + str(len(pathways)) + ' route to ' + str(row['Product 1 SMILES']))
                products_pathways_dict[product] = pathways

            new_df = get_new_df_from_pathways(new_row, products_pathways_dict[product])
            list_new_dfs.append(new_df)

    if len(list_new_dfs) != 0:
        new_df = pd.concat(list_new_dfs)
        new_df = select_only_positive_binary(new_df)
        new_df = select_only_single_positive_chemical_step(new_df)
        df = pd.concat([df, new_df])

    return df

def convert_reductive_amination(df, log):

    df_reaction = 'Reductive amination'
    reactions = ['Imine reduction', 'Imine formation',
                 'Iminium reduction', 'Iminium formation',
                 'Intramolecular iminium formation', 'Intramolecular imine formation']

    df = reaction_coversion(df, df_reaction, reactions, min_steps=2, max_steps=4, log=log)
    df = df.drop_duplicates()

    return df

def convert_reductive_amination_NH3(df, log):

    df_reaction = 'Reductive amination'
    reactions = ['Ketone amination', 'Aldehyde amination']

    df = reaction_coversion(df, df_reaction, reactions, min_steps=1, max_steps=1,
                            ignore_substrate_two=True, log=log)
    df = df.drop_duplicates()

    return df

def convert_alpha_beta_isomerisation(df, log):

    df_reaction = 'α- to β-amino isomerisation'
    reactions = ['α-amino amination', 'β-amino deamination']
    df = reaction_coversion(df, df_reaction, reactions, 2, ignore_substrate_two=True, log=log)
    df = df.drop_duplicates()

    return df



def do_conversions(df, log=True):
    df = convert_reductive_amination(df, log)
    df = convert_reductive_amination_NH3(df, log)
    df = convert_alpha_beta_isomerisation(df, log)
    return df