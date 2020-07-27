from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.retro.enzyme_identification import molecular_similarity
import pandas as pd
from rdkit.Chem import PandasTools
import os
from rdkit import RDConfig
from rdkit import Chem

from retrobiocat_web.app.biocatdb import bp
from rq import get_current_job



def filter_df_by_data_level(df, data_level):
    if data_level != 'All':
        if data_level == 'Categorical':
            df = df[df['Categorical'].notnull()]
        elif data_level == 'Quantitative':
            df1 = df[df['Specific activity (U/mg)'].notnull()]
            df2 = df[df['Conversion (%)'].notnull()]
            df = pd.concat([df1, df2])
        elif data_level == 'Specific Activity':
            df = df[df['Specific activity (U/mg)'].notnull()]
        elif data_level == 'Conversion':
            df = df[df['Conversion (%)'].notnull()]

    return df

def activity_columns_to_show(df):
    show_cat = False
    show_sa = False
    show_conv = False
    show_smiles_two = False
    show_kinetics = False

    for index, row in df.iterrows():
        if pd.notnull(row['Categorical']):
            show_cat = True
        if pd.notnull(row['Specific activity (U/mg)']):
            show_sa = True
        if pd.notnull(row['Conversion (%)']):
            show_conv = True
        if pd.notnull(row['Substrate 2 SMILES']):
            show_smiles_two = True
        if pd.notnull(row['kcat (min-1)']):
            show_kinetics = True

    return show_cat, show_sa, show_conv, show_smiles_two, show_kinetics

def only_required_activity_columns(df):

    show_cat, show_sa, show_conv, show_smiles_two, show_kinetics = activity_columns_to_show(df)

    if show_cat != True:
        df = df.drop(columns=['Categorical'])
    if show_sa != True:
        df = df.drop(columns=['Specific activity (U/mg)'])
    if show_conv != True:
        df = df.drop(columns=['Conversion (%)', 'Conversion time (hrs)'])
    if show_smiles_two != True:
        df = df.drop(columns=['Substrate 2 SMILES'])
        df = df.drop(columns=['Substrate 2 conc (mM)'])
        print('dropping smiles 2')
    if show_kinetics != True:
        df = df.drop(columns=['kcat (min-1)'])
        df = df.drop(columns=['KM (mM)'])
    return df

def remove_unwanted_columns(df):
    unwanted = ['Enz MW (Da)', 'Test Exceptions', 'Cascade num']
    df = df.drop(columns=unwanted)
    df = only_required_activity_columns(df)
    return df

def reorder_df_f(df):
    columns = ['Reaction',
               'Enzyme type',
               'Enzyme name',
               'Substrate 1 SMILES',
               'Substrate 2 SMILES',
               'Product 1 SMILES',
               'Binary',
               'Categorical',
               'Specific activity (U/mg)',
               'Conversion (%)',
               'Conversion time (hrs)',
               'Substrate 1 conc (mM)',
               'Substrate 2 conc (mM)',
               'Selectivity',
               'Data source',
               'Data source doi',
               'Temperature',
               'pH',
               'Solvent',
               'Other conditions',
               'Notes',
               'Reaction volume (ml)',
               'Biocatalyst Formulation',
               'Biocatalyst Concentration (mg/ml)',
               'Added by',
               'Auto Generated',
               'kcat (min-1)',
               'KM (mM)']

    columns = columns + (df.columns.drop(columns).tolist())
    df = df.reindex(columns=columns)

    return df

def get_headings(df):
    columns = list(df.columns)

    headings = ['Reaction',
               'Enzyme type',
               'Enzyme name',
               'Substrate 1',
               'Substrate 2',
               'Product 1',
               'Binary',
               'Categorical',
               'Specific activity (U/mg)',
               'Conversion (%)',
               'Conversion time (hrs)',
               'Substrate 1 conc (mM)',
               'Substrate 2 conc (mM)',
               'Selectivity',
               'Data source',
               'Temperature (celsius)',
               'pH',
               'Solvent',
               'Other conditions',
               'Notes',
               'Reaction volume (ml)',
               'Biocatalyst Formulation',
               'Biocatalyst Concentration (mg/ml)',
               'Added by',
               'Auto Generated',
               'kcat (min-1)',
               'KM (mM)',
               ]

    if 'Categorical' not in columns:
        headings.remove('Categorical')
    if 'Specific activity (U/mg)' not in columns:
        headings.remove('Specific activity (U/mg)')
    if 'Conversion (%)' not in columns:
        headings.remove('Conversion (%)')
        headings.remove('Conversion time (hrs)')
    if 'Substrate 2 SMILES' not in columns:
        headings.remove('Substrate 2')
        headings.remove('Substrate 2 conc (mM)')
    if 'kcat (min-1)' not in columns:
        headings.remove('kcat (min-1)')
        headings.remove('KM (mM)')

    return headings

def get_table(df):
    headings = get_headings(df)
    columns = list(df.columns)
    html_rows = []
    activity_columns = ['Categorical',
                        'Specific activity (U/mg)',
                        'Conversion (%)']

    for index, row in df.iterrows():
        single_row = []
        for i, col in enumerate(columns):
            content = {'content':'', 'style':''}
            if col == 'Data source':
                content['content'] = '<a href="' + str(row[columns[i+1]]) + '" target="_blank">' + str(row[col]) + '</a>'

            elif (col == 'Substrate 1 SMILES') or (col == 'Substrate 2 SMILES') or (col == 'Product 1 SMILES'):
                content['style'] = 'font-size: x-small; \n padding: 0px;'
                try:
                    img = str(Chem.MolFromSmiles(row[col]))
                    img = img[:-1] + 'width="100" height="100" >'
                except:
                    img = ''
                content['content'] = img
                #content['content'] += '\n <br>' #+ str(row[col])
                if (col == 'Product 1 SMILES') and ('Similarity' in columns):
                    content['content'] += '<br>' + 'Similarity score = ' + str(row['Similarity'])

            elif col == 'Other conditions' or col == 'Notes' or col == 'Added by':
                content['content'] = row[col]
                content['style'] = 'font-size: x-small;'

            elif col=='Binary':
                content['content']= row[col]
                if row[col] == 1:
                    content['style'] = "background-color:#C7FFCD;"
                else:
                    content['style'] = "background-color:#FFD4D4;"

            elif col in activity_columns:
                content['content'] = row[col]
                if pd.notnull(row[col]) == True:
                    if row['Categorical'] == 'High':
                        content['style'] = "background-color:#8DFF89;"
                    elif row['Categorical'] == 'Medium':
                        content['style'] = "background-color:#FDFF89;"
                    elif row['Categorical'] == 'Low':
                        content['style'] = "background-color:#FFC089;"
                    elif row['Categorical'] == 'None':
                        content['style'] = "background-color:#FF8989;"

            else:
                if pd.isnull(row[col]) != True:
                    content['content']= row[col]

                else:
                    content['content'] = ''

            if pd.isnull(content['content']) == True:
                content['content'] = ''

            if col not in ['Data source doi', 'Similarity']:
                single_row.append(content)

        html_rows.append(single_row)

    return headings, html_rows



def run_query(product, reactions, enzyme_names,
              data_level,
              num_choices,
              similarity_cutoff,
              max_hits,
              include_auto_data):

    scorer = molecular_similarity.SubstrateSpecificityScorer(print_log=False)

    print(reactions)
    print(enzyme_names)

    df = scorer.querySpecificityDf(product, reactions, enzyme_names,
                                   dataLevel=data_level,
                                   numEnzymes=num_choices,
                                   simCutoff=similarity_cutoff,
                                   numHits=max_hits,
                                   include_auto_generated=include_auto_data)

    if df is None:
        return ['No data'], ['No data']
    elif len(df.index) == 0:
        return ['No data'], ['No data']
    else:
        df = reorder_df_f(df)
        df = remove_unwanted_columns(df)
        return get_table(df)

def task_run_query(form_data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()
    print('Started')

    enzyme_names = list(form_data['enzymes'].split(", "))
    reactions = list(form_data['reactions'].split(", "))
    try:
        product = rdkit_smile(form_data['target_smiles'])
    except:
        return {'headings': ['SMILE NOT ACCEPTED'],
                'rows': ['SMILE NOT ACCEPTED']}
    similarity_cutoff = form_data['similarity']
    num_choices = form_data['num_choices']
    data_level = form_data['data_level']
    max_hits = form_data['max_hits']
    include_auto_data = bool(form_data['auto_data'])

    headings, rows = run_query(product, reactions, enzyme_names,
                               data_level, num_choices, similarity_cutoff,
                               max_hits, include_auto_data)

    return {'headings' : headings,
            'rows' : rows}