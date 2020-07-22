import pandas as pd
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from pathlib import Path
import numpy as np
from rdkit.Chem import rdFingerprintGenerator

COLUMNS = ['Reaction',
           'Enzyme type',
           'Enzyme name',
           'Data source',
           'Data source doi',
           'Cascade num',
           'Substrate 1 SMILES',
           'Substrate 2 SMILES',
           'Product 1 SMILES',
           'Temperature',
           'pH',
           'Solvent',
           'Other conditions',
           'Notes',
           'Reaction volume (ml)',
           'Biocatalyst Formulation',
           'Biocatalyst Concentration (mg/ml)',
           'kcat (min-1)',
           'KM (mM)',
           'Enz MW (Da)',
           'Substrate 1 conc (mM)',
           'Substrate 2 conc (mM)',
           'Specific activity (U/mg)',
           'Conversion (%)',
           'Conversion time (hrs)',
           'Categorical',
           'Binary',
           'Added by',
           'Selectivity',
           'Test Exceptions',
           'Auto Generated']

CATEGORIES = ['High', 'Medium', 'Low', 'None']
CAT_TYPE = pd.CategoricalDtype(categories=CATEGORIES, ordered=True)

CATS_CONVERSION = {'High': 60,
                   'Medium': 10,
                   'Low': 0.01,
                   'None': 0}

CATS_SA = {'High': 1,
           'Medium': 0.25,
           'Low': 0.05,
           'None': 0}

DEFAULT_CONC = 10

def set_columns(df):
    df = df.reindex(columns=COLUMNS)
    return df

def set_categorical_data_type(df, cat_type=CAT_TYPE):
    df['Categorical'] = df['Categorical'].astype(cat_type)
    return df

def load_biocatdb(excel_path):
    if '.xlsx' in excel_path:
        df = pd.read_excel(excel_path)
    elif '.csv' in excel_path:
        df = pd.read_csv(excel_path)
    else:
        print('Could not load file')
        return None

    df = set_columns(df)
    df = set_categorical_data_type(df)

    return df


def get_enzymes_with_data(df):
    enzymes = []
    for index, row in df.iterrows():
        if row['Enzyme type'] not in enzymes:
            enzymes.append(row['Enzyme type'])

    return enzymes

def get_reactions_with_data(df):
    reactions = []
    for index, row in df.iterrows():
        if row['Reaction'] not in reactions:
            reactions.append(row['Reaction'])

    return reactions

if __name__ == '__main__':
    pass