import pandas as pd
import numpy as np
import sqlite3
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile

def convert_to_rdkit(smi):
    try:
        new_smi = rdkit_smile(smi)
        return new_smi
    except:
        return None

def load_data(path, cols, sep, smi_col):
    print(f'Load path: {path}')
    try:
        data = pd.read_csv(path, sep=sep)
        data.columns = cols
        df_smiles = data[[smi_col]]
        df_smiles.rename(columns={smi_col: "SMILES"})

        print('Data loaded - processing smiles..')
        df_smiles["SMILES"].apply(convert_to_rdkit)
        df_smiles["SMILES"].replace('', np.nan, inplace=True)
        df_smiles = df_smiles.dropna()
        df_smiles = df_smiles.drop_duplicates()
    except:
        df_smiles = pd.DataFrame()
    return df_smiles




if __name__ == '__main__':
    ss_cols = ['SMILES', 'code']
    ss_path = 'http://files.docking.org/catalogs/source/sialbb.src.txt'
    df = load_data(ss_path, ss_cols, ' ', "SMILES")
    df.to_csv('ss.csv')

    '''
    zinc_cols = ["SMILES", "zinc", 'data_type', 'supplier', 'code']
    zinc_path = 'http://files.docking.org/bb/current/In-Stock.txt'
    df = load_data(zinc_path, zinc_cols, '\t', "SMILES")
    df.to_csv('zinc_in_stock.csv')
    
    emol_cols = ['isosmiles', 'version_id', 'parent_id']
    emol_path = "emols_buyable.smi"
    df = load_data(emol_path, emol_cols, ' ', 'isosmiles')

    sigma_path = 'http://files.docking.org/catalogs/source/sialbb.src.txt'
    sigma_cols = ['SMILES', 'code']
    df = load_data(sigma_path, sigma_cols, ' ', 'SMILES')

    molport_path = 'http://files.docking.org/catalogs/source/molportbb.src.txt'
    molport_cols = ['SMILES', 'code']
    df = load_data(molport_path, molport_cols, ' ', 'SMILES')
    '''










