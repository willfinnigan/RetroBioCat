import sqlite3
import pandas as pd


def make_db(name):
    conn = sqlite3.connect(name)
    c = conn.cursor()
    c.execute('''CREATE TABLE BUYABLE
                 ([id] INTEGER PRIMARY KEY,[SMILES] text)''')
    conn.commit()
    return name

def df_to_sql(db_path, df, table='BUYABLE'):
    conn = sqlite3.connect(db_path)
    df.to_sql(table, conn, if_exists='replace')

def create_index(db_path, name, table='BUYABLE'):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute(f"CREATE INDEX index_{name} ON {table}({name})")
    conn.commit()


if __name__ == '__main__':
    zinc_df = pd.read_csv('zinc_in_stock.csv')
    emol_df = pd.read_csv('salt_removed_emolecules.csv')
    mol_df = pd.read_csv('molport_building_blocks.csv')
    #ss_df = pd.read_csv('ss.csv')

    df = pd.concat([zinc_df, emol_df, mol_df])
    df_smi = df[['SMILES']]
    df_smi = df_smi.drop_duplicates()

    df.to_csv('building_blocks.csv')

    print(df_smi.info())
    print(df_smi.head())

    #db_name = 'building_blocks_small.db'
    #make_db(db_name)
    #df_to_sql(db_name, df_smi)
    #create_index(db_name, 'SMILES')





