import sqlite3
import pandas as pd


def make_db(path):
    conn = sqlite3.connect(path)
    c = conn.cursor()
    c.execute('''CREATE TABLE BUYABLE
                 ([id] INTEGER PRIMARY KEY,[SMILES] text)''')
    conn.commit()
    return path

def df_to_sql(db_path, df, table='BUYABLE'):
    conn = sqlite3.connect(db_path)
    df.to_sql(table, conn, if_exists='replace')

def create_index(db_path, name, table='BUYABLE'):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute(f"CREATE INDEX index_{name} ON {table}({name})")
    conn.commit()