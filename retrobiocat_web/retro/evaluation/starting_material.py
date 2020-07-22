import os
import pandas as pd
from pathlib import Path
import sqlite3
import time


data_folder = str(Path(__file__).parents[2]) + '/retro/data/buyability'


class StartingMaterialEvaluator():

    if os.path.exists(data_folder + '/building_blocks.db'):
        sqlite_db_path = data_folder + '/building_blocks.db'
    else:
        sqlite_db_path = data_folder + '/building_blocks_small.db'


    def __init__(self, alternative_db_path=None, print_log=False):
        self.print_log = print_log

        if alternative_db_path != None:
            self.conn = sqlite3.connect(alternative_db_path)
        elif os.path.exists(self.sqlite_db_path):
            self.conn = sqlite3.connect(self.sqlite_db_path)
        else:
            self.conn = False

    def eval(self, smiles):
        cursor = self.conn.cursor()
        query = f"SELECT EXISTS(SELECT 1 FROM BUYABLE WHERE SMILES='{smiles}');"

        cursor.execute(query)
        result = cursor.fetchall()[0][0]
        self._log(f"Query for {smiles} = {result}")
        return result

    def _log(self, msg):
        if self.print_log == True:
            print(msg)


if __name__ == '__main__':

    ev = StartingMaterialEvaluator()

    t0 = time.time()
    print(ev.eval('O=CC(=O)Cc1ccccc1'))
    t1 = time.time()
    print(f"Time  = {round(t1 - t0, 4)}")

    t0 = time.time()
    print(ev.eval('O=CC(=O)Cc1ccccc1'))
    t1 = time.time()
    print(f"Time  = {round(t1 - t0, 4)}")

    t0 = time.time()
    print(ev.eval('CCCCC=O'))
    t1 = time.time()
    print(f"Time  = {round(t1 - t0, 4)}")