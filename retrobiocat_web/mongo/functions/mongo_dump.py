from pathlib import Path
import subprocess as sp
from flask import current_app

def execute_mongo_dump():
    print("Executing mongo dump..")
    output_path = str(Path(__file__).parents[1]) + '/mongo_dump'
    command = f"mongodump --host='{current_app.config['MONGODB_HOST']}:{current_app.config['MONGODB_PORT']}' -o {output_path} &"
    print(f"CMD = {command}")
    sp.run(command, shell=True)
