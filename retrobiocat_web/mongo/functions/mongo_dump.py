from pathlib import Path
import subprocess as sp
from flask import current_app
from retrobiocat_web.app.app import db

def execute_mongo_dump():
    print("Executing mongo dump..")
    output_path = str(Path(__file__).parents[1]) + '/mongo_dump/mongo_dump.gz'
    if current_app.config['MONGODB_HOST'] == 'localhost':
        command = f'''mongodump --uri="mongodb://localhost:{current_app.config['MONGODB_PORT']}" --gzip --archive={output_path} --verbose'''
    else:
        command = f'''mongodump --uri="{current_app.config['MONGODB_HOST']}:{current_app.config['MONGODB_PORT']}" --gzip --archive={output_path} --verbose'''
    print(f"CMD = {command}")
    sp.run(command, shell=True)

def execute_mongo_restore(filename):
    db.connection.drop_database(current_app.config['MONGODB_DB'])
    print("Executing mongo restore..")

    if current_app.config['MONGODB_HOST'] == 'localhost':
        command = f'''mongorestore --uri="mongodb://localhost:{current_app.config['MONGODB_PORT']}" --db={current_app.config['MONGODB_DB']} --gzip --archive={filename} --verbose'''
    else:
        command = f'''mongorestore --uri="{current_app.config['MONGODB_HOST']}:{current_app.config['MONGODB_PORT']}" --db={current_app.config['MONGODB_DB']} --gzip --archive={filename} --verbose'''

    print(f"CMD = {command}")
    sp.run(command, shell=True)

