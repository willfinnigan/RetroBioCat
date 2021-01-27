from pathlib import Path
import requests
import os
import json
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.mongo.models.biocatdb_models import ActivityMol
from rdkit import Chem
from retrobiocat_web.app.retrobiocat.functions.get_images import moltosvg

DATA_FOLDER = str(Path(__file__).parents[0]) + '/uploaded_molecule_images'
ALLOWED_FILE_TYPES = ['.png', '.jpeg']

def save_uploaded_image(file, filename):
    save_path = f'{DATA_FOLDER}/{filename}'

    for filetype in ALLOWED_FILE_TYPES:
        if filename.endswith(filetype) == True:
            file.save(save_path)
            return True
    return False

def process_image(filename, osra_api_url, remove=True, log=False):
    save_path = f'{DATA_FOLDER}/{filename}'
    for filetype in ALLOWED_FILE_TYPES:
        if filename.endswith(filetype) == True:
            url = f'{osra_api_url}/process_image'
            files = {'image': open(save_path, 'rb')}
            r = requests.post(url, files=files)
            if remove == True:
                os.remove(save_path)

            text = r.text
            list_smi = json.loads(text)

            # if can load via rdkit, return the rdkit smile.
            list_rdkit_smi = []
            for smi in list_smi:
                rd_smi = rdkit_smile(smi)
                if rd_smi != None:
                    list_rdkit_smi.append(rd_smi)
                else:
                    list_rdkit_smi.append(smi)

            if log == True:
                print(f"-- Processed file: {filename} --")
                for smi in list_rdkit_smi:
                    print(smi)

            return list_rdkit_smi
    return False

def save_mol(smi, paper, name=''):
    try:
        mol = Chem.MolFromSmiles(smi)
        svg = moltosvg(mol)
    except:
        svg = ""

    new_mol = ActivityMol(smi=smi,
                          name=name,
                          paper=paper,
                          svg=svg)
    new_mol.save()


if __name__ == '__main__':
    process_image('test1.png', 'http://0.0.0.0:8080', remove=False, log=True)
