from rdkit import Chem
from retrobiocat_web.app.biocatdb.functions import images
import numpy as np


def process_activity_data(activity_data):
    for i, record in enumerate(activity_data):
        activity_data[i]['paper'] = str(activity_data[i]['paper'])
        if 'id' in activity_data[i]:
            activity_data[i]['_id'] = str(activity_data[i]['id'])
        else:
            activity_data[i]['_id'] = str(activity_data[i]['_id'])

        if 'added_by' in activity_data[i]:
            activity_data[i]['added_by'] = str(activity_data[i]['added_by'])

        for key in activity_data[i]:
            if activity_data[i][key] == True:
                activity_data[i][key] = "True"
            if activity_data[i][key] == False:
                activity_data[i][key] = "False"
            if activity_data[i][key] == np.nan:
                activity_data[i][key] = ""
            if type(activity_data[i][key]) == float:
                activity_data[i][key] = round(activity_data[i][key], 2)

    return activity_data

def smiles_to_svg(activity_data):
    for i, record in enumerate(activity_data):
        if activity_data[i]["substrate_1_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["substrate_1_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["substrate_1_smiles"] = url

        if activity_data[i]["substrate_2_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["substrate_2_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["substrate_2_smiles"] = url

        if activity_data[i]["product_1_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["product_1_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["product_1_smiles"] = url

    return activity_data

