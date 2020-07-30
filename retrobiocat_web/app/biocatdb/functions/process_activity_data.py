from rdkit import Chem
from retrobiocat_web.app.biocatdb.functions import images


def process_activity_data(activity_data):
    for i, record in enumerate(activity_data):
        activity_data[i]['paper'] = str(activity_data[i]['paper'])
        activity_data[i]['activity_id'] = str(activity_data[i]['activity_id'])

    return activity_data

def smiles_to_svg(activity_data):
    for i, record in enumerate(activity_data):
        if activity_data[i]["Substrate 1 SMILES"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["Substrate 1 SMILES"])
            url = images.moltosvg_url(mol)
            activity_data[i]["Substrate 1 SMILES"] = url

        if activity_data[i]["Substrate 2 SMILES"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["Substrate 2 SMILES"])
            url = images.moltosvg_url(mol)
            activity_data[i]["Substrate 2 SMILES"] = url

        if activity_data[i]["Product 1 SMILES"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["Product 1 SMILES"])
            url = images.moltosvg_url(mol)
            activity_data[i]["Product 1 SMILES"] = url

    return activity_data

