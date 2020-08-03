from rdkit import Chem
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, Paper


def check_column(data_dict, col):
    issues = []
    if col in data_dict:
        if data_dict[col] == "" or data_dict[col] is None:
            return [f"Row {data_dict['n']}: Must specificy {col}"]
    else:
        return [f"Row {data_dict['n']}: Must specificy {col}"]
    return []

def check_smiles_is_accepted_by_rdkit(data_dict, col):
    smi = data_dict[col]
    try:
        mol = Chem.MolFromSmiles(smi)
        smi = Chem.MolToSmiles(mol)
    except:
        return [f"Row {data_dict['n']}: SMILE {smi} not accepted by rdkit"]
    return []

def check_all_required_are_numbers(data_dict):
    cols = ['conversion', 'specific_activity', 'km', 'kcat', 'mw']
    for col in cols:
        if col in data_dict:
            if not isinstance(data_dict[col], (int, float, complex)) or isinstance(data_dict[col], bool):
                if data_dict[col] != '':
                    return [f"Row {data_dict['n']}: {col} data is type string"]
    return []


def initial_check_data(data_dict):
    issues = []

    required_cols = ['reaction', 'enzyme_name']
    for col in required_cols:
        issues += check_column(data_dict, col)

    smi_cols = ["substrate_1_smiles", "substrate_2_smiles", "product_1_smiles"]
    for col in smi_cols:
        if col in data_dict:
            issues += check_smiles_is_accepted_by_rdkit(data_dict, col)

    issues += check_all_required_are_numbers(data_dict)

    return issues


def check_all_have_binary(data_dict):
    if 'binary' not in data_dict:
        return [f"Row {data_dict['n']}: No activity data has been entered - no binary"]
    if data_dict['binary'] != 1 and data_dict['binary'] != 0:
        return [f"Row {data_dict['n']}: No activity data has been entered - no binary"]
    return []

def check_seqs_are_defined(data_dict, paper_seqs):
    issues = []

    if data_dict['enzyme_name'] not in paper_seqs:
        issues.append(f"Enzymes must be defined in sequence tab first - {data_dict['enzyme_name']}")

    return issues






