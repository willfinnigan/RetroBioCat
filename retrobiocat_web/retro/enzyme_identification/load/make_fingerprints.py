from rdkit import Chem
import pandas as pd
from rdkit.Chem import rdFingerprintGenerator
from retrobiocat_web.mongo.models.biocatdb_models import Molecule
from rdkit import DataStructs


fingerprint_options = {'mfp_default': ['morgan', {}],
                       'rdfp_default': ['rdkit', {}]}

default_fp_mode = 'rdfp_default'

def make_fp_generator(fp_type, settings):
    if fp_type == 'morgan':
        arguments = {'includeChirality': True}
        for arg in settings:
            arguments[arg] = settings[arg]
        fp_gen = rdFingerprintGenerator.GetMorganGenerator(**arguments)

    elif fp_type == 'atom_pair':
        arguments = {'includeChirality': True}
        for arg in settings:
            arguments[arg] = settings[arg]

        fp_gen = rdFingerprintGenerator.GetAtomPairGenerator(**arguments)

    elif fp_type == 'rdkit':
        arguments = {}
        for arg in settings:
            arguments[arg] = settings[arg]

        fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator(**arguments)

    elif fp_type == 'toplogical':
        arguments = {'includeChirality': True}
        for arg in settings:
            arguments[arg] = settings[arg]
        fp_gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(**arguments)
    else:
        fp_gen = False

    return fp_gen

def unique_smiles_list(df):

    big_list = pd.concat([df['Substrate 1 SMILES'], df['Substrate 2 SMILES'], df['Product 1 SMILES']])
    unique = list(big_list.unique())

    if '' in unique:
        unique.remove('')

    return unique

def make_fingerprint_df(df, fp_type, fp_settings):
    unique_smiles = unique_smiles_list(df)
    list_smi = []
    list_mol = []
    for i, smi in enumerate(unique_smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            list_smi.append(smi)
            list_mol.append(mol)
        except:
            pass

    fp_gen = make_fp_generator(fp_type, fp_settings)

    #fingerprint_df = pd.DataFrame({'smiles' : list_smi, 'mol' : list_mol, 'fp' : fp})
    fingerprint_df = pd.DataFrame({'smiles': list_smi, 'mol': list_mol})
    fingerprint_df['fp'] = fingerprint_df['mol'].apply(fp_gen.GetFingerprint)
    fingerprint_df.set_index('smiles')
    fingerprint_df.drop(columns=['mol'], inplace=True)
    #fingerprint_df.dropna(subset=['fp'], inplace=True)

    return fingerprint_df

def make_fingerprint_df_for_mongo(df):
    unique_smiles = unique_smiles_list(df)
    list_smi = []
    list_mol = []
    for i, smi in enumerate(unique_smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            list_smi.append(smi)
            list_mol.append(mol)
        except:
            pass

    fingerprint_df = pd.DataFrame({'smiles': list_smi, 'mol': list_mol})
    fingerprint_df.set_index('smiles')

    for option_name, options in fingerprint_options.items():
        fp_gen = make_fp_generator(options[0], options[1])
        fingerprint_df[option_name] = fingerprint_df['mol'].apply(fp_gen.GetFingerprint)

    fingerprint_df.drop(columns=['mol'], inplace=True)

    return fingerprint_df

def load_fp_df_from_mongo(mode):
    def convert_to_fp(bitstring):
        return DataStructs.CreateFromBitString(bitstring)

    query_result = Molecule.objects.as_pymongo().only('smiles', mode)
    fp_df = pd.DataFrame(list(query_result))

    fp_df[mode] = fp_df[mode].apply(convert_to_fp)
    fp_df.rename(columns={mode:'fp'}, inplace=True)

    return fp_df

if __name__ == "__main__":
    import time
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    t0 = time.time()
    fp_df2 = load_fp_df_from_mongo(default_fp_mode)
    t1 = time.time()
    print(fp_df2.head())
    print(f"Time to load fingerprints from mongo = {round(t1 - t0, 4)} seconds")