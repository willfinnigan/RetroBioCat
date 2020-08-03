import os
import yaml
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Paper, Activity, Molecule
from pathlib import Path
import mongoengine as db
from retrobiocat_web.retro.enzyme_identification.load import make_fingerprints, load_spec_df_excel
from rdkit import Chem
import pandas as pd
import numpy as np
from retrobiocat_web.mongo.models.user_models import User

EXCEL_PATH = str(Path(__file__).parents[2]) + '/data/substrate_specificity/biocatdb_2.xlsx'

def make_added_by_user_dict():
    added_by_dict = {}
    users = User.objects()
    for user in users:
        added_by_dict[f"{user.first_name} {user.last_name}"] = user
    return added_by_dict

def get_user(added_by, added_by_dict):
    for key in added_by_dict:
        if added_by in key:
            print(added_by_dict[key])
            return added_by_dict[key]

def check_is_float(to_test):
    if type(to_test) != float:
        try:
            to_test = float(to_test)
        except:
            return None
    if np.isnan(to_test) == True:
        return None

    return to_test

def check_is_nan(to_test):
    if str(to_test) == 'NaN' or str(to_test) == '' or str(to_test) == ' ':
        return None
    if type(to_test) != str:
        if isinstance(to_test, (int, float, complex)):
            if np.isnan(to_test) == True:
                return None
            else:
                to_test = str(to_test)
    return to_test

def get_mol(smile):
    try:
        mol = Chem.MolFromSmiles(smile)
    except:
        mol = None
    return mol

def get_smile(smile):
    smile = check_is_nan(smile)
    if smile == None:
        return ''
    return smile

def fp_molecules_to_db(fp_df):
    Molecule.drop_collection()

    print("Save fingerprints into molecule objects..")
    for i, row in fp_df.iterrows():
        if len(Molecule.objects(smiles=row['smiles'])) == 0:
            mol = Molecule()
        else:
            mol = Molecule.objects(smiles=row['smiles'])[0]

        mol.smiles = row['smiles']
        mol.mol = Chem.MolFromSmiles(row['smiles']).ToBinary()

        for option_name in make_fingerprints.fingerprint_options:
            mol[option_name] = row[option_name].ToBitString()

        mol.save()
    print("..done")

def load_df(excel_path):
    print("Load Spec DF..")
    spec_df = load_spec_df_excel.load_biocatdb(excel_path)
    return spec_df

def df_to_db(spec_df):
    added_by_dict = make_added_by_user_dict()

    print('Saving biocatdb_2 excel to mongodb..')
    for i, row in spec_df.iterrows():
        html_doi = row['Data source doi']
        doi = row['Data source doi']
        added_by_string = str(row['Added by'])

        list_html_to_remove = ['https://doi.org/', 'http://doi.org/', 'http://dx.doi.org/']
        for to_remove in list_html_to_remove:
            if to_remove in doi:
                doi = html_doi.replace(to_remove, '')

        if len(Paper.objects(doi=doi)) == 0:
            paper = Paper(short_citation=row['Data source'],
                          html=html_doi,
                          doi=doi)
            paper = paper.save()
            print(f"{row['Data source']} added")
        else:
            paper = Paper.objects(doi=doi)[0]

        if row['Enzyme type'] is not None and row['Enzyme type'] != '' and type(row['Enzyme type']) == str:
            if len(EnzymeType.objects(enzyme_type=row['Enzyme type'])) == 0:
                enz_type = EnzymeType(enzyme_type=row['Enzyme type'],
                                      description='')
                enz_type.save()

        if row['Enzyme name'] is not None and row['Enzyme name'] != '' and type(row['Enzyme name']) == str:
            if len(Sequence.objects(enzyme_name=row['Enzyme name'])) == 0:
                seq = Sequence(enzyme_name=check_is_nan(row['Enzyme name']),
                               enzyme_type=check_is_nan(row['Enzyme type']),
                               papers=[paper])
                seq.save()
            else:
                seq = Sequence.objects(enzyme_name=row['Enzyme name'])[0]
                if paper not in seq.papers:
                    seq.papers.append(paper)
                    seq = seq.save()

        if row['Binary'] == 1:
            binary = True
        else:
            binary = False

        if row['Auto Generated'] == 1:
            auto_gen = True
        else:
            auto_gen = False

        activity = Activity(enzyme_type=check_is_nan(row['Enzyme type']),
                            enzyme_name=check_is_nan(row['Enzyme name']),
                            reaction=check_is_nan(row['Reaction']),
                            short_citation=check_is_nan(row['Data source']),
                            html_doi=check_is_nan(row['Data source doi']),
                            added_by_string=added_by_string,
                            paper=paper,
                            cascade_num=check_is_nan(row['Cascade num']),
                            substrate_1_smiles=get_smile(row['Substrate 1 SMILES']),
                            substrate_2_smiles=get_smile(row['Substrate 2 SMILES']),
                            product_1_smiles=get_smile(row['Product 1 SMILES']),
                            temperature=check_is_nan(row['Temperature']),
                            ph=check_is_nan(row['pH']),
                            solvent=check_is_nan(row['Solvent']),
                            other_conditions=check_is_nan(row['Other conditions']),
                            notes=check_is_nan(row['Notes']),
                            reaction_vol=check_is_nan(row['Reaction volume (ml)']),
                            formulation=check_is_nan(row['Biocatalyst Formulation']),
                            biocat_conc=check_is_nan(row['Biocatalyst Concentration (mg/ml)']),
                            kcat=check_is_float(row['kcat (min-1)']),
                            km=check_is_float(row['KM (mM)']),
                            mw=check_is_float(row['Enz MW (Da)']),
                            substrate_1_conc=check_is_nan(row['Substrate 1 conc (mM)']),
                            substrate_2_conc=check_is_nan(row['Substrate 2 conc (mM)']),
                            specific_activity=check_is_float(row['Specific activity (U/mg)']),
                            conversion=check_is_float(row['Conversion (%)']),
                            conversion_time=check_is_float(row['Conversion time (hrs)']),
                            categorical=check_is_nan(row['Categorical']),
                            binary=binary,
                            selectivity=check_is_nan(row['Selectivity']),
                            auto_generated=auto_gen,
                            )

        activity.save()
    print('..done')

def spec_df_to_db(excel_path=EXCEL_PATH):
    spec_df = load_df(excel_path)
    df_to_db(spec_df)

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    #spec_df_to_db()
    print(make_added_by_user_dict())

