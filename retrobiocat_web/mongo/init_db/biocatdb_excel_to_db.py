from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Paper, Activity, Molecule
from pathlib import Path
import mongoengine as db
from retrobiocat_web.retro.enzyme_identification.load import make_fingerprints
from rdkit import Chem
import pandas as pd
import numpy as np
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.init_db.data_checks import check_is_float, check_is_nan, get_smile, get_mol

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
    if '.xlsx' in excel_path:
        df = pd.read_excel(excel_path)
    elif '.csv' in excel_path:
        df = pd.read_csv(excel_path)
    else:
        print('Could not load file')
        df = None

    return df

def df_to_db(spec_df):
    #added_by_dict = make_added_by_user_dict()

    print('Saving biocatdb_2 excel to mongodb..')
    for i, row in spec_df.iterrows():
        html_doi = str(row['html_doi'])
        doi = str(row['html_doi'])
        added_by_string = str(row['added_by'])

        list_html_to_remove = ['https://doi.org/', 'http://doi.org/', 'http://dx.doi.org/']
        for to_remove in list_html_to_remove:
            if to_remove in doi:
                doi = html_doi.replace(to_remove, '')

        if len(Paper.objects(doi=doi)) == 0:
            paper = Paper(short_citation=str(row['short_citation']),
                          html=html_doi,
                          doi=doi)
            paper = paper.save()
            print(f"{row['short_citation']} added")
        else:
            paper = Paper.objects(doi=doi)[0]

        if row['enzyme_type'] is not None and row['enzyme_type'] != '' and type(row['enzyme_type']) == str:
            if len(EnzymeType.objects(enzyme_type=row['enzyme_type'])) == 0:
                enz_type = EnzymeType(enzyme_type=row['enzyme_type'],
                                      description='')
                enz_type.save()

        if row['enzyme_name'] is not None and row['enzyme_name'] != '' and type(row['enzyme_name']) == str:
            if len(Sequence.objects(enzyme_name=row['enzyme_name'])) == 0:
                seq = Sequence(enzyme_name=check_is_nan(row['enzyme_name']),
                               enzyme_type=check_is_nan(row['enzyme_type']),
                               papers=[paper])
                seq.save()
            else:
                seq = Sequence.objects(enzyme_name=row['enzyme_name'])[0]
                if paper not in seq.papers:
                    seq.papers.append(paper)
                    seq = seq.save()

        if row['binary'] == 1:
            binary = True
        else:
            binary = False

        if row['auto_generated'] == 1:
            auto_gen = True
        else:
            auto_gen = False

        activity = Activity(enzyme_type=check_is_nan(row['enzyme_type']),
                            enzyme_name=check_is_nan(row['enzyme_name']),
                            reaction=check_is_nan(row['reaction']),
                            short_citation=check_is_nan(row['short_citation']),
                            html_doi=check_is_nan(row['html_doi']),
                            added_by_string=added_by_string,
                            paper=paper,
                            cascade_num=check_is_nan(row['cascade_num']),
                            substrate_1_smiles=get_smile(row['substrate_1_smiles']),
                            substrate_2_smiles=get_smile(row['substrate_2_smiles']),
                            product_1_smiles=get_smile(row['product_1_smiles']),
                            temperature=check_is_nan(row['temperature']),
                            ph=check_is_nan(row['ph']),
                            solvent=check_is_nan(row['solvent']),
                            other_conditions=check_is_nan(row['other_conditions']),
                            notes=check_is_nan(row['notes']),
                            reaction_vol=check_is_nan(row['reaction_vol']),
                            formulation=check_is_nan(row['formulation']),
                            biocat_conc=check_is_nan(row['biocat_conc']),
                            kcat=check_is_float(row['kcat']),
                            km=check_is_float(row['km']),
                            mw=check_is_float(row['mw']),
                            substrate_1_conc=check_is_nan(row['substrate_1_conc']),
                            substrate_2_conc=check_is_nan(row['substrate_2_conc']),
                            specific_activity=check_is_float(row['specific_activity']),
                            conversion=check_is_float(row['conversion']),
                            conversion_time=check_is_float(row['conversion_time']),
                            categorical=check_is_nan(row['categorical']),
                            binary=binary,
                            selectivity=check_is_nan(row['selectivity']),
                            auto_generated=auto_gen)

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

