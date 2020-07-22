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
from retrobiocat_web.retro.enzyme_identification import query_mongodb


def fp_molecules_to_db(fp_df):
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

def get_spec_df_from_mongo():
    print('Load all of specdf from mongo..')
    spec_df = query_mongodb.query_specificity_data(['All'], ['All'])
    return spec_df

def make_fp_db():
    Molecule.drop_collection()
    spec_df = get_spec_df_from_mongo()
    fp_df = make_fingerprints.make_fingerprint_df_for_mongo(spec_df)

    fp_molecules_to_db(fp_df)


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    make_fp_db()

