from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, Activity
from mongoengine.queryset.visitor import Q
import pandas as pd

COLUMNS = ['Reaction',
           'Enzyme type',
           'Enzyme name',
           'Data source',
           'Data source doi',
           'Cascade num',
           'Substrate 1 SMILES',
           'Substrate 2 SMILES',
           'Product 1 SMILES',
           'Temperature',
           'pH',
           'Solvent',
           'Other conditions',
           'Notes',
           'Reaction volume (ml)',
           'Biocatalyst Formulation',
           'Biocatalyst Concentration (mg/ml)',
           'kcat (min-1)',
           'KM (mM)',
           'Enz MW (Da)',
           'Substrate 1 conc (mM)',
           'Substrate 2 conc (mM)',
           'Specific activity (U/mg)',
           'Conversion (%)',
           'Conversion time (hrs)',
           'Categorical',
           'Binary',
           'Added by',
           'Selectivity',
           'Auto Generated',
           'paper',
           '_id']

def rename_mongo_columns(df):
    col_rename = {'reaction': 'Reaction',
                  'enzyme_type': 'Enzyme type',
                  'enzyme_name': 'Enzyme name',
                  'short_citation': 'Data source',
                  'html_doi': 'Data source doi',
                  'cascade_num': 'Cascade num',
                  'substrate_1_smiles': 'Substrate 1 SMILES',
                  'substrate_2_smiles': 'Substrate 2 SMILES',
                  'product_1_smiles': 'Product 1 SMILES',
                  'temperature': 'Temperature',
                  'ph': 'pH',
                  'solvent': 'Solvent',
                  'other_conditions': 'Other conditions',
                  'notes': 'Notes',
                  'reaction_vol': 'Reaction volume (ml)',
                  'formulation': 'Biocatalyst Formulation',
                  'biocat_conc': 'Biocatalyst Concentration (mg/ml)',
                  'kcat': 'kcat (min-1)',
                  'km': 'KM (mM)',
                  'mw': 'Enz MW (Da)',
                  'substrate_1_conc': 'Substrate 1 conc (mM)',
                  'substrate_2_conc': 'Substrate 2 conc (mM)',
                  'specific_activity': 'Specific activity (U/mg)',
                  'conversion': 'Conversion (%)',
                  'conversion_time': 'Conversion time (hrs)',
                  'categorical': 'Categorical',
                  'binary': 'Binary',
                  'added_by': 'Added by',
                  'selectivity': 'Selectivity',
                  'auto_generated': 'Auto Generated'}
    df.rename(columns=col_rename, inplace=True)
    df = df.reindex(columns=COLUMNS)
    return df

def query_specificity_data(listReactions, listEnzymes):
    if "All" in listEnzymes or len(listEnzymes) == 0 or listEnzymes == ['']:
        enzQ = Q()
    else:
        enzQ = Q(enzyme_type__in=listEnzymes)

    if "All" in listReactions or len(listReactions) == 0 or listReactions == ['']:
        reacQ = Q()
    else:
        reacQ = Q(reaction__in=listReactions)

    result = Activity.objects(enzQ & reacQ).as_pymongo()
    spec_df = pd.DataFrame(list(result))
    spec_df = rename_mongo_columns(spec_df)

    return spec_df

def get_reactions_in_db():
    reactions = Activity.objects().distinct('reaction')
    return reactions

def get_enzymes_in_db():
    enzymes = Activity.objects().distinct('enzyme_type')
    return enzymes

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    import time
    make_default_connection()

    t0 = time.time()
    df = query_specificity_data([], ['CAR'])
    t1 = time.time()
    print(f"Code took {round(t1-t0,3)} seconds to run")