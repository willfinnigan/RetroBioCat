from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, Activity
from mongoengine.queryset.visitor import Q
import pandas as pd

COLUMNS = ['reaction',
           'enzyme_type',
           'enzyme_name',
           'short_citation',
           'html_doi',
           'cascade_num',
           'substrate_1_smiles',
           'substrate_2_smiles',
           'product_1_smiles',
           'temperature',
           'ph',
           'solvent',
           'other_conditions',
           'notes',
           'reaction_vol',
           'formulation',
           'biocat_conc',
           'kcat',
           'km',
           'mw',
           'substrate_1_conc',
           'substrate_2_conc',
           'specific_activity',
           'conversion',
           'conversion_time',
           'categorical',
           'binary',
           'added_by',
           'selectivity',
           'auto_generated',
           'paper',
           '_id']

def organise_mongo_columns(df):
    #col_rename = {'_id': 'activity_id'}
    #df.rename(columns=col_rename, inplace=True)
    df = df.reindex(columns=COLUMNS)
    return df

def query_specificity_data(listReactions, listEnzymes, only_reviewed=False):
    if "All" in listEnzymes or len(listEnzymes) == 0 or listEnzymes == ['']:
        enzQ = Q()
    else:
        enzQ = Q(enzyme_type__in=listEnzymes)

    if "All" in listReactions or len(listReactions) == 0 or listReactions == ['']:
        reacQ = Q()
    else:
        reacQ = Q(reaction__in=listReactions)

    if only_reviewed == True:
        revQ = Q(reviewed=True)
    else:
        revQ = Q()

    result = Activity.objects(enzQ & reacQ & revQ).as_pymongo()
    spec_df = pd.DataFrame(list(result))
    spec_df = organise_mongo_columns(spec_df)

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
    print(df)