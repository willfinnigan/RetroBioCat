from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, Paper, Molecule
from retrobiocat_web.mongo.models.reaction_models import Reaction

def clear_databases():
    EnzymeType.drop_collection()
    Sequence.drop_collection()
    Paper.drop_collection()
    Reaction.drop_collection()
    Activity.drop_collection()
    Molecule.drop_collection()

if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    clear_databases()

