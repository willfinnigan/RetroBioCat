import os
import yaml
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from pathlib import Path
import mongoengine as db

YAML_PATH = str(Path(__file__).parents[2]) + '/data/rxn_yaml/reaction_rules.yaml'

def load_yaml_dict(yaml_path):
    with open(yaml_path) as file:
        yaml_dict = yaml.load(file, Loader=yaml.FullLoader)
    return yaml_dict

def load_into_mongo(yaml_dict):
    for rxn_name in yaml_dict:
        enzymes = []
        for enz in list(yaml_dict[rxn_name]['enzymes'].keys()):
            enzymes.append(enz)
            if len(EnzymeType.objects(enzyme_type=enz)) == 0:
                new_type = EnzymeType(enzyme_type=enz,
                                      description='')
                new_type.save()

        if len(Reaction.objects(name=rxn_name)) == 0:
            reaction = Reaction(name=rxn_name,
                                smarts=yaml_dict[rxn_name]['smarts'],
                                enzyme_types=enzymes,
                                cofactors=yaml_dict[rxn_name]['enzymes'],
                                positive_tests=yaml_dict[rxn_name]['positive_tests'],
                                negative_tests=yaml_dict[rxn_name]['negative_tests'],
                                type=yaml_dict[rxn_name]['type'])

            if 'experimental' in yaml_dict[rxn_name]:
                reaction.experimental = bool(yaml_dict[rxn_name]['experimental'])

            if 'two_step' in yaml_dict[rxn_name]:
                reaction.two_step = bool(yaml_dict[rxn_name]['two_step'])

            reaction.save()

def load_reaction_rules_into_mongo(yaml_file=YAML_PATH):
    yaml_dict = load_yaml_dict(yaml_file)
    load_into_mongo(yaml_dict)




if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    Reaction.drop_collection()
    load_reaction_rules_into_mongo()

    test_reaction = Reaction.objects.get(name='Primary alcohol oxidation')
    print(test_reaction.smarts)
