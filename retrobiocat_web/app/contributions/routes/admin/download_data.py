from retrobiocat_web.app.contributions import bp
from flask import render_template, flash, redirect, url_for, make_response
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import Sequence
from retrobiocat_web.retro.enzyme_identification import query_mongodb
import yaml
import collections
import json
import pandas as pd


def setup_yaml():
  represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
  yaml.add_representer(collections.OrderedDict, represent_dict_order)

setup_yaml()


def reaction_to_yaml_dict(yaml_dict, reaction):
    name = reaction.name
    yaml_dict[name] = {}

    yaml_dict[name]['smarts'] = list(reaction.smarts)
    yaml_dict[name]['enzymes'] = dict(reaction.cofactors)
    yaml_dict[name]['positive_tests'] = list(reaction.positive_tests)
    yaml_dict[name]['negative_tests'] = list(reaction.negative_tests)
    yaml_dict[name]['type'] = reaction.type
    if reaction.experimental == True:
        yaml_dict[name]['experimental'] = True

    yaml_json = json.dumps(yaml_dict)
    yaml_dict = json.loads(yaml_json)

    return yaml_dict

def yaml_dict_to_yaml(yaml_dict):
    yaml_json = json.dumps(yaml_dict)
    yaml_dict = json.loads(yaml_json)
    return yaml.dump(yaml_dict)

@bp.route('/download_rxn_yaml', methods=['GET'])
@roles_required('admin')
def download_rxn_yaml():
    all_reactions = Reaction.objects()
    yaml_dict = {}
    for reaction in all_reactions:
        print(reaction)
        yaml_dict = reaction_to_yaml_dict(yaml_dict, reaction)

    print(yaml_dict)
    yaml_json = json.dumps(yaml_dict)
    yaml_dict = json.loads(yaml_json)
    resp = make_response(yaml.dump(yaml_dict))
    resp.headers["Content-Disposition"] = "attachment; filename=rxns_yaml.yaml"
    return resp

@bp.route('/download_biocatdb', methods=['GET'])
@roles_required('admin')
def download_biocatdb():
    spec_df = query_mongodb.query_specificity_data(['All'], ['All'])
    resp = make_response(spec_df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=biocatdb.csv"
    return resp

@bp.route('/download_sequences', methods=['GET'])
@roles_required('admin')
def download_sequences():
    sequences = Sequence.objects().as_pymongo()
    df = pd.DataFrame(list(sequences))
    resp = make_response(df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=sequences.csv"
    return resp

