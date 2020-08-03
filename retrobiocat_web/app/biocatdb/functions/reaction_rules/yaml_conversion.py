import yaml
import collections
import json


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