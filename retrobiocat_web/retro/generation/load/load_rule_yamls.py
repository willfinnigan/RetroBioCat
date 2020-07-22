import yaml
from pathlib import Path
from retrobiocat_web.retro.rdchiral.main import rdchiralReaction

def load_yamls(yaml_path):
    with open(yaml_path) as file:
        yaml_dict = yaml.load(file, Loader=yaml.FullLoader)

    return yaml_dict

def load_rxns(yaml_dict):
    rxns = {}
    for name in yaml_dict:
        if name not in rxns:
            rxns[name] = []
        for rxn_string in yaml_dict[name]['smarts']:
            try:
                rxns[name].append(rdchiralReaction(rxn_string))
            except:
                print('Error loading reaction - ' + str(name))

    return rxns

def load_rxn_strings(yaml_dict):
    rxns = {}
    for name in yaml_dict:
        if name not in rxns:
            rxns[name] = []
        for rxn_string in yaml_dict[name]['smarts']:
            try:
                rxns[name].append(rxn_string)
            except:
                print('Error loading reaction - ' + str(name))
    return rxns

def load_rules_by_type(yaml_dict):
    rules_by_type = {}
    for name in yaml_dict:
        rule_type = yaml_dict[name]['type']
        if rule_type not in rules_by_type:
            rules_by_type[rule_type] = []
        rules_by_type[rule_type].append(name)
    return rules_by_type

def load_tests(yaml_dict):
    rxn_postitive_tests = {}
    rxn_negative_tests = {}

    for name in yaml_dict:
        rxn_postitive_tests[name] = yaml_dict[name]['positive_tests']
        rxn_negative_tests[name] = yaml_dict[name]['negative_tests']

    return rxn_postitive_tests, rxn_negative_tests

def load_reactions_and_enzymes(yaml_dict):
    reactions = set()
    enzymes = set()
    reaction_enzyme_map = dict()

    for name in yaml_dict:
        reactions.add(name)
        reaction_enzyme_map[name] = set()
        for enz in yaml_dict[name]['enzymes']:
            enzymes.add(enz)
            reaction_enzyme_map[name].add(enz)
        reaction_enzyme_map[name] = list(reaction_enzyme_map[name])

    return list(reactions), list(enzymes), reaction_enzyme_map

def load_cofactors(yaml_dict):
    reactionEnzymeCofactorDict = {}

    for name in yaml_dict:
        reactionEnzymeCofactorDict[name] = {}

        for enz in yaml_dict[name]['enzymes']:
            cofactor_minus = yaml_dict[name]['enzymes'][enz]['cofactors_minus']
            cofactor_plus = yaml_dict[name]['enzymes'][enz]['cofactors_plus']
            reactionEnzymeCofactorDict[name][enz] = [cofactor_plus, cofactor_minus]

    return reactionEnzymeCofactorDict

if __name__ == '__main__':
    yaml_path = str(Path(__file__).parents[3]) + '/data/rxn_yaml/reaction_rules.yaml'
    yaml_dict = load_yamls(yaml_path)
    rxns = load_rxns(yaml_dict)
    reactions, enzymes, reaction_enzyme_map = load_reactions_and_enzymes(yaml_dict)
    reactionEnzymeCofactorDict = load_cofactors(yaml_dict)

    print(yaml_dict)