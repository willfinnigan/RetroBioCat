from retrobiocat_web.app.biocatdb import bp
from retrobiocat_web.retro.rdchiral.main import rdchiralReaction
from flask import render_template, jsonify, request
import yaml
from retrobiocat_web.app.retrobiocat.functions.get_images import rxntosvg
from flask_security import roles_required
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from retrobiocat_web.app.biocatdb.functions.reaction_rules import yaml_conversion
from retrobiocat_web.app.biocatdb.functions.reaction_rules.reaction_tests import ReactionTester
from retrobiocat_web.retro.generation.network_generation.network import Network

from retrobiocat_web.mongo.models.biocatdb_models import Activity

def get_reactions():
    reactions = list(Reaction.objects().distinct('name'))
    reactions.sort()
    return reactions

def get_enzymes():
    enzymes = list(EnzymeType.objects().distinct('enzyme_type'))
    enzymes.sort()
    return enzymes

def get_rxn_smarts_yaml(rxn_yaml_dict):
    list_smarts = rxn_yaml_dict['smarts']
    return yaml.dump(list_smarts)

def get_enzyme_cofactor_yaml(rxn_yaml_dict):
    enzymes = rxn_yaml_dict['enzymes']
    return yaml.dump(enzymes)

def get_reaction_yaml_dict(reaction_name):
    reaction = Reaction.objects(name=reaction_name)[0]
    yaml_dict = yaml_conversion.reaction_to_yaml_dict({}, reaction)
    return yaml_dict

def get_positive_tests(rxn_yaml_dict):
    pos = rxn_yaml_dict['positive_tests']
    return yaml.dump(pos)

def get_negative_tests(rxn_yaml_dict):
    neg = rxn_yaml_dict['negative_tests']
    return yaml.dump(neg)


@bp.route('/rule_editor')
@roles_required('rxn_rules_admin')
def rule_editor():
    reactions = get_reactions()
    enzymes = get_enzymes()
    return render_template('rxn_rule_editor/new_rxn_rules_page.html', reactions=reactions, enzymes=enzymes)


@bp.route("/_load_rule", methods=["POST"])
@roles_required('rxn_rules_admin')
def load_rule():
    reaction_name = request.form['selected_rule']

    if reaction_name == 'Empty template':
        rxn_name = ""
        rxn_smarts = ""
        rxn_enzyme_cofactor = ""
        positive_tests = ""
        negative_tests = ""
        reaction_type = ""
        experimental = True
        two_step = False
        requires_absence_of_water = False
    else:
        yaml_dict = get_reaction_yaml_dict(reaction_name)
        rxn_name = reaction_name
        rxn_smarts = get_rxn_smarts_yaml(yaml_dict[rxn_name])
        rxn_enzyme_cofactor = get_enzyme_cofactor_yaml(yaml_dict[rxn_name])
        positive_tests = get_positive_tests(yaml_dict[rxn_name])
        negative_tests = get_negative_tests(yaml_dict[rxn_name])
        reaction_type = yaml_dict[rxn_name]['type']
        rxn = Reaction.objects(name=reaction_name)[0]
        experimental = rxn.experimental
        two_step = rxn.two_step
        requires_absence_of_water = rxn.requires_absence_of_water

    result = {"rxn_name": rxn_name,
              "rxn_smarts": rxn_smarts,
              "rxn_enzyme_cofactor": rxn_enzyme_cofactor,
              "positive_tests": positive_tests,
              "negative_tests": negative_tests,
              "reaction_type": reaction_type,
              "experimental": experimental,
              "two_step": two_step,
              "requires_absence_of_water": requires_absence_of_water}

    return jsonify(result=result)

@bp.route("/_visualise_smarts", methods=["POST"])
def visualise_smarts():
    smarts_yaml = request.form['smarts_yaml']
    try:
        smarts_list = yaml.safe_load(smarts_yaml)
    except:
        result = {'status': 'fail',
                  'msg': 'Could not load yaml'}
        return jsonify(result=result)

    try:
        rxn_list = []
        for sma in smarts_list:
            rxn_list.append(rdchiralReaction(sma))
        list_imgs = rxntosvg(rxn_list, rxnSize=(450, 150))
    except:
        result = {'status': 'fail',
                  'msg': 'Could not load rxn imgs'}
        return jsonify(result=result)

    result = {'status': 'success',
              'msg': '',
              'list_imgs': list_imgs}

    return jsonify(result=result)


@bp.route("/_save_reaction", methods=["POST"])
@roles_required('rxn_rules_admin')
def save_reaction():
    rxn_selection = request.form['rxn_selection']
    rxn_name = request.form['rxn_name']
    smarts_yaml = request.form['smarts_yaml']
    cofactors = request.form['cofactors']
    positive_tests = request.form['positive_tests']
    negative_tests = request.form['negative_tests']
    rxn_type = request.form['rxn_type']
    experimental = request.form['experimental']
    two_step = request.form['two_step']
    requires_absence_of_water = request.form['requires_absence_of_water']

    if experimental == 'false':
        experimental = False
    elif experimental == 'true':
        experimental = True

    if two_step == 'false':
        two_step = False
    elif two_step == 'true':
        two_step = True

    if requires_absence_of_water == 'false':
        requires_absence_of_water = False
    elif requires_absence_of_water == 'true':
        requires_absence_of_water = True

    if rxn_selection == 'Empty template':
        if len(Reaction.objects(name=rxn_name)) == 0 and rxn_name != '':
            reaction = Reaction(name=rxn_name)
        else:
            result = {'status': 'failed',
                      'msg': "Could not create new reaction - name already exists",
                      'issues' : []}
            return jsonify(result=result)
    else:
        reaction = Reaction.objects(name=rxn_selection)[0]


    if len(Reaction.objects(name=rxn_name)) == 0:
        reaction.name = rxn_name
    reaction.smarts = yaml.load(smarts_yaml, Loader=yaml.FullLoader)
    reaction.enzyme_types = list(yaml.load(cofactors, Loader=yaml.FullLoader).keys())
    reaction.cofactors = yaml.load(cofactors, Loader=yaml.FullLoader)
    reaction.positive_tests = yaml.load(positive_tests, Loader=yaml.FullLoader)
    reaction.negative_tests = yaml.load(negative_tests, Loader=yaml.FullLoader)
    reaction.type = rxn_type
    reaction.experimental = experimental
    reaction.two_step = two_step
    reaction.requires_absence_of_water = requires_absence_of_water

    reaction.save()

    refresh = 'False'
    if rxn_selection != reaction.name:
        refresh = 'True'
        activities = Activity.objects(reaction=rxn_selection)
        for act in activities:
            act.reaction = reaction.name
            act.save()

    print("reaction saved")
    result = {'status': 'info',
              'msg': "Reaction saved",
              'issues': [],
              'refresh': refresh}
    return jsonify(result=result)

@bp.route("/_test_reaction", methods=["POST"])
@roles_required('rxn_rules_admin')
def test_reaction():
    rxn_selection = request.form['rxn_selection']
    rxn_name = request.form['rxn_name']
    smarts_yaml = request.form['smarts_yaml']
    cofactors = request.form['cofactors']
    positive_tests = request.form['positive_tests']
    negative_tests = request.form['negative_tests']
    rxn_type = request.form['rxn_type']
    experimental = request.form['experimental']
    two_step = request.form['experimental']
    requires_absence_of_water = request.form['requires_absence_of_water']

    tester = ReactionTester()
    tester.run(rxn_selection, rxn_name, smarts_yaml, cofactors,
               positive_tests, negative_tests, rxn_type, experimental, two_step, requires_absence_of_water)

    result = {'status': tester.state,
              'msg': tester.get_msg(),
              'issues': tester.issues}
    return jsonify(result=result)

@bp.route("/_test_product_against_rules", methods=["POST"])
def test_product_against_rules():
    target_smiles = request.form['target_smiles']
    smarts = request.form['smarts']

    try:
        smarts_list = yaml.load(smarts, Loader=yaml.FullLoader)
    except:
        return jsonify(result={'status': 'fail'})

    try:
        rxn_list = []
        for sma in smarts_list:
            rxn_list.append(rdchiralReaction(sma))
    except:
        return jsonify(result={'status': 'fail'})

    try:
        network = Network()

        if request.form['combine_enantiomers'] == 'true':
            network.settings['combine_enantiomers'] = True
        else:
            network.settings['combine_enantiomers'] = False

        network.rxn_obj.rxns = {'test_smarts': rxn_list}
        network.rxns = {'test_smarts': rxn_list}
        network.generate(target_smiles, 1)
        nodes, edges = network.get_visjs_nodes_and_edges()
        print(nodes)
        result = {'status': 'success',
                  'nodes': nodes,
                  'edges': edges}
    except:
        result = {'status': 'fail'}

    return jsonify(result=result)


