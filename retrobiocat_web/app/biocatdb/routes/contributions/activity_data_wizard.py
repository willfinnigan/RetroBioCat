from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, Paper
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.app.app import user_datastore
import json
import numpy as np
from retrobiocat_web.app.biocatdb.functions.papers import papers_functions
from retrobiocat_web.app.biocatdb.functions.activity import check_activity_data, cascade_activity_data
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions.papers import paper_status
from retrobiocat_web.app.biocatdb.functions import check_permission
from retrobiocat_web.app.biocatdb.functions import sequence_table
from retrobiocat_web.app.biocatdb.functions.substrate_specificity import process_activity_data
from retrobiocat_web.app.biocatdb.routes.contributions.data_submission import get_activity_data, get_paper_data
from werkzeug.utils import secure_filename
from rq import get_current_job
from rq.job import Job
import string
from flask_security import roles_required, current_user, auth_required
from natsort import natsorted, ns
from rdkit import Chem
from retrobiocat_web.app.retrobiocat.functions.get_images import moltosvg

from retrobiocat_web.curation import structure_recognition
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.mongo.models.biocatdb_models import ActivityMol

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from retrobiocat_web.retro.generation.retrosynthesis_engine.retrosynthesis_engine import RuleApplicator
from retrobiocat_web.retro.rdchiral.main import rdchiralReaction
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg



def getletterfromindex(num):
    #produces a string from numbers so
    # from https://stackoverflow.com/questions/23199733/convert-numbers-into-corresponding-letter-using-python

    #1->a
    #2->b
    #26->z
    #27->aa
    #28->ab
    #52->az
    #53->ba
    #54->bb

    num2alphadict = dict(zip(range(1, 27), string.ascii_lowercase))
    outval = ""
    numloops = (num-1) //26

    if numloops > 0:
        outval = outval + getletterfromindex(numloops)

    remainder = num % 26
    if remainder > 0:
        outval = outval + num2alphadict[remainder]
    else:
        outval = outval + "z"
    return outval

def group_activity_data(activity_data):
    if len(activity_data) == 0:
        return []

    g = 0
    previous_substrate_1 = None
    previous_substrate_2 = None
    previous_product_1 = None
    previous_reaction = None
    activity_grouped = {}
    for data in activity_data:
        if data["substrate_1_smiles"] != previous_substrate_1 or \
                data["substrate_2_smiles"] != previous_substrate_2 or \
                data["product_1_smiles"] != previous_product_1 or \
                data["reaction"] != previous_reaction:
            g += 1

        if str(g) not in activity_grouped:
            activity_grouped[str(g)] = []
        activity_grouped[str(g)].append(data)

        previous_substrate_1 = data["substrate_1_smiles"]
        previous_substrate_2 = data["substrate_2_smiles"]
        previous_product_1 = data["product_1_smiles"]
        previous_reaction = data['reaction']

    return activity_grouped


@bp.route('/_upload_molecule_images',methods=['GET', 'POST'])
@roles_required('contributor')
def upload_molecule_images():
    issues = []
    saved_files = []
    paper_id = request.form['mol_paper_id_field']
    if request.method != 'POST':
        issues.append('Method is not POST')
    else:
        for name, file in request.files.items():
            filename = secure_filename(file.filename)
            saved = structure_recognition.save_uploaded_image(file, filename)
            if saved == True:
                saved_files.append(filename)

    start_num = int(request.form['current_num_mols']) + 1
    task = current_app.osra_queue.enqueue(task_process_images, saved_files, paper_id, start_num=start_num)
    task_id = task.get_id()

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f'Uploaded files {len(saved_files)} - now processing',
                  'issues': [],
                  'task_id': task_id}
    else:
        result = {'status': 'danger',
                  'msg': 'Problem uploading / processing files',
                  'issues': issues}

    return jsonify(result=result)

def task_process_images(list_filenames, paper_id, start_num=1):
    osra_url = current_app.config['OSRA_API_HOST']
    paper = Paper.objects(id=paper_id)[0]

    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    dict_smiles = {}
    n = start_num
    for i, filename in enumerate(list_filenames):
        new_smis = structure_recognition.process_image(filename, osra_url, log=True)
        for j, smi in enumerate(new_smis):
            name = f"{n}"
            dict_smiles[name] = smi
            n += 1

    for name, smi in dict_smiles.items():
        structure_recognition.save_mol(smi, paper, name=name)

    return dict_smiles


@bp.route("/process_images_status/<task_id>", methods=["GET"])
def process_images_status(task_id):
    task = current_app.osra_queue.fetch_job(task_id)
    task_id = task.get_id()
    task_status = task.get_status()

    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task_status == 'finished':
        result = task.result
    else:
        result = {}

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task_id,
                "task_status": task_status,
                "task_progress": progress,
                "result": result
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202


@bp.route('/_delete_paper_molecule', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_paper_molecule():
    id = request.form['id']
    paper_id = request.form['paper_id']

    paper_query = Paper.objects(id=paper_id).select_related()
    if len(paper_query) != 0:
        paper = paper_query[0]
        if check_permission.check_paper_permission(current_user.id, paper):

            mols = ActivityMol.objects(db.Q(id=id))

            if len(mols) != 0:
                mol = mols[0]
                mol.delete()

                result = {'status': 'success',
                          'msg': 'Molecule deleted',
                          'issues': []}
                return jsonify(result=result)

    result = {'status': 'danger',
              'msg': 'Could not delete molecule',
              'issues': []}
    return jsonify(result=result)

@bp.route('/_duplicate_paper_molecule', methods=['GET', 'POST'])
@roles_required('contributor')
def duplicate_paper_molecule():
    id = request.form['id']
    paper_id = request.form['paper_id']

    paper_query = Paper.objects(id=paper_id).select_related()
    if len(paper_query) != 0:
        paper = paper_query[0]
        if check_permission.check_paper_permission(current_user.id, paper):

            mols = ActivityMol.objects(db.Q(id=id))

            if len(mols) != 0:
                mol = mols[0]
                num_paper_mols = ActivityMol.objects(paper=paper).count()
                new_mol = ActivityMol(name=f"{num_paper_mols+1}",
                                      smi=mol.smi,
                                      svg=mol.svg,
                                      paper=paper)
                new_mol.save()
                mol_list = [(new_mol.name, new_mol.smi, new_mol.svg, str(new_mol.id))]

                result = {'status': 'success',
                          'msg': 'Molecule deleted',
                          'issues': [],
                          'mol_list': mol_list}
                return jsonify(result=result)

    result = {'status': 'danger',
              'msg': 'Could not delete molecule',
              'issues': []}
    return jsonify(result=result)


@bp.route('/_reorder_paper_molecules', methods=['GET', 'POST'])
@roles_required('contributor')
def reorder_paper_molecules():
    paper_id = request.form['paper_id']

    paper_query = Paper.objects(id=paper_id).select_related()
    if len(paper_query) != 0:
        paper = paper_query[0]
        if check_permission.check_paper_permission(current_user.id, paper):
            paper_mols = list(ActivityMol.objects(paper=paper))
            paper_mols = natsorted(paper_mols, key=lambda x: x.name)

            for i, mol in enumerate(list(paper_mols)):
                mol.name = str(i+1)
                mol.save()

            result = {'status': 'success',
                      'msg': 'Molecules reordered',
                      'issues': []}
            return jsonify(result=result)

    result = {'status': 'danger',
              'msg': 'Could not delete molecule',
              'issues': []}
    return jsonify(result=result)

@bp.route('/_update_paper_molecule', methods=['GET', 'POST'])
@roles_required('contributor')
def update_paper_molecule():
    paper_id = request.form['paper_id']
    id = request.form['id']
    smi = request.form['smi']
    name = request.form['name']
    paper_query = Paper.objects(id=paper_id).select_related()
    if len(paper_query) != 0:
        paper = paper_query[0]
        if check_permission.check_paper_permission(current_user.id, paper):
            act_mol = ActivityMol.objects(id=id)[0]

            if smi != act_mol.smi:
                try:
                    mol = Chem.MolFromSmiles(smi)
                    svg = moltosvg(mol)
                except:
                    svg = ""

                act_mol.svg = svg
                act_mol.smi = smi

            act_mol.name = name
            act_mol.save()

            result = {'status': 'success',
                      'msg': 'Molecule saved',
                      'issues': []}
            return jsonify(result=result)

    result = {'status': 'danger',
              'msg': 'Could not delete molecule',
              'issues': []}
    return jsonify(result=result)

def reverse_smarts(smarts):
    m = smarts.find('>>')
    start = smarts[m+2:]
    end = smarts[:m]
    new_smarts = start + '>>' + end
    return new_smarts


def make_svg_dict(products):
    svg_dict = {}
    for list_smis in products:
        for smi in list_smis:
            if smi not in svg_dict:
                mol = Chem.MolFromSmiles(smi)
                svg_dict[smi] = moltosvg(mol)
    return svg_dict

def get_reaction_svg(s1, s2, p):
    if s2 == "":
        smi_rxn = f"{s1}>>{p}"
    else:
        smi_rxn = f"{s1}.{s2}>>{p}"

    svg = smiles_rxn_to_svg(smi_rxn, rxnSize=(500,200))
    return svg


@bp.route('/_apply_reaction_fwd', methods=['GET', 'POST'])
@roles_required('contributor')
def apply_reaction_fwd():
    s1 = request.form['s1']
    s2 = request.form['s2']
    reaction_q = Reaction.objects(name=request.form['reaction'])
    if len(reaction_q) == 0:
        print(f"{request.form['reaction']} was not found in RetroBioCat")
        result = {'status': 'danger',
                  'msg': 'No reaction on RetroBioCat with this name',
                  'issues': [f"{request.form['reaction']} was not found in RetroBioCat"]}
        return jsonify(result=result)

    rdkit_rxns = {'rdkit_rxns': []}
    rdchiral_rxns = {'rd_chiral_rxns': []}
    try:
        for smarts in reaction_q[0].smarts:
            reverse = reverse_smarts(smarts)
            rdkit_rxns['rdkit_rxns'].append(rdChemReactions.ReactionFromSmarts(reverse))
            try:
                rdchiral_rxns['rd_chiral_rxns'].append(rdchiralReaction(reverse))
            except:
                pass
    except Exception as e:
        print('Error parsing reaction SMARTS')
        print(e)
        result = {'status': 'danger',
                  'msg': 'Error parsing reaction SMARTS',
                  'issues': []}
        return jsonify(result=result)

    rule_applicator = RuleApplicator(None)

    smi = f"{s1}"
    if s2 != "":
        smi += f".{s2}"


    try:
        products = rule_applicator.apply_rules(smi, rdchiral_rxns)['rd_chiral_rxns']
    except:
        products = []

    if products == []:
        try:
            products = rule_applicator.apply_rules_rdkit(smi, rdkit_rxns)['rdkit_rxns']
        except:
            products = []


    no_dupl = []
    for prod in products:
        if prod not in no_dupl:
            no_dupl.append(prod)
    products = no_dupl

    svg_dict = make_svg_dict(products)
    reaction_svgs = []
    for smi_list in products:
        reaction_svgs.append(get_reaction_svg(s1, s2, smi_list[0]))

    if len(products) == 0:
        result = {'status': 'danger',
                  'msg': 'No products',
                  'issues': []}
        return jsonify(result=result)

    result = {'status': 'success',
              'msg': '',
              'products': products,
              'reaction_svgs': reaction_svgs,
              'svg_dict': svg_dict,
              'issues': []}
    return jsonify(result=result)

@bp.route('/_apply_reaction_rev', methods=['GET', 'POST'])
@roles_required('contributor')
def apply_reaction_rev():
    p = request.form['p']
    reaction_q = Reaction.objects(name=request.form['reaction'])
    if len(reaction_q) == 0:
        print(f"{request.form['reaction']} was not found in RetroBioCat")
        result = {'status': 'danger',
                  'msg': 'No reaction on RetroBioCat with this name',
                  'issues': [f"{request.form['reaction']} was not found in RetroBioCat"]}
        return jsonify(result=result)

    rdkit_rxns = {'rdkit_rxns': []}
    rdchiral_rxns = {'rd_chiral_rxns': []}
    try:
        for smarts in reaction_q[0].smarts:
            rdkit_rxns['rdkit_rxns'].append(rdChemReactions.ReactionFromSmarts(smarts))
            try:
                rdchiral_rxns['rd_chiral_rxns'].append(rdchiralReaction(smarts))
            except:
                pass
    except Exception as e:
        print('Error parsing reaction SMARTS')
        print(e)
        result = {'status': 'danger',
                  'msg': 'Error parsing reaction SMARTS',
                  'issues': []}
        return jsonify(result=result)

    rule_applicator = RuleApplicator(None)

    smi = p

    try:
        products = rule_applicator.apply_rules(smi, rdchiral_rxns)['rd_chiral_rxns']
    except:
        products = []

    if products == []:
        try:
            products = rule_applicator.apply_rules_rdkit(smi, rdkit_rxns)['rdkit_rxns']
        except:
            products = []

    no_dupl = []
    for prod in products:
        if prod not in no_dupl:
            no_dupl.append(prod)
    products = no_dupl

    svg_dict = make_svg_dict(products)
    reaction_svgs = []
    for smi_list in products:
        if len(smi_list) == 1:
            reaction_svgs.append(get_reaction_svg(smi_list[0], '', p))
        elif len(smi_list) == 2:
            reaction_svgs.append(get_reaction_svg(smi_list[0], smi_list[1], p))
        else:
            reaction_svgs.append(get_reaction_svg('', '', p))

    if len(products) == 0:
        result = {'status': 'danger',
                  'msg': 'No products',
                  'issues': []}
        return jsonify(result=result)

    result = {'status': 'success',
              'msg': '',
              'products': products,
              'reaction_svgs': reaction_svgs,
              'svg_dict': svg_dict,
              'issues': []}

    return jsonify(result=result)





