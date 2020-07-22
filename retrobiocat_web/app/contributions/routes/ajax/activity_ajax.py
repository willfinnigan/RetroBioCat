from retrobiocat_web.app.contributions import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, Paper, Molecule
from retrobiocat_web.app.app import user_datastore
import json
from retrobiocat_web.app.contributions.functions.activity import check_activity_data, cascade_activity_data
import mongoengine as db
import numpy as np
from retrobiocat_web.retro.enzyme_identification import query_mongodb
from retrobiocat_web.retro.enzyme_identification.load import make_fingerprints
from retrobiocat_web.mongo.init_db.make_molecule_db import fp_molecules_to_db

def check_is_float(to_test):
    if type(to_test) != float:
        try:
            to_test = float(to_test)
        except:
            return None
    if np.isnan(to_test) == True:
        return None

    return to_test

def check_is_nan(to_test):
    if str(to_test) == 'NaN' or str(to_test) == '' or str(to_test) == ' ':
        return None
    if type(to_test) != str:
        if isinstance(to_test, (int, float, complex)):
            if np.isnan(to_test) == True:
                return None
            else:
                to_test = str(to_test)
    return to_test

@bp.route('/_save_activity_data', methods=['GET', 'POST'])
@roles_required('contributor')
def save_activity_data():
    user = user_datastore.get_user(current_user.id)
    paper = Paper.objects(id=request.form['paper_id'])[0]

    if (paper.owner != user) and (not current_user.has_role('super_contributor')):
        result = {'status': 'danger',
                  'msg': 'No access to modify this paper',
                  'issues': ['Paper not assigned to user, and not a super_contributor']}
        return jsonify(result=result)

    try:
        data = json.loads(request.form['data'])
        result = do_save(data, user, paper)
    except Exception as e:
        result = {'status': 'danger',
                  'msg': 'There was an error',
                  'issues': [str(e)]}
    return jsonify(result=result)

def do_save(data, user, paper):
    data = [x for x in data if x != {}]
    data = add_row_numbers(data)

    processed_data = []

    issues = []
    for data_dict in data:
        issues += check_activity_data.initial_check_data(data_dict)

        if len(issues) == 0:
            data_dict = cascade_activity_data.set_activities(data_dict)
            issues += check_activity_data.check_all_have_binary(data_dict)
            if len(issues) == 0:
                processed_data.append(data_dict)

    if len(issues) != 0:
        result = {'status': 'danger',
                  'msg': 'failed to save data',
                  'issues': issues}
    else:
        ids = [d['_id'] for d in processed_data if '_id' in d]
        ids = [i for i in ids if i != '']
        delete_other_data(ids, paper)

        for data_dict in processed_data:
            update_activity(data_dict, paper, user)

        current_app.task_queue.enqueue(task_update_fingerprints)

        result = {'status': 'success',
                  'msg': 'Data processed and passed tests',
                  'issues': issues,
                  'data': processed_data}
        flash("Data updated", 'success')

    return result

def task_update_fingerprints():
    # Get a list of smiles in db, and those with fingerprints
    spec_df = query_mongodb.query_specificity_data(['All'], ['All'])
    fp_df = make_fingerprints.make_fingerprint_df_for_mongo(spec_df)
    fp_smiles = list(fp_df['smiles'])
    current_smiles = list(Molecule.objects().distinct('smiles'))

    # Get the smiles which are new
    new_smiles = []
    for smi in fp_smiles:
        if smi not in current_smiles:
            new_smiles.append(smi)

    # Get smiles where fp no longer needed
    old_smiles = []
    for smi in current_smiles:
        if smi not in fp_smiles:
            old_smiles.append(smi)

    # Add new fingerprints to database
    new_fp_df = fp_df[fp_df['smiles'].isin(new_smiles)]
    fp_molecules_to_db(new_fp_df)

    # Delete old fingerprints
    mols = Molecule.objects(smiles__in=old_smiles)
    mols.delete()


def delete_other_data(ids, paper):
    p_Q = db.Q(paper=paper)
    data_id_Q = db.Q(id__nin=ids)
    Activity.objects(p_Q & data_id_Q).delete()

def update_activity(data_dict, paper, user):
    if data_dict.get('_id', '') != '':
        query = Activity.objects(id=data_dict['_id'])
    else:
        query = []

    if len(query) == 0:
        activity = Activity()
        activity.added_by = user
    else:
        activity = query[0]

    if user not in paper.edits_by:
        paper.edits_by.append(user)

    seq_types = {}
    if data_dict.get('enzyme_name') not in seq_types:
        seq = Sequence.objects(enzyme_name=data_dict.get('enzyme_name'))[0]
        seq_types['enzyme_name'] = seq.enzyme_type

    activity.enzyme_type = seq_types['enzyme_name']
    activity.enzyme_name = data_dict.get('enzyme_name')
    activity.reaction = data_dict.get('reaction')
    activity.short_citation = paper.short_citation
    activity.html_doi = paper.html
    activity.paper = paper

    activity.cascade_num = str(data_dict.get('cascade_num', ''))
    activity.substrate_1_smiles = str(data_dict.get('substrate_1_smiles', ''))
    activity.substrate_2_smiles = str(data_dict.get('substrate_2_smiles', ''))
    activity.product_1_smiles = str(data_dict.get('product_1_smiles', ''))
    activity.temperature = str(data_dict.get('temperature', ''))
    activity.ph = str(data_dict.get('ph', ''))

    activity.solvent = str(data_dict.get('solvent', ''))
    activity.other_conditions = str(data_dict.get('other_conditions', ''))
    activity.notes = str(data_dict.get('notes', ''))
    activity.reaction_vol = str(data_dict.get('reaction_vol', ''))
    activity.formulation = str(data_dict.get('formulation', ''))
    activity.biocat_conc = str(data_dict.get('biocat_conc', ''))
    activity.kcat = check_is_float(data_dict.get('kcat', None))
    activity.km = check_is_float(data_dict.get('km', None))
    activity.mw = check_is_float(data_dict.get('mw', None))

    activity.substrate_1_conc = str(data_dict.get('substrate_1_conc', ''))
    activity.substrate_2_conc = str(data_dict.get('substrate_2_conc', ''))

    activity.specific_activity = check_is_float(data_dict.get('specific_activity', None))
    activity.conversion = check_is_float(data_dict.get('conversion', None))
    activity.conversion_time = check_is_float(data_dict.get('conversion_time', None))
    activity.categorical = str(data_dict.get('categorical', None))
    if data_dict.get('binary', None) == 1:
        binary = True
    elif data_dict.get('binary', None) == 0:
        binary = False
    else:
        binary = None
    activity.binary = binary
    activity.selectivity = str(data_dict.get('selectivity', ''))
    activity.auto_generated = False

    activity.save()

def add_row_numbers(data):
    for i, data_dict in enumerate(data):
        data[i]['n'] = i+1
    return data
