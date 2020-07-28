from flask import render_template, jsonify, session
from retrobiocat_web.app.biocatdb.forms import SubstrateForm
from retrobiocat_web.app.biocatdb.routes.specificity_query_task import task_run_query
from flask import current_app
from rq.job import Job
from retrobiocat_web.retro.enzyme_identification import query_mongodb
from retrobiocat_web.app.biocatdb import bp
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.retro.enzyme_identification import molecular_similarity
import pandas as pd
from rdkit.Chem import PandasTools
import os
from rdkit import RDConfig
from rdkit import Chem
from requests.utils import quote

from retrobiocat_web.app.biocatdb import bp
from rq import get_current_job

from rdkit.Chem import rdDepictor
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

import mongoengine as db
from retrobiocat_web.mongo.models.biocatdb_models import Activity

import numpy as np

COLUMNS = ["Reaction",
           "Enzyme type",
           "Enzyme name",
           "Substrate 1 SMILES",
           "Substrate 2 SMILES",
           "Product 1 SMILES",
           "Temperature",
           "pH",
           "Solvent",
           "Other conditions",
           "Notes",
           "Reaction volume (ml)",
           "Biocatalyst Formulation",
           "Biocatalyst Concentration (mg/ml)",
           "kcat (min-1)",
           "KM (mM)",
           "Enz MW (Da)",
           "Substrate 1 conc (mM)",
           "Substrate 2 conc (mM)",
           "Specific activity (U/mg)",
           "Conversion (%)",
           "Conversion time (hrs)",
           "Selectivity",
           "Categorical",
           "Binary",
           'Data source',
           'Data source doi',
           'paper',
           'activity_id']


@bp.route('/substrate_specificity_form',  methods=['GET', 'POST'])
def substrate_specificity_form():
    print('Spec form')
    form = SubstrateForm()
    enzymes = query_mongodb.get_enzymes_in_db() + ['All']
    reactions = query_mongodb.get_reactions_in_db()
    task_id = ''

    if 'specificity_task_id' in session:
        old_task_id = session['specificity_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        form_data = form.data
        task = current_app.task_queue.enqueue(get_spec_data, form_data)
        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
                print('Deleted task - ' + str(old_task_id))
            except:
                print('No task to delete - ' + str(old_task_id))

        task_id = task.get_id()
        session['specificity_task_id'] = task_id
        print(task_id)

    return render_template('substrate_specificity/substrate_specificity_form.html', form=form, enzymes=enzymes, reactions=reactions, task_id=task_id)


@bp.route("/substrate_specificity_form_status/<task_id>", methods=["GET"])
def substrate_specificity_status(task_id):
    task = current_app.task_queue.fetch_job(task_id)
    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task.get_id(),
                "task_status": task.get_status(),
                "task_progress" : progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route("/substrate_specificity/<task_id>/", methods=["GET"])
def substrate_specificity(task_id):
    task = current_app.task_queue.fetch_job(task_id)
    activity_data = task.result

    return render_template('spec_2/spec_2.html', activity_data=activity_data)


def get_spec_data(form_data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()
    print('Started')

    enzyme_names = list(form_data['enzymes'].split(", "))
    reactions = list(form_data['reactions'].split(", "))

    if form_data['target_smiles'] != '':
        product = rdkit_smile(form_data['target_smiles'])
    else:
        product = form_data['target_smiles']

    similarity_cutoff = form_data['similarity']
    num_choices = form_data['num_choices']
    data_level = form_data['data_level']
    max_hits = form_data['max_hits']
    include_auto_data = bool(form_data['auto_data'])

    scorer = molecular_similarity.SubstrateSpecificityScorer(print_log=False)

    activity_df = scorer.querySpecificityDf(product, reactions, enzyme_names,
                                   dataLevel=data_level,
                                   numEnzymes=num_choices,
                                   simCutoff=similarity_cutoff,
                                   numHits=max_hits,
                                   include_auto_generated=include_auto_data)

    if activity_df is None:
        return []

    if len(activity_df.index) != 0:
        return []

    activity_df = activity_df[COLUMNS]
    activity_df = activity_df.round(2)
    activity_df.replace(np.nan, '', inplace=True)
    activity_df.replace(True, 'True', inplace=True)
    activity_df.replace(False, 'False', inplace=True)

    activity_data = activity_df.to_dict(orient='records')
    activity_data = process_activity_data(activity_data)
    activity_data = smiles_to_svg(activity_data)

    return activity_data

def process_activity_data(activity_data):
    for i, record in enumerate(activity_data):
        activity_data[i]['paper'] = str(activity_data[i]['paper'])
        activity_data[i]['activity_id'] = str(activity_data[i]['activity_id'])

    return activity_data

def moltosvg_url(mol,molSize=(150,150),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    url = "data:image/svg+xml;charset=utf-8," + quote(svg)

    return url

def smiles_to_svg(activity_data):
    for i, record in enumerate(activity_data):
        if activity_data[i]["Substrate 1 SMILES"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["Substrate 1 SMILES"])
            url = moltosvg_url(mol)
            activity_data[i]["Substrate 1 SMILES"] = url

        if activity_data[i]["Substrate 2 SMILES"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["Substrate 2 SMILES"])
            url = moltosvg_url(mol)
            activity_data[i]["Substrate 2 SMILES"] = url

        if activity_data[i]["Product 1 SMILES"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["Product 1 SMILES"])
            url = moltosvg_url(mol)
            activity_data[i]["Product 1 SMILES"] = url
    return activity_data

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    scorer = molecular_similarity.SubstrateSpecificityScorer(print_log=False)

    test_product = ''
    reactions = []
    enzymes = ['CAR']

    activity_df = scorer.querySpecificityDf(test_product, reactions, enzymes)

    if len(activity_df.index) != 0:
        activity_data = []
    else:
        #activity_df = activity_df[COLUMNS]
        activity_df = activity_df.round(2)
        activity_df.replace(np.nan, '', inplace=True)
        activity_df.replace(True, 'True', inplace=True)
        activity_df.replace(False, 'False', inplace=True)

        activity_data = activity_df.to_dict(orient='records')

    print(activity_data)

    print(Activity.objects(id='5f1723fd5e413e263902cf17'))