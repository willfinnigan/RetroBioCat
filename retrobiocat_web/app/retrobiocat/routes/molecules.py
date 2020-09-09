from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
from retrobiocat_web.app.retrobiocat import bp, forms
from flask import render_template, jsonify, session, request, make_response, flash, redirect, url_for
from retrobiocat_web.mongo.models.reaction_models import Issue
import mongoengine as db
from flask_security import roles_required, current_user, auth_required
from retrobiocat_web.app.app import user_datastore
import json
from distutils.util import strtobool
from retrobiocat_web.mongo.models.reaction_models import Issue, Reaction
from retrobiocat_web.mongo.models.comments import Comment
from retrobiocat_web.mongo.models.user_saves import MyMolecule
from retrobiocat_web.app.retrobiocat.functions import get_images
from rdkit import Chem
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile

@bp.route('/_save_my_molecule', methods=['GET', 'POST'])
@auth_required()
def save_my_molecule():
    user = user_datastore.get_user(current_user.id)
    smiles = request.form['smiles']
    name = request.form['name']

    if rdkit_smile(smiles) is None:
        result = {'status': 'danger',
                  'msg': 'Please enter a valid SMILES',
                  'issues': []}
        return jsonify(result=result)

    mol_query = MyMolecule.objects(db.Q(owner=user) & db.Q(smiles=smiles))

    if len(mol_query) != 0:
        my_mol = mol_query[0]
    else:
        my_mol = MyMolecule(owner=user,
                            smiles=smiles)

        try:
            mol = Chem.MolFromSmiles(smiles)
            my_mol.svg = get_images.moltosvg(mol,molSize=(100,100),kekulize=True)
        except:
            my_mol.svg = ''

    my_mol.name = name
    my_mol.save()

    result = {'status': 'success',
              'msg': 'Molecule saved',
              'issues': []}
    return jsonify(result=result)

@bp.route('/_load_my_molecule_info', methods=['GET', 'POST'])
@auth_required()
def load_my_molecule_info():
    user = user_datastore.get_user(current_user.id)
    smiles = request.form['smiles']

    mol_query = MyMolecule.objects(db.Q(owner=user) & db.Q(smiles=smiles))

    if len(mol_query) != 0:
        my_mol = mol_query[0]
        name = my_mol.name
    else:
        name = ''

    if smiles == '':
        img = ''
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            img = get_images.moltosvg(mol,molSize=(200,200),kekulize=True)
        except:
            img = '<p>Could not generate SMILES - may be an invalid SMILES string</p>'

    result = {'name': name,
              'img': img}
    return jsonify(result=result)

@bp.route('/_load_my_molecules', methods=['GET', 'POST'])
@auth_required()
def load_my_molecules():
    user = user_datastore.get_user(current_user.id)

    mols = MyMolecule.objects(db.Q(owner=user))

    mol_dict = {}
    for mol in mols:
        mol_dict[mol.smiles] = (mol.name, mol.svg)

    result = {'mol_dict': mol_dict}
    return jsonify(result=result)

@bp.route('/_delete_my_molecule', methods=['GET', 'POST'])
@auth_required()
def delete_my_molecule():
    user = user_datastore.get_user(current_user.id)
    smiles = request.form['smiles']

    mols = MyMolecule.objects(db.Q(owner=user) & db.Q(smiles=smiles))

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









