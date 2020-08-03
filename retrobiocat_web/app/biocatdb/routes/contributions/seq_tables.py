from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.app.app import user_datastore
import mongoengine as db


def get_enzyme_data(query):
    seq_fields = ['id', 'enzyme_type', 'enzyme_name', 'sequence', 'sequence_unavailable', 'accession', 'structure', 'mutant_of', 'notes', 'papers', 'owner', 'other_names']
    enzyme_data = list(Sequence.objects(query).only(*seq_fields).order_by('enzyme_type').as_pymongo())
    owners_dict = {}
    for i, data in enumerate(enzyme_data):
        enzyme_data[i]['_id'] = str(enzyme_data[i]['_id'])
        enzyme_data[i]['sequence_unavailable'] = str(enzyme_data[i]['sequence_unavailable']).replace('False', '')
        enzyme_data[i]['structure'] = str(enzyme_data[i]['structure']).replace('False', '')

        if 'papers' not in enzyme_data[i]:
            enzyme_data[i]['papers'] = 0
        else:
            enzyme_data[i]['papers'] = len(enzyme_data[i]['papers'])

        if 'owner' in enzyme_data[i]:
            owner_id = str(enzyme_data[i]['owner'])
            if owner_id not in owners_dict:
                owner = User.objects(id=enzyme_data[i]['owner'])[0]
                owners_dict[owner_id] = f"{owner.first_name} {owner.last_name}"
            enzyme_data[i]['owner'] = owners_dict[owner_id]
        else:
            enzyme_data[i]['owner'] = ''

    return enzyme_data

@bp.route('/my_sequences', methods=['GET', 'POST'])
@roles_required('contributor')
def my_sequences():
    user = user_datastore.get_user(current_user.id)

    query = db.Q(owner=user)
    enzyme_data = get_enzyme_data(query)
    enzyme_types = sorted(list(EnzymeType.objects().distinct("enzyme_type")))

    return render_template('edit_tables/edit_sequences.html',
                           seq_data=enzyme_data, seq_button_columns=['edit', 'delete', 'papers'],
                           seq_table_height='80vh', enzyme_types=enzyme_types, show_header_filters=True, include_owner=True)

@bp.route('/edit_sequences', methods=['GET', 'POST'])
@roles_required('super_contributor')
def edit_sequences():
    query = db.Q()
    enzyme_data = get_enzyme_data(query)
    enzyme_types = sorted(list(EnzymeType.objects().distinct("enzyme_type")))

    return render_template('edit_tables/edit_sequences.html',
                           seq_data=enzyme_data, seq_button_columns=['edit', 'merge', 'delete', 'papers'],
                           seq_table_height='80vh', enzyme_types=enzyme_types, show_header_filters=True, include_owner=True)