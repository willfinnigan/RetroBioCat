from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.app.app import user_datastore
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions import sequence_table


@bp.route('/my_sequences', methods=['GET', 'POST'])
@roles_required('contributor')
def my_sequences():
    user = user_datastore.get_user(current_user.id)

    query = db.Q(owner=user)
    enzyme_data = sequence_table.get_enzyme_data(query)
    enzyme_types = sorted(list(EnzymeType.objects().distinct("enzyme_type")))

    return render_template('edit_tables/edit_sequences.html',
                           seq_data=enzyme_data, seq_button_columns=['edit', 'delete', 'papers'],
                           seq_table_height='80vh', enzyme_types=enzyme_types, show_header_filters=True, include_owner=True)

@bp.route('/edit_sequences', methods=['GET', 'POST'])
@roles_required('super_contributor')
def edit_sequences():
    query = db.Q()
    enzyme_data = sequence_table.get_enzyme_data(query)
    enzyme_types = sorted(list(EnzymeType.objects().distinct("enzyme_type")))

    return render_template('edit_tables/edit_sequences.html',
                           seq_data=enzyme_data, seq_button_columns=['edit', 'merge', 'delete', 'papers'],
                           seq_table_height='80vh', enzyme_types=enzyme_types, show_header_filters=True, include_owner=True)