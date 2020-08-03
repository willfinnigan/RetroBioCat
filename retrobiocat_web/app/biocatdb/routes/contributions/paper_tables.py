from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.app.biocatdb.model_forms import PaperInfo
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Paper
from retrobiocat_web.app.biocatdb.functions.papers import papers_functions, papers_table
from retrobiocat_web.app.app import user_datastore
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField
import mongoengine as db




@bp.route('/edit_papers', methods=['GET', 'POST'])
@roles_required('super_contributor')
def edit_papers():
    papers_data = list(Paper.objects().only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=['delete', 'edit', 'link'],
                           show_owner=True)

@bp.route('/my_papers', methods=['GET', 'POST'])
@roles_required('contributor')
def my_papers():
    user = user_datastore.get_user(current_user.id)
    papers_data = list(Paper.objects(owner=user).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=['delete', 'edit', 'link'],
                           show_owner=True)

@bp.route('/papers_need_data', methods=['GET', 'POST'])
@roles_required('contributor')
def papers_that_need_data():
    user = user_datastore.get_user(current_user.id)
    q_no_user = db.Q(owner=None)
    q_no_data = db.Q(status__nin=['Complete - Awaiting review', 'Complete'])

    papers_data = list(Paper.objects(q_no_user & q_no_data).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data, show_owner=False)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=['self_assign'],
                           show_owner=False)
