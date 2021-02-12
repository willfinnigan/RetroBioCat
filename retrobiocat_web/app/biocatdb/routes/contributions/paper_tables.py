from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.app.biocatdb.model_forms import PaperInfo
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Paper, EnzymeType, Activity, Sequence
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
                           show_owner=True,
                           title="Super contributor access to all papers",
                           row_click_modal=False)

@bp.route('/my_papers', methods=['GET', 'POST'])
@roles_required('contributor')
def my_papers():
    user = user_datastore.get_user(current_user.id)
    papers_data = list(Paper.objects(owner=user).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data)

    papers_button_columns = ['edit']
    if user.has_role('paper_adder'):
        papers_button_columns.append('delete')

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=papers_button_columns,
                           show_owner=True,
                           title=f"Papers assigned to {user.first_name} {user.last_name}",
                           row_click_modal=False)

@bp.route('/enz_champ_papers/<enzyme_type>', methods=['GET'])
@roles_required('enzyme_champion')
def enzyme_champion_papers(enzyme_type):
    user = user_datastore.get_user(current_user.id)
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    if enzyme_type_obj not in user.enzyme_champion:
        flash('No access', 'danger')
        return redirect(url_for('main_site.home'))

    papers_data = list(Paper.objects(tags=enzyme_type).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=['delete', 'edit'],
                           show_owner=True,
                           title=f"Enzyme champion for {enzyme_type} papers",
                           row_click_modal=False)

@bp.route('/papers_need_data', methods=['GET', 'POST'])
@roles_required('contributor')
def papers_that_need_data():
    user = user_datastore.get_user(current_user.id)
    title = "Papers that require curating"

    args = request.args.to_dict()
    q_no_user = db.Q(owner=None)
    q_no_data = db.Q(status__nin=['Complete - Awaiting review', 'Complete'])
    if 'enzyme_type' in args:
        q_e_type = db.Q(tags=args['enzyme_type'])
        title += f" - {args['enzyme_type']}"
    else:
        q_e_type = db.Q()

    papers_data = list(Paper.objects(q_no_user & q_no_data & q_e_type).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data, show_owner=False)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=['self_assign'],
                           show_owner=False,
                           title=title,
                           row_click_modal=False)

@bp.route('/papers_with_orhpan_sequences', methods=['GET', 'POST'])
@roles_required('admin')
def papers_with_orhpan_sequences():
    title = "Papers with orphan sequences"

    activity_enzyme_names = list(set(Activity.objects().distinct('enzyme_name')))
    paper_ids = []
    for name in activity_enzyme_names:
        if len(Sequence.objects(enzyme_name=name)) == 0:
            act = Activity.objects(enzyme_name=name)[0]
            paper_ids.append(act.paper)

    papers_data = list(Paper.objects(id__in=paper_ids).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data, show_owner=False)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=['delete', 'edit'],
                           show_owner=True,
                           title=title,
                           row_click_modal=False)

@bp.route('/priority_papers', methods=['GET'])
def high_importance_papers():
    hi_q = db.Q(high_importance=True)
    assigned_q = db.Q(owner=None)
    hi_papers = Paper.objects(hi_q & assigned_q).select_related()
    enzyme_types = EnzymeType.objects()

    tags = []
    for paper in hi_papers:
        for tag in paper.tags:
            if tag not in tags:
                tags.append(str(tag))
    tags = sorted(tags)

    data_by_tag = {}
    for tag in tags:
        hi_q = db.Q(high_importance=True)
        tag_q = db.Q(tags=tag)

        papers_data = list(Paper.objects(hi_q & tag_q).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
        papers_data = papers_table.process_papers_dict(papers_data, show_owner=False)
        data_by_tag[tag] = papers_data

    enzyme_full_names = {}
    for enz_type in enzyme_types:
        enzyme_full_names[enz_type.enzyme_type] = enz_type.full_name

    return render_template('edit_tables/high_importance_papers.html',
                           data_by_tag=data_by_tag,
                           tags=tags,
                           enzyme_full_names=enzyme_full_names)
