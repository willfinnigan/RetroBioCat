from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, Paper
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.app.app import user_datastore
import json
import numpy as np
from retrobiocat_web.app.biocatdb.functions import papers_functions
from retrobiocat_web.app.biocatdb.functions.activity import check_activity_data, cascade_activity_data
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions import paper_status

def get_activity_data(paper):
    include = ['id', "reaction", "enzyme_name", "substrate_1_smiles", "substrate_2_smiles",
               "product_1_smiles", "temperature", "ph", "solvent", "other_conditions",
               "notes", "reaction_vol", "formulation", "biocat_conc", "kcat", "km",
               "mw", "substrate_1_conc", "substrate_2_conc", "specific_activity",
               "conversion", "conversion_time", "selectivity", "categorical", "binary",
               ]
    q_paper = db.Q(paper=paper)
    q_not_automated = db.Q(auto_generated__ne=True)

    activity_data = list(Activity.objects(q_paper & q_not_automated).only(*include).as_pymongo())
    for i, data in enumerate(activity_data):
        activity_data[i]['_id'] = str(activity_data[i]['_id'])
        if activity_data[i]['binary']:
            activity_data[i]['binary'] = 1
        else:
            activity_data[i]['binary'] = 0

    return activity_data

def get_paper_data(paper, user):
    self_assigned = ''
    other_user = ''
    if paper.owner == user:
        self_assigned = 'checked'
    elif paper.owner is not None and paper.owner != '':
        other_user = 'disabled'

    paper_owner_name = 'None'
    if paper.owner is not None:
        paper_owner_name = f"{paper.owner.first_name} {paper.owner.last_name}, {paper.owner.affiliation}"

    paper_dict = {'short_cit': paper.short_citation,
                  'doi': paper.doi,
                  'date': paper.date,
                  'title': paper.title,
                  'journal': paper.journal,
                  'authors': papers_functions.list_to_string(paper.authors),
                  'tags': papers_functions.list_to_string(paper.tags),
                  'self_assigned': self_assigned,
                  'disable_for_other_user': other_user,
                  'id': paper.id,
                  'paper_owner_name': paper_owner_name}

    return paper_dict

def get_enzyme_data(paper):
    seq_fields = ['id', 'enzyme_type', 'enzyme_name', 'sequence', 'sequence_unavailable', 'accession', 'structure', 'mutant_of', 'notes', 'papers', 'owner']
    enzyme_data = list(Sequence.objects(papers=paper).only(*seq_fields).as_pymongo())
    owners_dict = {}
    for i, data in enumerate(enzyme_data):
        enzyme_data[i]['_id'] = str(enzyme_data[i]['_id'])
        enzyme_data[i]['sequence_unavailable'] = str(enzyme_data[i]['sequence_unavailable']).replace('False', '')
        enzyme_data[i]['structure'] = str(enzyme_data[i]['structure']).replace('False', '')
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



def get_status(paper):
    paper_progress_text, paper_progress = paper_status.paper_metadata_status(paper)
    sequence_progress_text, sequence_progress = paper_status.sequences_status(paper)
    activity_progress_text, activity_progress = paper_status.activity_status(paper)
    status, status_colour = paper_status.get_status(paper_progress, sequence_progress, activity_progress)

    paper.status = status
    paper.save()

    status_dict = {'review_checked': '',
                   'review_disabled': 'disabled',
                   'status': status,
                   'status_colour': status_colour,
                   'paper_progress': paper_progress,
                   'paper_progress_text': paper_progress_text,
                   'sequences_progress': sequence_progress,
                   'sequences_progress_text': sequence_progress_text,
                   'activity_progress': activity_progress,
                   'activity_progress_text': activity_progress_text}

    if current_user.has_role('super_contributor'):
        status_dict['review_disabled'] = ''

    if paper.reviewed == True:
        status_dict['review_checked'] = 'checked'

    return status_dict

@bp.route('/submission_main_page/<paper_id>', methods=['GET'])
@roles_required('contributor')
def submission_main_page(paper_id):
    user = user_datastore.get_user(current_user.id)
    paper_query = Paper.objects(id=paper_id).select_related()
    if len(paper_query) == 0:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("biocatdb.launch_add_paper"))

    paper = paper_query[0]

    if (paper.owner != user) and (not current_user.has_role('super_contributor')):
        flash('No access to edit this entry', 'fail')
        return redirect(url_for("biocatdb.launch_add_paper"))

    paper_data = get_paper_data(paper, user)
    activity_data = get_activity_data(paper)
    reactions = list(Reaction.objects().distinct('name'))
    enzyme_names = list(Sequence.objects(papers=paper).distinct('enzyme_name'))
    enzyme_types = list(EnzymeType.objects().distinct('enzyme_type'))
    enzyme_data = get_enzyme_data(paper)
    status_dict = get_status(paper)

    return render_template('data_submission/submission_main_page.html',
                           paper=paper_data,
                           activity_data=activity_data,
                           seq_data=enzyme_data, seq_button_columns=['edit', 'remove', 'papers'],
                           status=status_dict,
                           seq_table_height='60vh', enzyme_types=enzyme_types, show_header_filters=False, include_owner=True,
                           reactions=reactions, enzyme_names=enzyme_names+['Chemical'],
                           doi=paper.doi)
