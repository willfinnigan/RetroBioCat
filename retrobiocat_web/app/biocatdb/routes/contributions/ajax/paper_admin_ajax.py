from retrobiocat_web.app.biocatdb import bp
from flask import flash, url_for, request, jsonify
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, Activity, Paper, EnzymeType
from retrobiocat_web.app.app import user_datastore
from retrobiocat_web.app.biocatdb.functions.papers import papers_functions, papers_crossref
import mongoengine as db
from distutils.util import strtobool
from retrobiocat_web.app.biocatdb.functions import check_permission


@bp.route('/_admin_set_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_set_owner():
    paper = Paper.objects(id=request.form['paper_id'])[0]
    new_owner_id = request.form['new_owner_id']
    print(new_owner_id)
    new_owner = User.objects(id=new_owner_id)[0]

    paper.owner = new_owner
    paper.save()

    result = {'status': 'success',
              'msg': 'Paper owner updated',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_admin_activity_to_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_activity_to_owner():
    paper = Paper.objects(id=request.form['paper_id']).select_related()[0]
    activities = Activity.objects(paper=paper)

    for activity in activities:
        activity.added_by = paper.owner
        activity.save()

    result = {'status': 'success',
              'msg': 'Activity added by updated',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_admin_unassigned_seqs_to_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_unassigned_seqs_to_owner():
    paper = Paper.objects(id=request.form['paper_id']).select_related()[0]
    seqs = Sequence.objects(db.Q(papers=paper) & db.Q(owner=None))

    for seq in seqs:
        seq.owner = paper.owner
        seq.added_by = paper.owner
        seq.save()

    result = {'status': 'success',
              'msg': 'Unassigned sequences assigned to paper owner',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_admin_all_seqs_to_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_all_seqs_to_owner():
    paper = Paper.objects(id=request.form['paper_id']).select_related()[0]
    seqs = Sequence.objects(db.Q(papers=paper))

    for seq in seqs:
        seq.owner = paper.owner
        seq.save()

    result = {'status': 'success',
              'msg': 'All sequences assigned to paper owner',
              'issues': []}

    return jsonify(result=result)