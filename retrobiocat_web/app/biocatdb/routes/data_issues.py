from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
from retrobiocat_web.app.biocatdb import bp, forms
from flask import render_template, jsonify, session, request, make_response, flash, redirect, url_for
import mongoengine as db
from flask_security import roles_required, current_user, auth_required
from retrobiocat_web.app.app import user_datastore
import json
from distutils.util import strtobool
from retrobiocat_web.mongo.models.biocatdb_models import ActivityIssue, Activity
from retrobiocat_web.mongo.models.comments import Comment

# Reaction issues table
# Clicking on it loads the activity page
# Display open issues on activity page.  Or display all issues.

def generate_rxn_svg(activity, rxnSize=(300, 70)):
    reaction_smiles = ""
    if activity.substrate_1_smiles is not None:
        reaction_smiles += f"{activity.substrate_1_smiles}"
    if activity.substrate_2_smiles is not None:
        reaction_smiles += f".{activity.substrate_2_smiles}>>"
    if activity.product_1_smiles is not None:
        reaction_smiles += f"{activity.product_1_smiles}"

    reaction_svg = smiles_rxn_to_svg(reaction_smiles, rxnSize=rxnSize)
    return reaction_svg

def get_issue_data(issues):
    issues_data = []
    for issue in issues:
        data = {}
        data['_id'] = str(issue.id)
        data['date'] = str(issue.date)
        data['paper'] = str(issue.activity.paper.short_citation)
        data['reaction_svg'] = generate_rxn_svg(issue.activity)
        data['raised_by'] = f"{issue.raised_by.first_name} {issue.raised_by.last_name}, {issue.raised_by.affiliation}"
        data['status'] = issue.status
        if len(issue.comments) == 0:
            data['first_comment'] = ''
        else:
            data['first_comment'] = issue.comments[0].text
        issues_data.append(data)

    return issues_data

@bp.route('/activity_data_issues_table')
def activity_data_issues_table():
    data_issues = ActivityIssue.objects().only('id', 'activity', 'date', 'raised_by', 'status', 'comments').select_related()
    issues_data = get_issue_data(data_issues)

    return render_template('issues/activity_issues_table.html', activity_data_issues_data=issues_data)


@bp.route('/_raise_activity_data_issue', methods=['GET', 'POST'])
def raise_activity_data_issue():

    user = user_datastore.get_user(current_user.id)
    activity_id = request.form['activity_id']
    comment = request.form['comment']
    activity = Activity.objects(id=activity_id)[0].select_related()

    comment_obj = Comment(owner=user,
                          text=comment)
    comment_obj.save()

    if len(ActivityIssue.objects(activity=activity)) == 0:
        new_issue = ActivityIssue(activity=activity,
                                  raised_by=user,
                                  comments=[comment_obj])
        new_issue.save()
        result = {'status': 'success',
                  'msg': 'Issue has been created',
                  'issues': []}
        return jsonify(result=result)

    else:
        result = {'status': 'danger',
                  'msg': 'An issue has already been raised for this data',
                  'issues': ['Please see the activity data issues page for more information']}
        return jsonify(result=result)


@bp.route('/activity_data_issue/<issue_id>')
def activity_data_issue(issue_id):
    print(issue_id)
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
    else:
        user = None
    issue = ActivityIssue.objects(id=issue_id).select_related()[0]

    reaction_svg = generate_rxn_svg(issue.activity, rxnSize=(600, 150))

    paper = issue.activity.paper.short_citation

    comments = []
    for comment in issue.comments:
        can_edit = False
        can_delete = False
        if current_user.has_role('rxn_rules_admin') or comment.owner == user:
            can_edit = True
            can_delete = True

        new_comment = {'user': f"{comment.owner.first_name} {comment.owner.last_name}, {comment.owner.affiliation}",
                       'date': comment.date.strftime("%d/%m/%Y, %H:%M:%S"),
                       'comment': comment.text,
                       'comment_id': str(comment.id),
                       'can_edit': can_edit,
                       'can_delete': can_delete
                       }
        comments.append(new_comment)


    return render_template('issues/activity_issue_page.html',
                           reaction_svg=reaction_svg,
                           raised_by=f"{issue.raised_by.first_name} {issue.raised_by.last_name}, {issue.raised_by.affiliation}",
                           date=issue.date.strftime("%d/%m/%Y"),
                           status=issue.status,
                           paper=paper,
                           issue_id=str(issue.id),
                           activity_id=str(issue.activity.id),
                           enzyme_name=issue.activity.enzyme_name,
                           paper_id=str(issue.activity.paper.id),
                           comments=comments)


@bp.route('/_open_close_activity_issue', methods=['GET', 'POST'])
@roles_required('rxn_rules_admin')
def open_close_activity_issue():
    id = request.form['issue_id']
    open_close = request.form['open_close']

    issue = ActivityIssue.objects(id=id)[0]
    if open_close == 'Open':
        issue.status = 'Open'
    elif open_close == 'Close':
        issue.status = 'Closed'
    issue.save()

    flash(f'{open_close} issue complete', 'success')
    result = {'status': 'success',
              'msg': f'{open_close} issue complete',
              'issues': []}
    return jsonify(result=result)





@bp.route('/_delete_activity_issue', methods=['GET', 'POST'])
@roles_required('rxn_rules_admin')
def delete_activity_issue():
    issue_id = request.form['issue_id']

    issue = ActivityIssue.objects(id=issue_id)[0]

    for comment in issue.comments:
        comment.delete()

    issue.delete()

    flash('Issue deleted', 'success')
    result = {'status': 'success',
              'msg': 'Issue deleted',
              'issues': []}
    return jsonify(result=result)


