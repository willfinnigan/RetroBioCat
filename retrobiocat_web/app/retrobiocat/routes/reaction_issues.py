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


def get_issue_data(issues):
    reaction_issues_data = []
    for issue in issues:
        data = {}
        data['_id'] = str(issue.id)
        data['date'] = str(issue.date)
        data['reaction'] = issue.reaction.name
        data['issue_reaction_svg'] = issue.issue_reaction_svg
        data['raised_by'] = f"{issue.raised_by.first_name} {issue.raised_by.last_name}, {issue.raised_by.affiliation}"
        data['status'] = issue.status
        data['public'] = str(issue.public)
        if len(issue.comments) == 0:
            data['first_comment'] = ''
        else:
            data['first_comment'] = issue.comments[0].text
        reaction_issues_data.append(data)

    return reaction_issues_data


@bp.route('/reaction_issues_table')
def reaction_issues_table():
    if current_user.has_role('rxn_rules_admin'):
        is_public = db.Q()
        display_public = 'true'
    else:
        is_public = db.Q(public=True)
        display_public = 'false'
    issues = Issue.objects(is_public).only('id', 'date', 'reaction', 'issue_reaction_svg', 'raised_by', 'status', 'public', 'comments').select_related()
    reaction_issues_data = get_issue_data(issues)

    return render_template('reaction_issues/reactions_issues_table.html',
                            reaction_issues_data=reaction_issues_data,
                            display_public=display_public)

@bp.route('/my_reaction_issues_table')
def my_reaction_issues_table():
    user = user_datastore.get_user(current_user.id)

    issues = Issue.objects(raised_by=user).only('id', 'date', 'reaction', 'issue_reaction_svg', 'raised_by', 'status', 'public', 'comments').select_related()
    reaction_issues_data = get_issue_data(issues)
    display_public = 'true'

    return render_template('reaction_issues/reactions_issues_table.html',
                            reaction_issues_data=reaction_issues_data,
                            display_public=display_public)

@bp.route('/_load_reaction_issue_info', methods=['GET', 'POST'])
def load_reaction_issue_info():
    substrates = json.loads(request.form['parents'])
    products = json.loads(request.form['children'])

    reaction_smiles = f"{substrates[0]}"
    if len(substrates) > 1:
        reaction_smiles += f".{substrates[1]}"
    reaction_smiles += f">>{products[0]}"
    query_reaction_svg = smiles_rxn_to_svg(reaction_smiles, rxnSize=(600, 150))

    result = {'query_reaction_svg': query_reaction_svg}
    return jsonify(result=result)

@bp.route('/_submit_reaction_issue', methods=['GET', 'POST'])
@auth_required()
def submit_reaction_issue():
    substrates = json.loads(request.form['parents'])
    products = json.loads(request.form['children'])
    reaction_name = request.form['reaction']
    comment = request.form['comment']
    public = bool(strtobool(request.form['public']))
    print(public)

    reaction_smiles = f"{substrates[0]}"
    if len(substrates) > 1:
        reaction_smiles += f".{substrates[1]}"
    reaction_smiles += f">>{products[0]}"
    query_reaction_svg = smiles_rxn_to_svg(reaction_smiles, rxnSize=(300, 70))

    user = user_datastore.get_user(current_user.id)

    reaction = Reaction.objects(name=reaction_name)[0]
    comment_obj = Comment(owner=user,
                          text=comment)
    comment_obj.save()

    issue = Issue(reaction=reaction,
                  issue_reaction_smiles=reaction_smiles,
                  issue_reaction_svg=query_reaction_svg,
                  raised_by=user,
                  status='Open',
                  comments=[comment_obj],
                  public=public)
    issue.save()

    print(f"Issue saved for {issue.reaction}")

    result = {'status': 'success',
              'msg': 'Issue raised',
              'issues': []}
    return jsonify(result=result)


@bp.route('/reaction_issue/<issue_id>')
def reaction_issue(issue_id):
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
    else:
        user = None
    issue = Issue.objects(id=issue_id).select_related()[0]

    reaction_svg = smiles_rxn_to_svg(issue.issue_reaction_smiles, rxnSize=(600, 150))
    if issue.public == True:
        public = 'Public'
    else:
        public = 'Not public'

    num_issues = len(Issue.objects(db.Q(reaction=issue.reaction) & db.Q(status='Open')))

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

    return render_template('reaction_issues/issue_page.html',
                           reaction_name=issue.reaction.name,
                           reaction_svg=reaction_svg,
                           reaction_smiles=issue.issue_reaction_smiles,
                           raised_by=f"{issue.raised_by.first_name} {issue.raised_by.last_name}, {issue.raised_by.affiliation}",
                           date=issue.date.strftime("%d/%m/%Y"),
                           status=issue.status,
                           public=public,
                           num_issues=num_issues,
                           issue_id=str(issue.id),
                           comments=comments)


@bp.route('/_open_close_reaction_issue', methods=['GET', 'POST'])
@roles_required('rxn_rules_admin')
def open_close_reaction_issue():
    id = request.form['issue_id']
    open_close = request.form['open_close']

    issue = Issue.objects(id=id)[0]
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

@bp.route('/_delete_reaction_issue', methods=['GET', 'POST'])
@roles_required('rxn_rules_admin')
def delete_reaction_issue():
    issue_id = request.form['issue_id']
    print(issue_id)

    issue = Issue.objects(id=issue_id)[0]

    for comment in issue.comments:
        comment.delete()

    issue.delete()


    flash('Issue deleted', 'success')
    result = {'status': 'success',
              'msg': 'Comment deleted',
              'issues': []}
    return jsonify(result=result)



