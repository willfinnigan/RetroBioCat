from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
from retrobiocat_web.app.retrobiocat import bp, forms
from flask import render_template, jsonify, session, request, make_response, flash, redirect, url_for
from retrobiocat_web.mongo.models.reaction_models import Issue
import mongoengine as db
from flask_security import roles_required, current_user, auth_required
from retrobiocat_web.app.app import user_datastore
import json
from distutils.util import strtobool
from retrobiocat_web.mongo.models.reaction_models import Issue, Reaction, ReactionSuggestion
from retrobiocat_web.mongo.models.comments import Comment
import yaml

@bp.route('/suggest_new_reaction/')
def suggest_new_reaction():
    page_title = "Suggest a new reaction"
    suggestion_id = ""
    reaction_name = ""
    reaction_smarts = ""
    reaction_details = ""
    owner = ""
    comments = []
    status = ''
    can_save = False
    can_delete = False
    enable_comments = False
    if current_user.is_authenticated:
        can_save = True

    return render_template('reaction_suggestions/reaction_suggestion.html',
                           page_title=page_title,
                           suggestion_id=suggestion_id,
                           reaction_name=reaction_name,
                           reaction_smarts=reaction_smarts,
                           reaction_details=reaction_details,
                           owner=owner,
                           comments=comments,
                           status=status,
                           can_save=can_save,
                           can_delete=can_delete,
                           enable_comments=enable_comments)


@bp.route('/reaction_suggestion/<suggestion_id>')
def reaction_suggestion(suggestion_id):
    r_sug = ReactionSuggestion.objects(id=suggestion_id).select_related()[0]
    reaction_name = r_sug.name
    page_title = f"Reaction suggestion for {reaction_name}"
    reaction_smarts = str(yaml.dump(list(r_sug.smarts)))
    reaction_details = r_sug.details
    owner = f"{r_sug.owner.first_name} {r_sug.owner.last_name}, {r_sug.owner.affiliation}"
    status = r_sug.status

    can_save = False
    can_delete = False
    enable_comments = True
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
        if r_sug.owner == user or user.has_role('rxn_rules_admin'):
            can_save = True
            can_delete = True
    else:
        user = None

    comments = []
    for comment in r_sug.comments:
        comment_can_edit = False
        comment_can_delete = False
        if current_user.has_role('rxn_rules_admin') or comment.owner == user:
            comment_can_edit = True
            comment_can_delete = True

        new_comment = {'user': f"{comment.owner.first_name} {comment.owner.last_name}, {comment.owner.affiliation}",
                       'date': comment.date.strftime("%d/%m/%Y, %H:%M:%S"),
                       'comment': comment.text,
                       'comment_id': str(comment.id),
                       'can_edit': comment_can_edit,
                       'can_delete': comment_can_delete
                       }
        comments.append(new_comment)


    return render_template('reaction_suggestions/reaction_suggestion.html',
                           page_title=page_title,
                           suggestion_id=suggestion_id,
                           reaction_name=reaction_name,
                           reaction_smarts=reaction_smarts,
                           reaction_details=reaction_details,
                           owner=owner,
                           comments=comments,
                           status=status,
                           can_save=can_save,
                           can_delete=can_delete,
                           enable_comments=enable_comments)



@bp.route('/_save_reaction_suggestion', methods=['GET', 'POST'])
@auth_required()
def save_reaction_suggestion():
    user = user_datastore.get_user(current_user.id)
    suggestion_id = request.form['suggestion_id']
    reaction_name = request.form['reaction_name']
    reaction_smarts = str(request.form['reaction_smarts'])
    details = request.form['details']

    smarts_list = []
    if len(reaction_smarts) != 0:
        try:
            smarts_list = yaml.safe_load(reaction_smarts)
        except:
            result = {'status': 'danger',
                      'msg': 'Could not load SMARTS yaml',
                      'issues': []}
            return jsonify(result=result)

    if suggestion_id == '':

        r_suggestion = ReactionSuggestion(name=reaction_name,
                                          owner=user,
                                          details=details,
                                          smarts=smarts_list,
                                          comments=[])
        r_suggestion.save()

        result = {'status': 'success',
                  'msg': 'Reaction suggestion saved',
                  'issues': [],
                  'suggestion_id': str(r_suggestion.id)}
        return jsonify(result=result)

    else:
        r_suggestion = ReactionSuggestion.objects(id=suggestion_id)[0]
        r_suggestion.name = reaction_name
        r_suggestion.smarts = smarts_list
        r_suggestion.details = details
        r_suggestion.save()

    result = {'status': 'success',
              'msg': 'Reaction suggestion saved',
              'issues': [],
              'suggestion_id': str(r_suggestion.id)}
    return jsonify(result=result)


@bp.route('/_open_close_reaction_suggestion', methods=['GET', 'POST'])
@roles_required('rxn_rules_admin')
def open_close_reaction_suggestion():
    id = request.form['suggestion_id']
    open_close = request.form['open_close']

    r_suggestion = ReactionSuggestion.objects(id=id)[0]
    if open_close == 'Open':
        r_suggestion.status = 'Open'
    elif open_close == 'Close':
        r_suggestion.status = 'Closed'
    r_suggestion.save()

    flash(f'{open_close} issue complete', 'success')
    result = {'status': 'success',
              'msg': f'{open_close} issue complete',
              'issues': []}
    return jsonify(result=result)


@bp.route('/reaction_suggestions_table')
def reaction_suggestions_table():
    suggestions = ReactionSuggestion.objects().only('id', 'date', 'name', 'owner', 'status', 'details').select_related()

    reaction_suggestions_data = []
    for sug in suggestions:
        data = {}
        data['_id'] = str(sug.id)
        data['date'] = str(sug.date)
        data['reaction'] = sug.name
        data['owner'] = f"{sug.owner.first_name} {sug.owner.last_name}, {sug.owner.affiliation}"
        data['status'] = sug.status
        data['details'] = sug.details
        reaction_suggestions_data.append(data)

    return render_template('reaction_suggestions/reactions_suggestions_table.html',
                           reaction_suggestions_data=reaction_suggestions_data)

@bp.route('/_delete_reaction_suggestion', methods=['GET', 'POST'])
def delete_reaction_suggestion():

    suggestion_id = request.form['suggestion_id']
    sug = ReactionSuggestion.objects(id=suggestion_id)[0]

    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
        if sug.owner == user or user.has_rule('rxn_rules_admin'):
            sug.delete()
            for comment in sug.comments:
                comment.delete()

            flash('Issue deleted', 'success')
            result = {'status': 'success',
                      'msg': 'Comment deleted',
                      'issues': []}
            return jsonify(result=result)


    result = {'status': 'danger',
              'msg': 'No access to delete this suggestion',
              'issues': []}
    return jsonify(result=result)








