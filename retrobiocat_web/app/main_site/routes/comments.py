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

@bp.route('/_submit_comment', methods=['GET', 'POST'])
@auth_required()
def submit_comment():
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
    else:
        user = None

    parent_type = request.form['parent_type']
    parent_id = request.form['parent_id']
    comment_text = request.form['comment']
    comment_id = request.form['comment_id']

    print(comment_id)

    if comment_id != '':
        comment_obj = Comment.objects(id=comment_id)[0]
        if current_user.has_role('rxn_rules_admin') or comment_obj.owner == user:
            comment_obj.text = comment_text
            comment_obj.save()
        else:
            result = {'status': 'danger',
                      'msg': 'Could not edit comment - no access',
                      'issues': []}
            return jsonify(result=result)

    else:
        comment_obj = Comment(owner=user,
                              text=comment_text)
        comment_obj.save()

        if parent_type == 'issue':
            issue = Issue.objects(id=parent_id)[0]
            issue.comments.append(comment_obj)
            issue.save()


    result = {'status': 'success',
              'msg': 'Comment submitted',
              'issues': []}
    return jsonify(result=result)


@bp.route('/_delete_comment', methods=['GET', 'POST'])
@auth_required()
def delete_comment():
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
    else:
        user = None

    comment_id = request.form['comment_id']
    comment_obj = Comment.objects(id=comment_id)[0]
    if comment_obj.owner == user or current_user.has_role('rxn_rules_admin'):
        comment_obj.delete()

        result = {'status': 'success',
                  'msg': 'Comment deleted',
                  'issues': []}
        flash('Comment deleted', 'success')
        return jsonify(result=result)

    else:

        result = {'status': 'danger',
                  'msg': 'Could not delete - You are not the owner of this comment and do not have access',
                  'issues': []}
        flash('You are not the owner of this comment and do not have access', 'danger')
        return jsonify(result=result)


@bp.route('/_load_comment_info', methods=['GET', 'POST'])
def load_comment_info():
    comment_id = request.form['comment_id']

    if comment_id != '':
        comment_obj = Comment.objects(id=comment_id)[0]
        comment_txt = comment_obj.text
    else:
        comment_txt = ''

    result = {'comment_text': comment_txt}
    return jsonify(result=result)

