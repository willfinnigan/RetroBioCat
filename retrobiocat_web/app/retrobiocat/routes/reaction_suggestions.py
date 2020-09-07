from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
from retrobiocat_web.app.retrobiocat import bp, forms
from flask import render_template, jsonify, session, request, make_response, flash, redirect, url_for
from retrobiocat_web.mongo.models.reaction_models import Issue
import mongoengine as db
from flask_security import roles_required, current_user, auth_required
from retrobiocat_web.app.app import user_datastore
import json
from distutils.util import strtobool
from retrobiocat_web.mongo.models.reaction_models import Issue, Reaction, Comment

@bp.route('/suggest_new_reaction/')
def suggest_new_reaction():
    page_title = "Suggest a new reaction"
    suggestion_id = ""
    reaction_name = ""
    reaction_smarts = ""
    why_this_reaction = ""
    comments = []
    can_edit_suggestion = True

    return render_template('reaction_suggestions/suggest_new_reaction.html',
                           page_title=page_title,
                           suggestion_id=suggestion_id,
                           reaction_name=reaction_name,
                           reaction_smarts=reaction_smarts,
                           why_this_reaction=why_this_reaction,
                           comments=comments,
                           can_edit_suggestion=can_edit_suggestion)
