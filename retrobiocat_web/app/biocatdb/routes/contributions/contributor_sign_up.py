from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.app.biocatdb.model_forms import PaperInfo
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.mongo.models.biocatdb_models import Paper
from retrobiocat_web.app.biocatdb.functions.papers import papers_functions, papers_crossref
from retrobiocat_web.app.app import user_datastore
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField
from retrobiocat_web.app.biocatdb.forms import ContributorSignup

@bp.route('/contributor_sign_up', methods=['GET', 'POST'])
def contributor_sign_up():

    form = ContributorSignup()

    if form.validate_on_submit() == True:
        if current_user.is_authenticated:
            user = user_datastore.get_user(current_user.id)
            contributor_role = Role.objects(name='contributor')[0]
            user_datastore.add_role_to_user(user, contributor_role)
            return redirect(url_for("biocatdb.contributor_is_signed_up"))

    return render_template('contributor_sign_up/contributor_signup.html', form=form)

@bp.route('/contributor_is_signed_up', methods=['GET'])
def contributor_is_signed_up():

    return render_template('contributor_sign_up/contributor_is_signed_up.html')