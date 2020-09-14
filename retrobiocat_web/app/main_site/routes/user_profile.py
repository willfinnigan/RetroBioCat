from retrobiocat_web.app.user_model_forms import UserProfileForm
from flask import render_template, redirect, url_for, request, jsonify
from retrobiocat_web.app.main_site import bp
from retrobiocat_web.app.app import user_datastore, db
from flask_security import current_user, login_required

@bp.route('/user_profile', methods=['GET', 'POST'])
@login_required
def user_profile():
    if not current_user.is_authenticated:
        redirect(url_for('main_site.home'))

    form = UserProfileForm()
    user = user_datastore.get_user(current_user.id)

    if form.validate_on_submit() == True:
        user.first_name = form.data['first_name']
        user.last_name = form.data['last_name']
        user.affiliation = form.data['affiliation']
        user.email_opt_in = form.data['email_opt_in']
        user.save()
        msg = "User profile updated"
        return render_template('msg.html', msg=msg)

    else:
        form.first_name.data = user['first_name']
        form.last_name.data = user['last_name']
        form.affiliation.data = user['affiliation']
        form.email_opt_in.data = user['email_opt_in']
        return render_template('user_profile.html', form=form)

@bp.route('/delete_account', methods=['GET', 'POST'])
@login_required
def delete_account():
    delete_check = request.form['delete']
    if delete_check == 'DELETE':
        user = user_datastore.get_user(current_user.id)
        user.delete()
        return jsonify({'result':'done'})
    else:
        return jsonify({'result':'failed'})