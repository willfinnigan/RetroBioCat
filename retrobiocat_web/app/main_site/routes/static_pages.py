from retrobiocat_web.app.main_site import bp
from flask import render_template
from retrobiocat_web import __version__ as website_version

@bp.route('/', methods=['GET'])
def home():
    return render_template('home.html',
                           rbc_version=website_version)

@bp.route('/automated_cascade_design', methods=['GET'])
def automated_cascade_design():
    return render_template('automated_cascade_design.html')

@bp.route('/cookie_policy', methods=['GET'])
def cookie_policy():
    return render_template('cookie_policy.html')

@bp.route('/terms', methods=['GET'])
def terms():
    return render_template('terms.html')