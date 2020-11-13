from retrobiocat_web.app.main_site import bp
from flask import render_template
from retrobiocat_web import __version__ as website_version
import datetime
__version__ = datetime.date.today().strftime("%d%m%y")

@bp.route('/', methods=['GET'])
def home():
    version = f"{website_version}"
    return render_template('home.html',
                           rbc_version=version)

@bp.route('/automated_cascade_design', methods=['GET'])
def automated_cascade_design():
    return render_template('automated_cascade_design.html')

@bp.route('/cookie_policy', methods=['GET'])
def cookie_policy():
    return render_template('cookie_policy.html')

@bp.route('/terms', methods=['GET'])
def terms():
    return render_template('terms.html')

@bp.route('/retrosynthesis_help', methods=['GET'])
def retrosynthesis_help():
    return render_template('retrosynthesis_help.html')

@bp.route('/instructions_for_data_contribution', methods=['GET'])
def instructions_for_data_contribution():
    return render_template('instructions_for_data_contribution.html')