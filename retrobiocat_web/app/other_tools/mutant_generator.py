from retrobiocat_web.app.other_tools import bp
from flask import render_template
from retrobiocat_web import __version__ as website_version


@bp.route('/mutant_generator', methods=['GET'])
def mutant_generator():
    return render_template('mutant_generator.html')