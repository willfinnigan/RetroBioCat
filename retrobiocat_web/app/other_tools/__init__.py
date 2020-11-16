from flask import Blueprint

bp = Blueprint('other_tools',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/other_tools/static'
               )

from retrobiocat_web.app.other_tools import mutant_generator