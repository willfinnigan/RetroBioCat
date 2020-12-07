from flask import Blueprint

bp = Blueprint('db_analysis',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/db_analysis/static'
               )

from retrobiocat_web.app.db_analysis.routes import bioinformatics, ssn

