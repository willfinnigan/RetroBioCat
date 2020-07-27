from flask import Blueprint

bp = Blueprint('biocatdb',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/biocatdb/static'
               )

from retrobiocat_web.app.biocatdb.routes import spec_2
