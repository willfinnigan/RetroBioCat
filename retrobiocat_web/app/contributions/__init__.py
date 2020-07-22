from flask import Blueprint

bp = Blueprint('contributions',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/contributions/static'
               )

from retrobiocat_web.app.contributions.routes import data_submission, seq_tables, paper_tables, add_new_paper, progress
from retrobiocat_web.app.contributions.routes.ajax import sequences_ajax, activity_ajax, papers_ajax, enzyme_types_ajax
from retrobiocat_web.app.contributions.routes.admin import add_or_edit_enzyme_type, download_data, edit_rxn_rules, init_mongodb

