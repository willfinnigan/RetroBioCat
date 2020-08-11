from flask import Blueprint

bp = Blueprint('biocatdb',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/biocatdb/static'
               )

from retrobiocat_web.app.biocatdb.routes.contributions import data_submission, seq_tables, paper_tables, add_new_paper, progress
from retrobiocat_web.app.biocatdb.routes.contributions.ajax import sequences_ajax, activity_ajax, papers_ajax, enzyme_types_ajax
from retrobiocat_web.app.biocatdb.routes.contributions.admin import add_or_edit_enzyme_type, download_data, edit_rxn_rules, init_mongodb

from retrobiocat_web.app.biocatdb.routes.query_db import sequence_query, substrate_specificity, substrate_specificity_ajax, paper_query, leaderboards
