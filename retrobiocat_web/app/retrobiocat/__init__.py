from flask import Blueprint

bp = Blueprint('retrobiocat',
               __name__,
               template_folder='templates',
               static_folder='static',
               static_url_path='/retrobiocat/static')

from retrobiocat_web.app.retrobiocat.routes.network_explorer import make_network, get_network, \
    save_network, add_steps, change_network_options, keep_session_open
from retrobiocat_web.app.retrobiocat.routes.pathway_explorer import pathway, pathway_options_and_reorder
from retrobiocat_web.app.retrobiocat.routes.reactions import reaction_routes
from retrobiocat_web.app.retrobiocat.routes import get_node_information, list_network_saves, node_modal_info
from retrobiocat_web.app.retrobiocat.routes import reaction_issues, reaction_suggestions
