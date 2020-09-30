from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Molecule, Tag, EnzymeType
from retrobiocat_web.analysis import embl_restfull, all_by_all_blast, make_ssn
from rq.registry import StartedJobRegistry
import datetime
import mongoengine as db

@bp.route('/ssn/<enzyme_type>', methods=['GET'])
def ssn(enzyme_type):
    nodes, edges = make_ssn.get_nodes_and_edges(enzyme_type)

    edges_options = {'smooth': False}
    physics_options = {'stabilization': {'enabled': True,
                                         'iterations': 100},
                       "barnesHut": {"gravitationalConstant": -15000,
                                     "springLength": 300,
                                     "springConstant": 0.1,
                                     "damping": 0.2,
                                     "avoidOverlap": 0.5},
                       "maxVelocity": 49,
                       "minVelocity": 0.75,
                       "solver": "barnesHut"}

    interaction_options = {'tooltipDelay': 0,
                           'hideEdgesOnDrag': True}


    return render_template('ssn/ssn.html',
                           nodes=nodes, edges=edges,
                           edges_options=edges_options,
                           physics_options=physics_options,
                           interaction_options=interaction_options)
