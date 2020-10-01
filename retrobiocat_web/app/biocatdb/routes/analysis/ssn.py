from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Molecule, Tag, EnzymeType
from retrobiocat_web.analysis import embl_restfull, all_by_all_blast, make_ssn
from rq.registry import StartedJobRegistry
import datetime
import mongoengine as db

@bp.route('/ssn/<enzyme_type>', methods=['GET'])
@roles_required('admin')
def ssn(enzyme_type):
    identity = 0.7
    nodes, edges = make_ssn.get_nodes_and_edges(enzyme_type, identity)

    edges_options = {'smooth': False}
    physics_options = {'stabilization': {'enabled': True,
                                         'iterations': 50},
                       "repulsion": {
                            "centralGravity": 0.1,
                            "nodeDistance": 500
                        },
                       "maxVelocity": 49,
                       "minVelocity": 0.75,
                       "solver": "repulsion"}

    interaction_options = {'tooltipDelay': 0,
                           'hideEdgesOnDrag': True}


    return render_template('ssn/ssn.html',
                           nodes=nodes, edges=edges,
                           edges_options=edges_options,
                           physics_options=physics_options,
                           interaction_options=interaction_options)
