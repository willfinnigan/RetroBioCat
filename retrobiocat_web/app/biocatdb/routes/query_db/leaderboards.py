from flask import render_template, jsonify, session, request, redirect, url_for
from retrobiocat_web.app.biocatdb import bp
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions import sequence_table
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Paper, Sequence, Activity
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.app.biocatdb.forms import SequenceSearch
import itertools


@bp.route("/leaderboard", methods=["GET"])
def leaderboard():
    contributor_role = Role.objects(name='contributor')[0]
    contributors = User.objects(roles=contributor_role)

    papers_dict = {}
    sequence_dict = {}
    activity_dict = {}
    for user in contributors:
        username = f"{user.first_name} {user.last_name}, {user.affiliation}"
        num_papers = len(Paper.objects(owner=user))
        num_sequences = len(Sequence.objects(owner=user))
        papers_dict[username] = num_papers
        sequence_dict[username] = num_sequences

    papers_dict = {k: v for k, v in sorted(papers_dict.items(), key=lambda item: item[1], reverse=True)}
    papers_dict = {k: v for k, v in papers_dict.items() if v != 0}
    sequence_dict = {k: v for k, v in sorted(sequence_dict.items(), key=lambda item: item[1], reverse=True)}
    sequence_dict = {k: v for k, v in sequence_dict.items() if v != 0}



    return render_template('leaderboard.html',
                           top_papers=papers_dict,
                           top_sequences=sequence_dict)
