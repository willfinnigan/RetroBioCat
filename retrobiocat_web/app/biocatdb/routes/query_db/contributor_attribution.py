from flask import render_template, jsonify, session, request, redirect, url_for
from retrobiocat_web.app.biocatdb import bp
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions import sequence_table
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Paper, Sequence, Activity
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.app.biocatdb.forms import SequenceSearch
import itertools
from retrobiocat_web.app.biocatdb.functions import calculate_paper_progress


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

@bp.route("/enzyme_teams", methods=["GET"])
def enzyme_teams():
    teams = []
    team_info = {}

    contributor_role = Role.objects(name='contributor')[0]
    contributors = User.objects(roles=contributor_role).select_related()

    for user in contributors:
        user_string = f"{user.first_name} {user.last_name}"
        affiliation_string = f"{user.affiliation}"
        for enzyme_type in user.enzyme_teams:
            if enzyme_type.enzyme_type not in teams:
                teams.append(enzyme_type.enzyme_type)
                team_info[enzyme_type.enzyme_type] = {'full_name': enzyme_type.full_name,
                                                      'enzyme_champions': [],
                                                      'team_members': [],
                                                      'progress': calculate_paper_progress.get_enzyme_paper_progress(enzyme_type)}

            if enzyme_type in user.enzyme_champion:
                team_info[enzyme_type.enzyme_type]['enzyme_champions'].append([user_string, affiliation_string])
            else:
                team_info[enzyme_type.enzyme_type]['team_members'].append([user_string, affiliation_string])

    teams.sort()

    return render_template('enzyme_teams.html', teams=teams, team_info=team_info)
