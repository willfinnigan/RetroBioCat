from retrobiocat_web.app.contributions import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.app.contributions.model_forms import PaperInfo
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Paper, EnzymeType
from retrobiocat_web.app.contributions.functions import papers_functions, papers_crossref
from retrobiocat_web.app.app import user_datastore
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField
import mongoengine as db

@bp.route('/paper_progress', methods=['GET', 'POST'])
@roles_required('experimental')
def paper_progress():

    enzyme_types = sorted(list(EnzymeType.objects().distinct('enzyme_type')))
    enzyme_type_progress_dict = {}

    for enz_type in enzyme_types:
        num_papers = len(Paper.objects(tags=enz_type))
        num_complete_papers = len(Paper.objects(db.Q(tags=enz_type) & (db.Q(status='Complete') | db.Q(status='Complete - Awaiting review'))))

        if num_papers == 0:
            pass
        elif num_complete_papers == 0:
            enzyme_type_progress_dict[enz_type] = ["0%", f"{num_complete_papers} out of {num_papers}", 'bg-danger']
        else:
            pc_complete = round((num_complete_papers / num_papers)*100,1)

            if pc_complete > 80:
                colour = 'bg-success'
            elif pc_complete > 40:
                colour = 'bg-warning'
            else:
                colour = 'bg-danger'

            enzyme_type_progress_dict[enz_type] = [f"{pc_complete}%",
                                                   f"{num_complete_papers} out of {num_papers}",
                                                   colour]



    return render_template('paper_progress/paper_progress.html', enzyme_types_dict=enzyme_type_progress_dict)
