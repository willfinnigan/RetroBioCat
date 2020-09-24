from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.app.biocatdb.model_forms import PaperInfo
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Paper, EnzymeType
from retrobiocat_web.app.biocatdb.functions.papers import papers_functions, papers_crossref
from retrobiocat_web.app.app import user_datastore
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField
import mongoengine as db

@bp.route('/paper_progress', methods=['GET', 'POST'])
def paper_progress():

    enzyme_types = EnzymeType.objects().order_by('enzyme_type')
    enzyme_type_progress_dict = {}

    for enz_type_obj in enzyme_types:
        enz_type = enz_type_obj.enzyme_type
        num_papers = len(Paper.objects(tags=enz_type))
        num_complete_papers = len(Paper.objects(db.Q(tags=enz_type) & (db.Q(status='Complete') | db.Q(status='Complete - Awaiting review'))))

        if num_papers == 0:
            pass
        elif num_complete_papers == 0:
            enzyme_type_progress_dict[enz_type] = ["0%", f"{num_complete_papers} out of {num_papers}", 'bg-danger', enz_type_obj.full_name]
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
                                                   colour,
                                                   enz_type_obj.full_name]

            print(enzyme_type_progress_dict[enz_type])


    return render_template('paper_progress/paper_progress.html', enzyme_types_dict=enzyme_type_progress_dict)
