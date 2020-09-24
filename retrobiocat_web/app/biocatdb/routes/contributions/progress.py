from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.mongo.models.biocatdb_models import Paper, EnzymeType
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions import calculate_paper_progress


@bp.route('/paper_progress', methods=['GET', 'POST'])
def paper_progress():

    enzyme_types = EnzymeType.objects().order_by('enzyme_type')
    enzyme_type_progress_dict = {}

    for enz_type_obj in enzyme_types:
        enz_type = enz_type_obj.enzyme_type
        enzyme_type_progress_dict[enz_type] = calculate_paper_progress.get_enzyme_paper_progress(enz_type_obj)

    return render_template('paper_progress/paper_progress.html', enzyme_types_dict=enzyme_type_progress_dict)
