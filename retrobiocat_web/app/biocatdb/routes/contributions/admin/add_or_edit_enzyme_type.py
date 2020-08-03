from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, request, jsonify
from retrobiocat_web.app.biocatdb.model_forms import EnzymeTypeForm
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity
from retrobiocat_web.mongo.models.reaction_models import Reaction

@bp.route('/add_enzyme_type', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def add_enzyme_type():
    form = EnzymeTypeForm()

    if form.validate_on_submit() == True:
        if form.description.data == '':
            form.description.data = None
        if form.other_abbreviations.data == '':
            other_abbreviations = None
        else:
            other_abbreviations = form.other_abbreviations.data.split(', ')

        enz_type = EnzymeType(enzyme_type=form.enzyme_type.data,
                              full_name=form.full_name.data,
                              other_abbreviations=other_abbreviations,
                              description=form.description.data)
        enz_type.save()

        flash("Data added successfully", 'success')
        return redirect(url_for('biocatdb.add_enzyme_type'))

    return render_template('enzyme_type/add_enzyme_type.html', form=form)

@bp.route('/edit_enzyme_types', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def edit_enzyme_types():
    headings = ['Name', 'Full name', 'Other abbreviations', 'Description', 'Num rules']
    enzyme_types = EnzymeType.objects().distinct("enzyme_type")
    enzyme_types.sort()

    enz_type_dict_list = EnzymeType.objects().order_by('enzyme_type').as_pymongo()
    renamed_enz_type_dict_list = []

    for enz_type_dict in enz_type_dict_list:
        new_enz_type_dict = {}
        enz_type = enz_type_dict.get('enzyme_type')
        new_enz_type_dict['Name'] = enz_type
        new_enz_type_dict['Full name'] = enz_type_dict.get('full_name', '')
        new_enz_type_dict['Other abbreviations'] = enz_type_dict.get('other_abbreviations', '')
        new_enz_type_dict['Description'] = enz_type_dict.get('description', '')
        new_enz_type_dict['Num rules'] = len(Reaction.objects(enzyme_types=enz_type))
        renamed_enz_type_dict_list.append(new_enz_type_dict)

    return render_template('enzyme_type/edit_enzyme_types.html',
                           headings=headings, rows=renamed_enz_type_dict_list, enzyme_types=enzyme_types)