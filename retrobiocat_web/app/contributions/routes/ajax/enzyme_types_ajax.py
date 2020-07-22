from retrobiocat_web.app.contributions import bp
from flask import request, jsonify
from flask_security import roles_required
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity
from retrobiocat_web.mongo.models.reaction_models import Reaction

def change_enzyme_type_name(enz_type, new_name):

    print("Updating enzyme types..")

    old_name = enz_type.enzyme_type
    for seq in Sequence.objects(enzyme_type=old_name):
        seq.enzyme_type = new_name
        seq.save()
    for reaction in Reaction.objects(enzyme_types=old_name):
        reaction.enzyme_types.remove(old_name)
        reaction.enzyme_types.append(new_name)
        reaction.cofactors[new_name] = reaction.cofactors.pop(old_name)
        reaction.save()
    for activity in Activity.objects(enzyme_type=old_name):
        activity.enzyme_type = new_name
        activity.save()

    print('..done')

@bp.route('/_load_enzyme_type_data', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def load_enzyme_type_data():
    name = request.form['enzyme_type']
    enz_type = EnzymeType.objects(enzyme_type=name)[0]

    result = {'name': enz_type.enzyme_type,
              'description': enz_type.description}

    return jsonify(result=result)

@bp.route('/_save_enzyme_type_changes', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def save_enzyme_type_changes():
    original_name = request.form['original_name']
    new_name = request.form['new_name']
    description = request.form['description']

    enz_type = EnzymeType.objects(enzyme_type=original_name)[0]
    enz_type.description = description

    if new_name != original_name:
        if new_name not in list(EnzymeType.objects().distinct('enzyme_type')):
            change_enzyme_type_name(enz_type, new_name)
            enz_type.enzyme_type = new_name
        else:
            result = {'status': 'danger',
                      'msg': 'New name is already taken',
                      'issues': []}
            return jsonify(result=result)

    enz_type.save()

    result = {'status': 'success',
              'msg': 'Changes saved',
              'issues': []}
    return jsonify(result=result)

@bp.route('/_merge_enzyme_type', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def merge_enzyme_types():
    to_merge = request.form['to_merge']
    merge_with = request.form['merge_with']

    if to_merge != merge_with:
        enz_type = EnzymeType.objects(enzyme_type=to_merge)[0]
        change_enzyme_type_name(enz_type, merge_with)
        enz_type.delete()
        result = {'status': 'success',
                  'msg': 'Enzyme type merged',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': "Can't merge with self",
                  'issues': []}

    print(result)

    return jsonify(result=result)

@bp.route('/_delete_enzyme_type', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def delete_enzyme_types():
    to_delete = request.form['to_delete']

    enz_type = EnzymeType.objects(enzyme_type=to_delete)[0]
    seqs = Sequence.objects(enzyme_type=to_delete)
    reacs = Reaction.objects(enzyme_types=to_delete)
    acts = Activity.objects(enzyme_type=to_delete)

    status = 'success'
    msg = 'Enzyme type deleted'
    issues = []

    if len(seqs) != 0:
        status = 'danger'
        msg = 'Could not delete'
        for seq in seqs:
            issues.append(f'Enzyme type is present in sequence: {seq.enzyme_name}')

    if len(reacs) != 0:
        status = 'danger'
        msg = 'Could not delete'
        for reac in reacs:
            issues.append(f'Enzyme type is present in reaction: {reac.name}')

    if len(acts) != 0:
        status = 'danger'
        msg = 'Could not delete'
        papers = []
        for act in acts:
            if act.short_citation not in papers:
                papers.append(act.short_citation)
        for paper in papers:
            issues.append(f"Enzyme type is recorded in activity data for {paper}")


    if status == 'success':
        enz_type.delete()

    result = {'status': status,
              'msg': msg,
              'issues': issues}

    return jsonify(result=result)

