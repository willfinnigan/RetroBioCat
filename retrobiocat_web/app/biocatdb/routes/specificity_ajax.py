from flask import render_template, jsonify, session, current_app, request
from retrobiocat_web.app.biocatdb import bp
from retrobiocat_web.mongo.models.biocatdb_models import Activity, Sequence
from retrobiocat_web.app.biocatdb.functions import images


@bp.route('/_load_single_activity_data', methods=['GET', 'POST'])
def load_single_activity_data():
    activity_id = request.form['activity_id']

    activity = Activity.objects(id=activity_id)[0].select_related()

    seq = Sequence.objects(enzyme_name=activity.enzyme_name)

    sub1 = images.smitosvg_url(activity.substrate_1_smiles)

    if activity.product_1_smiles == '' or activity.product_1_smiles is None:
        prod = ''
    else:
        prod = images.smitosvg_url(activity.product_1_smiles)

    if activity.substrate_2_smiles == '' or activity.substrate_2_smiles is None:
        sub2 = ''
    else:
        sub2 = images.smitosvg_url(activity.substrate_2_smiles)

    if activity.added_by is not None:
        added_by = f"{activity.added_by.first_name} {activity.added_by.last_name}, {activity.added_by.affiliation}"
    else:
        added_by = ''

    if activity.substrate_1_conc is None:
        sub1_conc = ''
    else:
        sub1_conc = activity.substrate_1_conc
        if sub1_conc.replace('.', '').isnumeric():
            sub1_conc += ' mM'

    if activity.substrate_2_conc is None:
        sub2_conc = ''
    else:
        sub2_conc = activity.substrate_2_conc
        if sub2_conc.replace('.', '').isnumeric():
            sub2_conc += ' mM'

    if activity.biocat_conc is None:
        biocat_conc = ''
    else:
        biocat_conc = activity.biocat_conc
        if biocat_conc.replace('.', '').isnumeric():
            biocat_conc += ' mg/ml'

    result = {'enzyme_name': activity.enzyme_name,
              'enzyme_type': activity.enzyme_type,
              'reaction': activity.reaction,
              'sub1_img': sub1,
              'sub2_img': sub2,
              'prod_img': prod,
              'added_by': added_by,
              'sub1_conc': sub1_conc,
              'sub2_conc': sub2_conc,
              'biocat_conc': biocat_conc,
              'formulation': activity.formulation}

    for key in result:
        if result[key] is None:
            result[key] = ''

    print(result)

    return jsonify(result=result)

