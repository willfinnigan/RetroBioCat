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

    result = {'enzyme_name': activity.enzyme_name,
              'enzyme_type': activity.enzyme_type,
              'reaction': activity.reaction,
              'sub1_img': sub1,
              'sub2_img': sub2,
              'prod_img': prod,
              'added_by': added_by}

    print(result)

    return jsonify(result=result)

