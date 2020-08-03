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

    if activity.temperature is None:
        temperature = ''
    else:
        temperature = activity.temperature
        if temperature.replace('.','').isnumeric():
            temperature += ' <sup>o</sup>C'

    if activity.ph is None:
        ph = ''
    else:
        ph = activity.ph
        if ph.replace('.','').isnumeric():
            ph = 'pH ' + ph

    if activity.reaction_vol is None:
        volume = ''
    else:
        volume = activity.reaction_vol
        if volume.replace('.', '').isnumeric():
            volume += ' mL scale'

    if activity.binary == True:
        active = 'Active'
    else:
        active = 'Not active'

    if activity.specific_activity is None:
        sa = ''
    else:
        sa = str(round(activity.specific_activity,2))
        sa += ' &#956;mol / min / mg'

    if activity.conversion is None:
        conv = ''
    else:
        conv = str(int(activity.conversion)) + " % conversion"
        if activity.conversion_time is not None:
            conv += f" in {str(int(activity.conversion_time))} hours"

    if activity.kcat is None:
        kinetics = ''
    else:
        kinetics = f"kcat: {str(round(activity.kcat,2))} min<sup>-1</sup> km: {str(round(activity.km,2))} mM"

    result = {'short_cit': activity.paper.short_citation,
              'enzyme_name': activity.enzyme_name,
              'enzyme_type': activity.enzyme_type,
              'reaction': activity.reaction,
              'sub1_img': sub1,
              'sub2_img': sub2,
              'prod_img': prod,
              'added_by': added_by,
              'sub1_conc': sub1_conc,
              'sub2_conc': sub2_conc,
              'biocat_conc': biocat_conc,
              'formulation': activity.formulation,
              'selectivity': activity.selectivity,
              'temperature': temperature,
              'ph': ph,
              'solvent': activity.solvent,
              'volume': volume,
              'other_conditions': activity.other_conditions,
              'notes': activity.notes,
              'active': active,
              'category': activity.categorical,
              'sa': sa,
              'conv': conv,
              'kinetics': kinetics}

    for key in result:
        if result[key] is None:
            result[key] = ''
    return jsonify(result=result)

