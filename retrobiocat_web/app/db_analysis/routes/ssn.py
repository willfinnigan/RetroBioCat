from retrobiocat_web.app.db_analysis import bp
from flask import render_template, request, jsonify, session, current_app
from flask_security import roles_required
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, UniRef50, SSN_record, Sequence
import mongoengine as db
from rq.job import Job
from rq import get_current_job
from retrobiocat_web.analysis.make_ssn import SSN, SSN_Visualiser, SSN_quickload
from retrobiocat_web.app.db_analysis.forms import SSN_Form
from retrobiocat_web.analysis import retrieve_uniref_info
import json
from distutils.util import strtobool
import datetime

@bp.route('/ssn_page/<task_id>/', methods=['GET'])
def ssn_page(task_id):
    task = current_app.network_queue.fetch_job(task_id)
    result = task.result
    node_one = result['nodes'][0]
    start_pos = {'x': node_one['x'], 'y': node_one['y']}

    return render_template('ssn/ssn.html',
                           nodes=result['nodes'],
                           edges=result['edges'],
                           alignment_score=result['alignment_score'],
                           start_pos=start_pos,
                           enzyme_type=result['enzyme_type'],
                           hide_mutants=result['hide_mutants'],
                           only_biocatdb=result['only_biocatdb'])

def task_get_ssn(enzyme_type, score, hide_mutants, only_biocatdb):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    ssn = SSN(enzyme_type)

    job.meta['progress'] = 'ssn loaded'
    job.save_meta()

    if only_biocatdb is False and str(score) in ssn.db_object.precalculated_vis:
        nodes = ssn.db_object.precalculated_vis[str(score)]
        # Need to filter mutants here
    else:
        ssn.load(include_mutants=not hide_mutants, only_biocatdb=only_biocatdb)
        vis = SSN_Visualiser(enzyme_type, log_level=1)
        nodes, edges = vis.visualise(ssn, score)

    edges = []

    result = {'nodes': nodes,
              'edges': edges,
              'alignment_score': score,
              'enzyme_type': enzyme_type,
              'hide_mutants': hide_mutants,
              'only_biocatdb': only_biocatdb}

    return result

@bp.route("/ssn_status/<task_id>", methods=["GET"])
def ssn_status(task_id):
    task = current_app.network_queue.fetch_job(task_id)
    task_id = task.get_id()
    task_status = task.get_status()
    seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
    print(seconds_since_active)
    if seconds_since_active > 600:
        print('Job no longer active')
        task_status = 'failed'

    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task_id,
                "task_status": task_status,
                "task_progress": progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route('/ssn_form', methods=['GET', 'POST'])
def ssn_form():
    form = SSN_Form()
    form.set_choices()

    task_id = ''

    if 'ssn_task_id' in session:
        old_task_id = session['ssn_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        enzyme_type = form.data['enzyme_type']
        min_score = form.data['alignment_score']
        hide_mutants = form.data['hide_mutants']
        only_biocatdb = form.data['only_biocatdb']
        task = current_app.network_queue.enqueue(task_get_ssn, enzyme_type, min_score, hide_mutants, only_biocatdb)

        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
            except:
                pass

        task_id = task.get_id()
        session['ssn_task_id'] = task_id

    return render_template('ssn/ssn_form.html', form=form, task_id=task_id)

@bp.route("/_ssn_object", methods=["POST"])
def ssn_object():
    enzyme_type = request.form['enzyme_type']
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    ssn_obj = SSN_record.objects(enzyme_type=enzyme_type_obj)[0]

    num_biocatdb = Sequence.objects(enzyme_type=enzyme_type).count()
    num_uniref = UniRef50.objects(enzyme_type=enzyme_type_obj).count()

    precalc_choices = {}
    for score in ssn_obj.num_at_alignment_score:
        clusters = ssn_obj.num_at_alignment_score[score]
        idt = ssn_obj.identity_at_alignment_score[score]

        choice_text = f"{score}, {clusters} clusters, avg identity {idt[0]} ± {idt[1]}"
        precalc_choices[score] = choice_text

    result = {'status': ssn_obj.status,
              'num_biocatdb': num_biocatdb,
              'num_uniref': num_uniref,
              'precalculated': precalc_choices}
    return jsonify(result=result)


@bp.route("/_load_uniref_data", methods=["POST"])
def load_uniref_data():
    name = request.form['name']
    enzyme_type = request.form['enzyme_type']
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    et = db.Q(enzyme_type=enzyme_type_obj)
    nq = db.Q(enzyme_name=name)

    query = UniRef50.objects(et & nq)
    seq = query[0]
    protein_name = seq.protein_name
    organism = seq.tax

    uniprot_id = retrieve_uniref_info.strip_uniref_name(name)
    if uniprot_id[0:2] == 'UP':
        uniprot_id = ""

    ref_parser = retrieve_uniref_info.UniRef_Parser()
    ref_parser.load_xml(name)
    uni90, uni100, uniprot = ref_parser.get_uniref_members()
    cluster_id = ref_parser.get_cluster_name()
    num_uni90 = len(uni90)
    num_uni100 = len(uni100)
    num_uniprot = len(list(uniprot.keys()))

    if uniprot_id != "":
        prot_parser = retrieve_uniref_info.UniProt_Parser()
        prot_parser.load_xml(uniprot_id)
        pfams = prot_parser.get_pfams()
    else:
        pfams = []

    result = {'rep_seq_name': protein_name,
              'rep_seq_organism': organism,
              'rep_seq_uniprot_id': uniprot_id,
              'cluster_id': cluster_id,
              'num_uni90': num_uni90,
              'num_uni100': num_uni100,
              'num_uniprot': num_uniprot,
              'pfam_object': pfams}
    return jsonify(result=result)

@bp.route("/_edge_ajax", methods=["POST"])
def edge_ajax():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    selected_node = request.form['selected_node']

    ql = SSN_quickload(enzyme_type, log_level=0)
    ql.load_df()
    edges = ql.get_edges(selected_node, alignment_score)

    result = {'edges': edges}
    return jsonify(result=result)

@bp.route("/_connected_nodes_ajax", methods=["POST"])
def connected_nodes_ajax():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    selected_nodes = json.loads(request.form['selected_nodes'])

    ql = SSN_quickload(enzyme_type, log_level=1)
    ql.load_df()
    nodes = ql.get_connected_nodes(selected_nodes, alignment_score)

    result = {'nodes': nodes}
    return jsonify(result=result)

@bp.route("/_get_clusters", methods=["POST"])
def get_clusters():
    enzyme_type = request.form['enzyme_type']
    alignment_score = int(request.form['alignment_score'])
    only_biocatdb = bool(strtobool(request.form['only_biocatdb']))
    hide_mutants = bool(strtobool(request.form['hide_mutants']))

    ssn = SSN(enzyme_type)
    ssn.load(include_mutants=not hide_mutants, only_biocatdb=only_biocatdb)
    without_uniref, with_uniref = ssn.get_clusters(alignment_score)

    result = {'with_uniref': with_uniref,
              'without_uniref': without_uniref}
    return jsonify(result=result)

