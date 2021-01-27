from flask import render_template, jsonify, request, current_app, session
from retrobiocat_web.app.retrobiocat import bp, forms
from rq import get_current_job
import networkx as nx
import json
from flask import current_app
from rq.job import Job
from retrobiocat_web.retro.generation.network_generation.network import Network
from retrobiocat_web.retro.generation.pathway_generation.pathway import Pathway
from retrobiocat_web.retro.generation.pathway_generation.best_first_search import BFS
from retrobiocat_web.retro.generation.pathway_generation.pathway_scoring import PathwayEvaluator
from retrobiocat_web.retro.generation.pathway_generation.group_pathways import group_pathways
import datetime


@bp.route('/pathway_explorer_form', methods=['GET', 'POST'])
def pathway_explorer_form():
    form = forms.PathwayExploreForm()
    task_id = ''

    if 'pathway_task_id' in session:
        old_task_id = session['pathway_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        form_data = form.data
        task = current_app.pathway_queue.enqueue(task_get_pathways, form_data)
        if old_task_id != None:
            try:
                old_job = Job.fetch(old_task_id, connection=current_app.redis)
                old_job.delete()
            except:
                pass

        task_id = task.get_id()
        session['pathway_task_id'] = task_id

    return render_template('pathway_explorer_form/pathway_explorer_form.html', form=form, task_id=task_id)

@bp.route("/pathway_explorer_status/<task_id>", methods=["GET"])
def pathway_explorer_status(task_id):
    task = current_app.pathway_queue.fetch_job(task_id)
    task_id = task.get_id()
    task_status = task.get_status()
    seconds_since_active = 0
    try:
        seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
        print(seconds_since_active)
    except:
        pass
    if seconds_since_active > 180 and task_status != 'finished':
        print('Job no longer active')
        print(task_status)
        task_status = 'failed'
    if seconds_since_active > 600:
        print('Job no longer active')
        print(task_status)
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
            }
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route("/pathway_explorer/<task_id>/", methods=["GET"])
def pathway_explorer(task_id):
    #nodes, edges, max_varient = get_visjs_pathway(task_id, 1, 1)
    pathway_settings = json.loads(current_app.redis.get(task_id + '__pathway_settings'))
    pathway_data = json.loads(current_app.redis.get(f"{task_id}__1"))
    nodes, edges, max_varient = pathway_data[0]

    return render_template('pathway_explorer/pathway_explorer.html',
                           nodes=nodes,
                           edges=edges,
                           max_varient=max_varient,
                           options=pathway_settings['options'],
                           weight_complexity=pathway_settings['weight_complexity'],
                           weight_num_enzymes=pathway_settings['weight_num_enzymes'],
                           weight_starting=pathway_settings['weight_starting'],
                           weight_known_enzymes=pathway_settings['weight_known_enzymes'],
                           weight_diversity=pathway_settings['weight_diversity'],
                           task_id=task_id)

#ajax call used by pathway_explorer
@bp.route('/_next_pathway', methods=['GET', 'POST'])
def next_pathway():
    pathway_num = int(request.form['pathway_num'])
    varient_num = int(request.form['varient_num'])
    task_id = request.form['task_id']

    pathway_data = json.loads(current_app.redis.get(f"{task_id}__{pathway_num}"))
    nodes, edges, max_varient = pathway_data[varient_num-1]
    current_app.redis.expire(f"{task_id}__network", 60 * 60)

    result = {'nodes': nodes,
              'edges': edges,
              'max_varient': max_varient
              }

    return jsonify(result=result)

def task_get_pathways(form_data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    network = Network(print_log=not current_app.config['PRODUCTION'],
                      include_experimental=form_data['include_experimental'],
                      include_two_step=form_data['include_two_step'],
                      include_requires_absence_of_water=bool(form_data['include_requires_absence_of_water']))

    network.update_settings({"remove_simple": bool(form_data['remove_small']),
                             "combine_enantiomers": bool(form_data['combine_enantiomers']),
                             'max_nodes': int(form_data['max_nodes']),
                             'similarity_score_threshold': float(form_data['sub_thres']),
                             'colour_reactions': form_data['colour_reactions'],
                             "colour_arrows": form_data['colour_edges'],
                             "show_negative_enzymes": form_data['show_neg_enz'],
                             "only_postitive_enzyme_data": not form_data['show_neg_enz'],
                             'only_reviewed_activity_data': bool(form_data["only_reviewed"])})

    if form_data["specificity_scoring_mode"] == 'Product + substrates (slower)':
        network.update_settings({'specificity_score_substrates': True})

    network.generate(form_data['target_smiles'], form_data['number_steps'], calculate_scores=False)

    job.meta['progress'] = 'network_generated'
    job.save_meta()

    network.calculate_scores()

    job.meta['progress'] = 'network_scored'
    job.save_meta()

    network_data = {'graph_dict': json.dumps(nx.to_dict_of_lists(network.graph)),
                    'target_smiles': str(network.target_smiles),
                    'network_options': json.dumps(network.settings),
                    'attr_dict': json.dumps(network.attributes_dict())}

    current_app.redis.mset({f"{job.id}__network": json.dumps(network_data)})
    current_app.redis.expire(f"{job.id}__network", 60 * 60)

    bfs = BFS(network=network,
              max_pathways=form_data['max_pathways'],
              max_pathway_length=form_data['number_steps'],
              min_weight=float(form_data['min_weight']),
              print_log=not current_app.config['PRODUCTION'])
    bfs.run()
    pathways = bfs.get_pathways()

    job.meta['progress'] = 'pathways_generated'
    job.save_meta()

    package_all_pathways(job.id, pathways)

    pathway_evaluator = evaluate_pathways(pathways, [form_data['weight_num_enzymes'],
                                                     form_data['weight_complexity'],
                                                     form_data['weight_starting'],
                                                     form_data['weight_known_enzymes'],
                                                     form_data['weight_diversity']])

    package_evaluated_pathways(pathway_evaluator.df, job.id)
    package_visjs_pathways(job.id)

    job.meta['progress'] = 'pathways_scored'
    job.save_meta()

    options = {}

    if form_data['hierarchical'] == True:
        options.update({"layout": {"improvedLayout": 'true',
                                   'hierarchical': {'direction': 'DU',
                                                    "sortMethod": "hubsize",
                                                    "nodeSpacing": 200,
                                                    "treeSpacing": 400}}})

    pathway_settings = {'weight_num_enzymes': form_data['weight_num_enzymes'],
                        'weight_complexity': form_data['weight_complexity'],
                        'weight_starting': form_data['weight_starting'],
                        'weight_known_enzymes': form_data['weight_known_enzymes'],
                        'weight_diversity': form_data['weight_diversity'],
                        'options': options}
    current_app.redis.mset({f"{job.id}__pathway_settings": json.dumps(pathway_settings)})
    current_app.redis.expire(job.id, 60 * 60)

def evaluate_pathways(pathways, weights):
    pathways = group_pathways(pathways)
    pathway_evaluator = PathwayEvaluator(pathways)
    pathway_evaluator.weights = {'Normalised num enzyme steps': weights[0],
                                 'Normalised change in complexity': weights[1],
                                 'starting_material': weights[2],
                                 'postitive_enzyme_steps_score': weights[3],
                                 'Normalised Cofactor Stoichiometry': 0}
    pathway_evaluator.diversity_weight = weights[4]
    pathway_evaluator.evaluate()
    return pathway_evaluator



def package_evaluated_pathways(pathway_evaluator_df, task_id):
    evaluated_pathways = []
    for index, row in pathway_evaluator_df.iterrows():
        pathway_data = [row['Pathway'].list_nodes]
        for pathway in row['Pathway'].other_varients:
            pathway_data.append(pathway.list_nodes)
        evaluated_pathways.append(pathway_data)

    current_app.redis.mset({f"{task_id}__evaluated_pathways": json.dumps(evaluated_pathways)})
    current_app.redis.expire(f"{task_id}__evaluated_pathways", 60 * 60)

def package_all_pathways(task_id, pathways):
    all_pathways = []
    all_scores = []
    for pathway in pathways:
        all_pathways.append(pathway.list_nodes)
        all_scores.append(pathway.scores.scores_dict())

    current_app.redis.mset({f"{task_id}__all_pathways": json.dumps((all_pathways, all_scores))})
    current_app.redis.expire(f"{task_id}__all_pathways", 60 * 60)

def package_visjs_pathways(task_id, max_vis=100):
    network_data = json.loads(current_app.redis.get(task_id + '__network'))
    graph = nx.from_dict_of_lists(json.loads(network_data['graph_dict']), create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=network_data['target_smiles'],
                      print_log=not current_app.config['PRODUCTION'])
    network.update_settings(json.loads(network_data['network_options']))
    network.add_attributes(json.loads(network_data['attr_dict']))

    evaluated_pathways = json.loads(current_app.redis.get(f"{task_id}__evaluated_pathways"))

    for i, pathway_varients in enumerate(evaluated_pathways):
        if i > max_vis:
            break
        pathway_vis_js_data = []
        max_var = len(pathway_varients)
        for nodes in pathway_varients:
            pathway = Pathway(nodes, network, calc_scores=False)
            nodes, edges = pathway.get_visjs_nodes_and_edges()
            pathway_vis_js_data.append((nodes, edges, max_var))
        current_app.redis.mset({f"{task_id}__{i+1}": json.dumps(pathway_vis_js_data)})
        current_app.redis.expire(f"{task_id}__{i+1}", 60 * 60)

def get_visjs_pathway(task_id, pathway_id, varient):
    network_data = json.loads(current_app.redis.get(task_id + '__network'))
    pathway_data = json.loads(current_app.redis.get(task_id + f'__{pathway_id}'))
    pathway_nodes = pathway_data[varient-1]

    graph = nx.from_dict_of_lists(json.loads(network_data['graph_dict']), create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=network_data['target_smiles'], print_log=not current_app.config['PRODUCTION'])
    network.update_settings(json.loads(network_data['network_options']))
    network.add_attributes(json.loads(network_data['attr_dict']))

    pathway = Pathway(pathway_nodes, network, calc_scores=False)

    nodes, edges = pathway.get_visjs_nodes_and_edges()
    max_var = len(pathway_data)

    return nodes, edges, max_var

