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
    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task.get_id(),
                "task_status": task.get_status(),
                "task_progress" : progress
            }
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

@bp.route("/pathway_explorer/<task_id>/", methods=["GET"])
def pathway_explorer(task_id):
    result = json.loads(current_app.redis.get(task_id))

    return render_template('pathway_explorer/pathway_explorer.html',
                           nodes=result['nodes'],
                           edges=result['edges'],
                           options=result['options'],
                           weight_complexity=result['weight_complexity'],
                           weight_num_enzymes=result['weight_num_enzymes'],
                           weight_starting=result['weight_starting'],
                           weight_known_enzymes=result['weight_known_enzymes'],
                           weight_diversity=result['weight_diversity'],
                           max_varient=result['max_varient'],
                           task_id=task_id)

#ajax call used by pathway_explorer
@bp.route('/_next_pathway', methods=['GET', 'POST'])
def next_pathway():
    pathway_num = int(request.form['pathway_num'])
    varient_num = int(request.form['varient_num'])
    task_id = request.form['task_id']
    data = json.loads(current_app.redis.get(task_id))

    graph_dict = json.loads(data['graph_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])
    attr_dict = json.loads(data['attr_dict'])

    all_pathways_data = json.loads(data['pathways_data'])
    pathway_data = all_pathways_data[pathway_num][varient_num]
    max_varient = len(all_pathways_data[pathway_num])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles, print_log=not current_app.config['PRODUCTION'])
    network.update_settings(network_options)
    network.add_attributes(attr_dict)

    pathway = Pathway(pathway_data['nodes'], network, calc_scores=False)

    nodes, edges = pathway.get_visjs_nodes_and_edges()

    result = {'nodes': nodes,
              'edges': edges,
              'max_varient' : max_varient
              }

    return jsonify(result=result)

def task_get_pathways(form_data):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    network = Network(print_log=not current_app.config['PRODUCTION'],
                      include_experimental=form_data['include_experimental'],
                      include_two_step=form_data['include_two_step'])

    network.update_settings({"remove_simple": bool(form_data['remove_small']),
                             "combine_enantiomers": bool(form_data['combine_enantiomers']),
                             'max_nodes': int(form_data['max_nodes']),
                             'similarity_score_threshold': float(form_data['sub_thres']),
                             'colour_substrates': form_data['colour_substrates'],
                             'colour_reactions': form_data['colour_reactions'],
                             "colour_arrows": form_data['colour_edges'],
                             "show_negative_enzymes": form_data['show_neg_enz'],
                             "only_postitive_enzyme_data" : not form_data['show_neg_enz']})

    if form_data["specificity_scoring_mode"] == 'Product + substrates (slower)':
        network.update_settings({'specificity_score_substrates' : True})

    network.generate(form_data['target_smiles'], form_data['number_steps'], calculate_scores=False)

    job.meta['progress'] = 'network_generated'
    job.save_meta()

    network.calculate_scores()

    job.meta['progress'] = 'network_scored'
    job.save_meta()

    bfs = BFS(network=network, max_pathways=form_data['max_pathways'], min_weight=float(form_data['min_weight']), print_log=not current_app.config['PRODUCTION'])
    bfs.run()
    pathways = bfs.get_pathways()

    job.meta['progress'] = 'pathways_generated'
    job.save_meta()

    pathways = group_pathways(pathways)

    pathway_evaluator = PathwayEvaluator(pathways)
    pathway_evaluator.weights = {'Normalised num enzyme steps': form_data['weight_num_enzymes'],
                                 'Normalised change in complexity': form_data['weight_complexity'],
                                 'starting_material': form_data['weight_starting'],
                                 'postitive_enzyme_steps_score': form_data['weight_known_enzymes'],
                                 'Normalised Cofactor Stoichiometry': 0}
    pathway_evaluator.diversity_weight = form_data['weight_diversity']
    pathway_evaluator.evaluate()

    pathways_data = []
    for index, row in pathway_evaluator.df.iterrows():
        varients_data = []
        pathway_varients = [row['Pathway']]
        for pathway in row['Pathway'].other_varients:
            pathway_varients.append(pathway)

        for pathway_var in pathway_varients:
            varients_data.append({'nodes': pathway_var.list_nodes,
                                  'scores_dict': pathway_var.scores.scores_dict()
                                  })
        pathways_data.append(varients_data)

    pathway = Pathway(pathways_data[0][0]['nodes'], network, calc_scores=False)

    job.meta['progress'] = 'pathways_scored'
    job.save_meta()

    nodes, edges = pathway.get_visjs_nodes_and_edges()

    #options = {'physics': {'enabled': 1, 'wind': {'x': 3, 'y': 0}}}
    options = {}

    if form_data['hierarchical'] == True:
        options.update({'layout': {'hierarchical': {'direction': 'directionInput'}}})

    network.clear_cofactor_data()

    result =  {
                'nodes':nodes,
                'edges':edges,
                'options':options,
                'graph_dict':json.dumps(nx.to_dict_of_lists(network.graph)),
                'target_smiles':str(network.target_smiles),
                'network_options':json.dumps(network.settings),
                'max_pathways':len(pathways_data),
                'pathways_data':json.dumps(pathways_data),
                'weight_complexity':form_data['weight_complexity'],
                'weight_num_enzymes':form_data['weight_num_enzymes'],
                'weight_starting':form_data['weight_starting'],
                'weight_known_enzymes':form_data['weight_known_enzymes'],
                'weight_diversity' : form_data['weight_diversity'],
                'attr_dict':json.dumps(network.attributes_dict()),
                'max_varient':len(pathways_data[0]),
                'substrate_nodes':pathway.substrates}

    current_app.redis.mset({job.id: json.dumps(result)})
    time_to_expire = 1*24*60*60  # 1 days * 24 hours * 60 mins * 60 seconds
    current_app.redis.expire(job.id, time_to_expire)
