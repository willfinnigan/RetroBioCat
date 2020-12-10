from retrobiocat_web.analysis.all_by_all_blast import AllByAllBlaster
from flask import current_app
from retrobiocat_web.analysis.make_ssn import SSN, SSN_Cluster_Precalculator, SSN_Visualiser

def task_expand_ssn(enzyme_type, log_level=1, max_num=200):
    current_app.app_context().push()

    aba_blaster = AllByAllBlaster(enzyme_type, log_level=log_level)
    aba_blaster.make_blast_db()

    ssn = SSN(enzyme_type, aba_blaster=aba_blaster, log_level=log_level)
    ssn.load()
    ssn.set_status('Checking SSN')
    ssn.remove_nonexisting_seqs()
    ssn.remove_seqs_marked_with_no_alignments()

    biocatdb_seqs = ssn.nodes_not_present(only_biocatdb=True, max_num=max_num)
    if len(biocatdb_seqs) != 0:
        ssn.clear_position_information()
        ssn.set_status('Adding and aligning BioCatDB sequences')
        ssn.add_multiple_proteins(biocatdb_seqs)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    need_alignments = ssn.nodes_need_alignments(max_num=max_num)
    if len(need_alignments) != 0:
        ssn.clear_position_information()
        ssn.set_status('Aligning sequences in SSN')
        ssn.add_multiple_proteins(need_alignments)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    not_present = ssn.nodes_not_present(max_num=max_num)
    if len(not_present) != 0:
        ssn.clear_position_information()
        ssn.set_status('Adding UniRef sequences which are not yet present')
        ssn.add_multiple_proteins(not_present)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)

        return

    if ssn.db_object.precalculated_vis == {} and len(ssn.graph.nodes) >= 20:
        ssn.set_status('Precalculating visualisations')
        current_app.preprocess_queue.enqueue(precalculate_job, enzyme_type)

    else:
        ssn.set_status('Complete')
        print(f'- SSN CONSTRUCTION FOR {enzyme_type} IS COMPLETE -')
        ssn.db_object.save()

def new_expand_ssn_job(enzyme_type):
    ssn = SSN(enzyme_type)
    if ssn.db_object.status != 'Complete':
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)

def precalculate_job(enzyme_type):
    ssn = SSN(enzyme_type)
    ssn.load()

    ssn_precalc = SSN_Cluster_Precalculator(ssn)

    num_nodes = len(list(ssn.graph.nodes))
    if num_nodes > 3000:
        num = 1
    elif num_nodes > 1000:
        num = 5
    else:
        num = 20

    if len(list(ssn.db_object.precalculated_vis.keys())) == 0:
        print('No existing % identity data, starting at alignment score 40')
        ssn_precalc.start = 40
        current_num_clusters = 0
    else:
        start_list = [int(s) for s in list(ssn.db_object.identity_at_alignment_score.keys())]
        current_num_clusters = max(list(ssn.db_object.num_at_alignment_score.values()))
        ssn_precalc.start = max(start_list) + 5

    print(f"Start = {ssn_precalc.start}")
    precalculated_nodes, cluster_numbers, identity_at_score = ssn_precalc.precalulate(num=num, current_num_clusters=current_num_clusters)

    if len(precalculated_nodes) == 0:
        print('Precalc complete - checking SSN again')
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)
    else:
        ssn.db_object.precalculated_vis.update(precalculated_nodes)
        ssn.db_object.num_at_alignment_score.update(cluster_numbers)
        ssn.db_object.identity_at_alignment_score.update(identity_at_score)
        ssn.db_object.save()
        current_app.preprocess_queue.enqueue(precalculate_job, enzyme_type)

def remove_sequence(enzyme_type, enzyme_name):
    ssn = SSN(enzyme_type)
    ssn.load()

    if len(list(ssn.graph.nodes)) != 0:
        if enzyme_name in list(ssn.graph.nodes):
            ssn.graph.nodes.remove(enzyme_name)
            ssn.save()
            current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)