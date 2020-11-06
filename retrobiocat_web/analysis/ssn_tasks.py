from retrobiocat_web.analysis.all_by_all_blast import AllByAllBlaster
from flask import current_app
from retrobiocat_web.analysis.make_ssn import SSN, SSN_Cluster_Precalculator

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
        ssn.set_status('Adding and aligning BioCatDB sequences')
        ssn.clear_position_information()
        ssn.add_multiple_proteins(biocatdb_seqs)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    need_alignments = ssn.nodes_need_alignments(max_num=max_num)
    if len(need_alignments) != 0:
        ssn.set_status('Aligning sequences in SSN')
        ssn.clear_position_information()
        ssn.add_multiple_proteins(need_alignments)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)
        return

    not_present = ssn.nodes_not_present(max_num=max_num)
    if len(not_present) != 0:
        ssn.set_status('Adding UniRef sequences which are not yet present')
        ssn.clear_position_information()
        ssn.add_multiple_proteins(not_present)
        ssn.save()
        current_app.alignment_queue.enqueue(new_expand_ssn_job, enzyme_type)

        return

    if ssn.db_object.identity_at_alignment_score == {}:
        ssn.set_status('Precalculating identity at alignment')
        current_app.preprocess_queue.enqueue(new_precalculate_identity_at_alignment_job, enzyme_type)

    elif ssn.db_object.pos_at_alignment_score == {}:
        ssn.set_status('Precalculating cluster positions')
        current_app.preprocess_queue.enqueue(new_precalculate_job, enzyme_type)

    else:
        ssn.set_status('Complete')
        print(f'- SSN CONSTRUCTION FOR {enzyme_type} IS COMPLETE -')
        ssn.save()

def new_expand_ssn_job(enzyme_type):
    ssn = SSN(enzyme_type)
    if ssn.db_object.status != 'Complete':
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)

def new_precalculate_job(enzyme_type):
    ssn = SSN(enzyme_type)
    ssn.load()
    ssn_precalc = SSN_Cluster_Precalculator(ssn)

    num_nodes = len(list(ssn.graph.nodes))
    if num_nodes > 7500:
        num = 1
    elif num_nodes > 5000:
        num = 4
    else:
        num = 20

    if len(list(ssn.db_object.num_at_alignment_score.keys())) == 0:
        ssn_precalc.start = 10
        current_num_clusters = 0
    else:
        start_list = [int(s) for s in list(ssn.db_object.num_at_alignment_score.keys())]
        current_num_clusters = max(list(ssn.db_object.num_at_alignment_score.values()))
        ssn_precalc.start = max(start_list) + 5

    num_at_alignment_score, pos_at_alignment_score = ssn_precalc.precalulate(num=num, current_num_clusters=current_num_clusters)

    if num_at_alignment_score == {}:
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)

    else:
        ssn.db_object.pos_at_alignment_score.update(pos_at_alignment_score)
        ssn.db_object.num_at_alignment_score.update(num_at_alignment_score)
        ssn.db_object.save()
        current_app.preprocess_queue.enqueue(new_precalculate_job, enzyme_type)

def new_precalculate_identity_at_alignment_job(enzyme_type):
    ssn = SSN(enzyme_type)
    ssn.load()

    ssn_precalc = SSN_Cluster_Precalculator(ssn)

    num_nodes = len(list(ssn.graph.nodes))
    if num_nodes > 7500:
        num = 1
    elif num_nodes > 5000:
        num = 5
    else:
        num = 20

    if len(list(ssn.db_object.identity_at_alignment_score.keys())) == 0:
        print('No existing % identity data, starting at alignment score 10')
        ssn_precalc.start = 10
    else:
        start_list = [int(s) for s in list(ssn.db_object.identity_at_alignment_score.keys())]
        ssn_precalc.start = max(start_list) + 5

    identity_at_alignment_score = ssn_precalc.precalculate_identity_at_alignment(num=num)

    if identity_at_alignment_score == {}:
        current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)
    else:
        ssn.db_object.identity_at_alignment_score.update(identity_at_alignment_score)
        ssn.db_object.save()
        current_app.preprocess_queue.enqueue(new_precalculate_identity_at_alignment_job, enzyme_type)


def remove_sequence(enzyme_type, enzyme_name):
    ssn = SSN(enzyme_type)
    ssn.load()

    if len(list(ssn.graph.nodes)) != 0:
        if enzyme_name in list(ssn.graph.nodes):
            ssn.graph.nodes.remove(enzyme_name)
            ssn.save()
            current_app.alignment_queue.enqueue(task_expand_ssn, enzyme_type)