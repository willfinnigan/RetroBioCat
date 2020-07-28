from retrobiocat_web.app.contributions import bp
from flask import request, jsonify
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, Activity
from retrobiocat_web.app.app import user_datastore
from distutils.util import strtobool
from retrobiocat_web.app.contributions.functions import paper_status

def seqs_of_type(enzyme_type):
    sequences = Sequence.objects(enzyme_type=enzyme_type).distinct('enzyme_name')
    sequences.sort()

    seq_array = {}
    for seq in sequences:
        seq_array[seq] = seq

    result = {'sequences': seq_array}
    return result

@bp.route('/_sequences_of_same_type', methods=['GET', 'POST'])
@roles_required('contributor')
def get_sequences_of_same_type():
    enzyme_type = Sequence.objects(enzyme_name=request.form['enzyme_name'])[0].enzyme_type
    result = seqs_of_type(enzyme_type)

    return jsonify(result=result)

@bp.route('/_sequences_of_type', methods=['GET', 'POST'])
@roles_required('contributor')
def get_sequences_of_type():
    enzyme_type = request.form['enzyme_type']
    result = seqs_of_type(enzyme_type)

    return jsonify(result=result)

@bp.route('/_save_edited_sequence', methods=['GET', 'POST'])
@roles_required('contributor')
def save_edited_sequence():
    original_name = request.form['original_name']
    enzyme_name = request.form['enzyme_name']
    enzyme_type = request.form['enzyme_type']
    sequence = request.form['sequence']
    sequence_unavailable = bool(strtobool(request.form['sequence_unavailable']))
    accession = request.form['accession']
    structure = bool(strtobool(request.form['structure']))
    mutant_of = request.form['mutant_of']
    notes = request.form['notes']
    other_names = request.form['other_names']

    status = 'success'
    msg = 'Sequence edited'
    issues = []

    seq = Sequence.objects(enzyme_name=original_name)[0]
    user = user_datastore.get_user(current_user.id)

    seq.enzyme_type = enzyme_type
    seq.sequence_unavailable = sequence_unavailable
    seq.accession = accession
    seq.structure = structure
    seq.notes = notes
    seq.mutant_of = mutant_of
    seq.other_names = other_names.split(', ')

    if seq.added_by is None:
        seq.added_by = user
    elif user not in seq.edits_by:
        seq.edits_by.append(user)

    if original_name != enzyme_name:
        if len(Sequence.objects(enzyme_name=enzyme_name)) == 0:
            seq.enzyme_name = enzyme_name
        else:
            status = 'danger'
            msg = 'Could not update sequence'
            issues.append('New name is already a sequence')

    self_assigned = bool(strtobool(request.form['self_assigned']))

    if self_assigned == True:
        seq.owner = user
    elif self_assigned == False and seq.owner == user:
        seq.owner = None

    amino_acids_list = ['X', 'V', 'G', 'F', 'E', 'N', 'P', 'Q', 'M', 'K', 'T', 'S', 'W', 'A', 'R', 'D', 'L', 'Y', 'H',
                        'I', 'C', '*']

    sequence = sequence.replace('\n', '')

    aa_check = True
    for letter in sequence:
        if letter.upper() not in amino_acids_list:
            aa_check = False

    if aa_check == False:
        status = 'danger'
        msg = 'Could not update sequence'
        issues.append('Protein sequence contains non amino acid characters')
    else:
        seq.sequence = sequence

    if (seq.owner != None) and (not current_user.has_role('super_contributor')):
        if seq.owner != user:
            status = 'danger'
            msg = 'Could not update sequence'
            issues.append('User does not have access to edit this sequence')

    if status == 'success':
        seq.save()
        update_seq_papers_status(seq.enzyme_name)

    result = {'status': status,
              'msg': msg,
              'issues': issues}
    return jsonify(result=result)

@bp.route('/_load_sequence_data', methods=['GET', 'POST'])
@roles_required('contributor')
def load_sequence_data():
    name = request.form['name']

    seq = Sequence.objects(enzyme_name=name)[0]

    sequences_same_type = Sequence.objects(enzyme_type=seq.enzyme_type).distinct('enzyme_name')
    sequences_same_type.sort()

    seq_array = {}
    for seq_same_type in sequences_same_type:
        seq_array[seq_same_type] = seq_same_type

    user = user_datastore.get_user(current_user.id)
    other_user = False
    self_assigned = False
    can_edit = True
    if seq.owner == user:
        self_assigned = True
    else:
        if seq.owner != '' and seq.owner is not None:
            other_user = True
            if not user.has_role('super_contributor'):
                can_edit = False

    other_names = ''
    for i, name in enumerate(seq.other_names):
        other_names += name
        if (len(seq.other_names) > 1) and (i < len(seq.other_names)-1):
            other_names += ', '


    result = {'enzyme_type': seq.enzyme_type,
              'enzyme_name': seq.enzyme_name,
              'sequence': seq.sequence,
              'sequence_unavailable': seq.sequence_unavailable,
              'accession': seq.accession,
              'structure': seq.structure,
              'mutant_of': seq.mutant_of,
              'sequences': seq_array,
              'notes': seq.notes,
              'can_edit': can_edit,
              'self_assigned': self_assigned,
              'owner_is_another_user': other_user,
              'other_names': other_names}

    return jsonify(result=result)

@bp.route('/_delete_sequence', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_sequence():
    to_delete = request.form['to_delete']

    seq = Sequence.objects(enzyme_name=to_delete)[0]
    acts = Activity.objects(enzyme_name=to_delete)

    status = 'success'
    msg = 'Sequence deleted'
    issues = []

    if len(acts) != 0:
        status = 'danger'
        msg = 'Could not delete'
        papers = []
        for act in acts:
            if act.short_citation not in papers:
                papers.append(act.short_citation)
        for paper in papers:
            issues.append(f"Sequence is recorded in activity data for {paper}")

    mutants = Sequence.objects(mutant_of=to_delete)
    if len(mutants) != 0:
        status = 'danger'
        msg = 'Could not delete'
        for mut in mutants:
            issues.append(f"Sequence is a parent of mutant {mut.enzyme_name}")

    if status == 'success':
        seq.delete()

    result = {'status': status,
              'msg': msg,
              'issues': issues}

    return jsonify(result=result)

@bp.route('/_merge_seq', methods=['GET', 'POST'])
@roles_required('contributor')
def merge_sequences():
    to_merge = request.form['to_merge']
    merge_with = request.form['merge_with']
    status = 'success'
    msg = 'Sequences merged'
    issues = []

    if to_merge != merge_with:
        seq = Sequence.objects(enzyme_name=to_merge)[0]
        seq_merge = Sequence.objects(enzyme_name=merge_with)[0]
        if seq.enzyme_type == seq_merge.enzyme_type:
            for paper in seq.papers:
                seq_merge.papers.append(paper)

            acts = Activity.objects(enzyme_name=to_merge)
            for act in acts:
                act.enzyme_name = seq_merge.enzyme_name
                act.save()
            seq.delete()
            seq_merge.other_names.append(to_merge)
            seq_merge.save()
            update_seq_papers_status(merge_with)
        else:
            status = 'danger'
            msg = 'Could not merge sequences'
            issues.append('Enzyme types must be the same')
    else:
        status = 'danger'
        msg = 'Could not merge sequences'
        issues.append('Cannot merge with self')

    result = {'status': status,
              'msg': msg,
              'issues': issues}

    return jsonify(result=result)

@bp.route('/_load_sequence_papers', methods=['GET', 'POST'])
@roles_required('contributor')
def load_sequence_papers():
    user = user_datastore.get_user(current_user.id)
    enzyme_name = request.form['name']
    seq = Sequence.objects(enzyme_name=enzyme_name).select_related()[0]

    papers_list = []
    for paper in seq.papers:
        paper_dict = {}
        paper_dict['_id'] = str(paper.id)
        paper_dict['short_citation'] = paper.short_citation
        paper_dict['doi'] = paper.doi
        paper_dict['title'] = paper.title

        if user.has_role('super_contributor') or paper.owner == user:
            paper_dict['can_edit'] = "True"
        else:
            paper_dict['can_edit'] = "False"

        papers_list.append(paper_dict)

    result = {'papers': papers_list}

    return jsonify(result=result)

def update_seq_papers_status(enzyme_name):
    seq = Sequence.objects(enzyme_name=enzyme_name).select_related()[0]
    for paper in seq.papers:
        paper_progress_text, paper_progress = paper_status.paper_metadata_status(paper)
        sequence_progress_text, sequence_progress = paper_status.sequences_status(paper)
        activity_progress_text, activity_progress = paper_status.activity_status(paper)
        status, status_colour = paper_status.get_status(paper_progress, sequence_progress, activity_progress)

        paper.status = status
        paper.save()
