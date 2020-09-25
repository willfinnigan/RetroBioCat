from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, Paper


def paper_metadata_status(paper):
    progress = 100
    required = []
    if paper.short_citation is None or paper.short_citation == '':
        progress -= 10
        required.append('Short citation')
    if paper.title is None or paper.title == '':
        progress -= 10
        required.append('Title')
    if paper.authors is None or paper.authors == [] or paper.authors == ['']:
        progress -= 10
        required.append('Authors')
    if paper.journal is None or paper.journal == '':
        progress -= 10
        required.append('Journal')
    if paper.date is None:
        progress -= 10
        required.append('Date')
    if paper.tags is None or paper.tags == [] or paper.tags == ['']:
        progress -= 10
        required.append('Tags')

    if progress == 100:
        progress_text = 'Complete'
    else:
        progress_text = f'Required: {required}'

    return progress_text, str(progress)+"%"

def sequences_status(paper):
    progress = 0
    progress_text = 'Enzymes need to be added to paper'

    seqs = Sequence.objects(papers=paper)

    if len(seqs) == 0:
        return progress_text, str(progress)+'%'

    else:
        progress = 50
        progress_text = 'Protein sequences need to be added (or marked as unavailable)'

        for seq in seqs:
            if (seq.sequence is None or seq.sequence == '') and (seq.sequence_unavailable is False):
                return progress_text, str(progress)+'%'

        progress = 100
        progress_text = 'Complete'
        return progress_text, str(progress)+'%'

def activity_status(paper):
    progress = 5
    progress_text = "Please add activity data"

    activity = Activity.objects(paper=paper)

    if len(activity) != 0:
        progress = 90
        progress_text = "Activity data added, please ensure all activity data (including negative data) is included.  Awaiting review."

        if paper.reviewed == True:
            progress = 100
            progress_text = 'Complete'

        return progress_text, str(progress) + '%'

    else:
        return progress_text, str(progress)+'%'

def get_status(paper_progress, sequence_progress, activity_progress, paper):
    paper_progress = int(paper_progress[:-1])
    sequence_progress = int(sequence_progress[:-1])
    activity_progress = int(activity_progress[:-1])

    if paper_progress == 100 and sequence_progress == 100 and activity_progress == 100:
        status, colour = 'Complete', 'green'
    elif paper_progress == 100 and sequence_progress == 100 and activity_progress == 90:
        status, colour = 'Complete - Awaiting review', 'green'
    elif paper_progress != 100 and sequence_progress == 100 and activity_progress >= 90:
        status, colour = 'Requires paper metadata', 'green'
    else:
        if sequence_progress == 50:
            status, colour = 'Enzymes need protein sequences', 'orange'
        else:
            status, colour = 'Data required', 'red'

    if paper.has_issues == True:
        status, colour = 'Issues need to be resolved', 'red'

    return status, colour

def update_status(paper):
    paper_progress_text, paper_progress = paper_metadata_status(paper)
    sequence_progress_text, sequence_progress = sequences_status(paper)
    activity_progress_text, activity_progress = activity_status(paper)
    status, status_colour = get_status(paper_progress, sequence_progress, activity_progress, paper)

    paper.status = status
    paper.save()
