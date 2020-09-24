from retrobiocat_web.mongo.models.biocatdb_models import Paper, EnzymeType
import mongoengine as db


def get_enzyme_paper_progress(enzyme_type_obj):
    enz_type = enzyme_type_obj.enzyme_type
    num_papers = len(Paper.objects(tags=enz_type))
    num_complete_papers = len(Paper.objects(db.Q(tags=enz_type) & (db.Q(status='Complete') | db.Q(status='Complete - Awaiting review'))))

    if num_papers == 0:
        progress = ["0%", f"{num_complete_papers} out of {num_papers}", 'bg-danger', enzyme_type_obj.full_name]
    elif num_complete_papers == 0:
        progress = ["0%", f"{num_complete_papers} out of {num_papers}", 'bg-danger', enzyme_type_obj.full_name]
    else:
        pc_complete = round((num_complete_papers / num_papers) * 100, 1)

        if pc_complete > 80:
            colour = 'bg-success'
        elif pc_complete > 40:
            colour = 'bg-warning'
        else:
            colour = 'bg-danger'
        progress = [f"{pc_complete}%", f"{num_complete_papers} out of {num_papers}", colour, enzyme_type_obj.full_name]

    return progress
