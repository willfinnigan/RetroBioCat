from flask import render_template, jsonify, session, request, redirect, url_for
from retrobiocat_web.app.biocatdb import bp
import mongoengine as db
from retrobiocat_web.app.biocatdb.functions import sequence_table
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Paper, Sequence, Activity
from retrobiocat_web.app.biocatdb.forms import PapersSearch
from retrobiocat_web.app.biocatdb.functions.papers import papers_table


def filter_papers_by_enzyme_name(papers, enzyme_name):
    seqs = Sequence.objects(enzyme_name=enzyme_name).select_related()
    seq_papers = [paper for paper in [seq.papers for seq in seqs]][0]
    filtered_papers = [paper for paper in papers if paper in seq_papers]
    return filtered_papers

def filter_papers_by_reaction(papers, reaction):
    activities = Activity.objects(reaction=reaction).select_related()
    act_papers = [act.paper for act in activities]
    act_papers = list(set(act_papers))
    filtered_papers = [paper for paper in papers if paper in act_papers]
    return filtered_papers


@bp.route("/papers_search", methods=["GET", "POST"])
def papers_search():
    form = PapersSearch()
    form.set_choices()

    if form.validate_on_submit() == True:
        form_data = form.data
        enzyme_type = form_data['enzyme_type']
        enzyme_name = form_data['enzyme_name']
        reaction = form_data['reaction']
        print(enzyme_type)
        if enzyme_type == 'All':
            enzyme_type = None
        if enzyme_name == 'All':
            enzyme_name = None
        if reaction == 'All':
            reaction = None

        return redirect(url_for("biocatdb.show_papers", enzyme_type=enzyme_type, enzyme_name=enzyme_name, reaction=reaction))

    return render_template('papers_query/papers_query.html', form=form)


# "localhost:5000/sequences?enzyme_type=CAR&paper_id=24934239"
@bp.route("/papers", methods=["GET"])
def show_papers():

    args = request.args.to_dict()

    if 'enzyme_type' in args:
        enzyme_type_query = db.Q(tags=args['enzyme_type'])
    else:
        enzyme_type_query = db.Q()

    papers = Paper.objects(enzyme_type_query)

    if 'enzyme_name' in args:
        papers = filter_papers_by_enzyme_name(papers, args['enzyme_name'])

    if 'reaction' in args:
        papers = filter_papers_by_reaction(papers, args['reaction'])

    paper_ids = [paper.id for paper in papers]
    papers_data = list(Paper.objects(id__in=paper_ids).only(*papers_table.PAPERS_TABLE_FIELDS).order_by('-status').as_pymongo())
    papers_data = papers_table.process_papers_dict(papers_data)

    return render_template('edit_tables/edit_papers.html',
                           papers_data=papers_data, papers_table_height='80vh',
                           papers_button_columns=[],
                           show_owner=True)








if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    papers = Paper.objects()
    enzyme_name = 'NiCAR'

    filtered_papers = filter_papers_by_enzyme_name(papers, enzyme_name)
    print(filtered_papers)

    reaction = 'Carboxylic acid reduction'
    filtered_papers_2 = filter_papers_by_reaction(papers, reaction)
    print(filtered_papers_2)