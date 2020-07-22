from pymed import PubMed
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Sequence, EnzymeType
import datetime
import mongoengine as db


def query_pubmed(doi):
    pubmed = PubMed(tool="RetroBioCat", email="william.finnigan@manchester.ac.uk")
    query = f'{doi}[Location ID]'
    results = list(pubmed.query(query, max_results=1))

    title = ''
    authors_list = []
    date = ''
    journal = ''
    cite_mini = ''

    if (len(results) != 0):
        article = results[0]
        if article.doi == doi:
            title = article.title
            journal = article.journal
            authors_list = []
            for author_dict in article.authors:
                authors_list.append(f"{author_dict['firstname']} {author_dict['lastname']}")
            date = article.publication_date
            cite_mini = f"{article.authors[0]['lastname']} et al, {date.year}"
    return title, authors_list, journal, date, cite_mini

def data_from_db_entry(paper):
    title = paper.title
    authors_list = paper.authors
    journal = paper.journal
    date = paper.date
    cite_mini = paper.short_citation
    tags = list_to_string(paper.tags)
    return title, authors_list, journal, date, cite_mini, tags

def list_to_string(list_field):
    list_string = ""
    for list_item in list_field:
        if list_string != "":
            list_string += ", "
        list_string += list_item
    return list_string


def tag_paper_with_enzyme_types(paper):
    all_types = EnzymeType.objects()
    for enz_type in all_types:
        if enz_type in paper.tags:
            paper.tags.remove(enz_type)

    seqs = Sequence.objects(papers=paper)
    for seq in seqs:
        if (seq.enzyme_type not in paper.tags) and (seq.enzyme_type.lower() != 'chemical'):
            paper.tags.append(seq.enzyme_type)

    paper.save()

def can_self_assign(user):
    q_user = db.Q(owner=user)
    q_no_data = db.Q(status__nin=['Complete - Awaiting review', 'Complete'])
    q_has_data = db.Q(status__in=['Complete - Awaiting review', 'Complete'])
    num_papers_need_data = len(Paper.objects(q_user & q_no_data))
    num_papers_with_data = len(Paper.objects(q_user & q_has_data))

    if num_papers_need_data > num_papers_with_data:
        return False
    else:
        return True


if __name__ == '__main__':
    doi = "10.1002/cctc.201601249"
    title, authors_list, journal, date, cite_mini = query_pubmed(doi)
    print(type(date))
    print(authors_list)