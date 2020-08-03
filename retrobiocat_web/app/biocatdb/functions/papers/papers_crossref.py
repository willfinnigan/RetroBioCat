from crossref.restful import Works
import datetime


def query_crossref(doi):
    works = Works()
    result = works.doi(doi)
    return result

def date_to_datetime(date):
    if len(date) == 3:
        dt = datetime.date(year=date[0], month=date[1], day=date[2])
    elif len(date) == 2:
        dt = datetime.date(year=date[0], month=date[1], day=1)
    elif len(date) == 1:
        dt = datetime.date(year=date[0], month=1, day=1)
    else:
        dt = None

    return dt

def get_authors_list(authors):
    authors_list = []
    for author_dict in authors:
        name = f"{author_dict['given']} {author_dict['family']}"
        authors_list.append(name)
    return authors_list

def select_a_date(result):
    if 'created' in result:
        date = result['created']['date-parts'][0]
    elif 'deposited' in result:
        date = result['deposited']['date-parts'][0]
    elif 'published-online' in result:
        date = result['published-online']['date-parts'][0]
    elif 'published-print' in result:
        date = result['published-print']['date-parts'][0]
    elif 'issued' in result:
        date = result['issued']['date-parts'][0]
    else:
        date = ''
    return date


def get_relevant_info(result):
    try:
        journal = result['container-title'][0]
        journal_short = result['short-container-title'][0]
        authors = result['author']
        authors_list = get_authors_list(authors)
        title = result['title'][0]
        date = select_a_date(result)
        date = date_to_datetime(date)

        if date != None:
            cite_mini = f"{authors[0]['family']} et al, {date.year}, {journal_short}"
        else:
            cite_mini = f"{authors[0]['family']} et al, {journal_short}"

    except Exception as e:
        print(e)
        return '', '', '', '', ''

    return title, authors_list, journal, date, cite_mini


def get_metadata_from_crossref(doi):
    result = query_crossref(doi)
    if result == None:
        return '', '', '', '', ''
    else:
        return get_relevant_info(result)

def print_result(result):
    for key in result:
        print()
        print(f"Key: {key}  Result: {result[key]}")

if __name__ == '__main__':

    #doi = '10.1016/j.molcatb.2013.10.023'
    #doi = '10.1021/acscatal.8b02386'
    doi = "10.1002/cctc.201601249"
    #doi = '10.1039/c7ob02299a'
    #doi = '10.1016/S0014-5793(00)01992-X'  # issue
    #doi = '10.1002/cctc.201300842'  # issue
    #doi = '10.1023/A:1005479910765'  # issue

    print(get_metadata_from_crossref(doi))

    print_result(query_crossref(doi))


