from retrobiocat_web.mongo.models.user_models import User


PAPERS_TABLE_FIELDS = ['id', 'short_citation', 'title', 'doi', 'status', 'tags', 'owner']

def process_papers_dict(papers_data, show_owner=True):
    owners_dict = {}
    for i, data in enumerate(papers_data):
        papers_data[i]['_id'] = str(papers_data[i]['_id'])

        if show_owner == False:
            if 'owner' in papers_data[i]:
                papers_data[i].pop('owner')

        else:
            if 'owner' in papers_data[i]:
                owner_id = str(papers_data[i]['owner'])
                if owner_id not in owners_dict:
                    owner = User.objects(id=papers_data[i]['owner'])[0]
                    owners_dict[owner_id] = f"{owner.first_name} {owner.last_name}"
                papers_data[i]['owner'] = owners_dict[owner_id]
            else:
                papers_data[i]['owner'] = '-'

    return papers_data