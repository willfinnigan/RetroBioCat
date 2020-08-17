from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Sequence

def get_enzyme_data(query):
    seq_fields = ['id', 'enzyme_type', 'enzyme_name', 'sequence', 'sequence_unavailable', 'accession', 'pdb', 'mutant_of', 'notes', 'papers', 'owner', 'other_names']
    enzyme_data = list(Sequence.objects(query).only(*seq_fields).order_by('enzyme_type').as_pymongo())
    owners_dict = {}
    for i, data in enumerate(enzyme_data):
        enzyme_data[i]['_id'] = str(enzyme_data[i]['_id'])
        enzyme_data[i]['sequence_unavailable'] = str(enzyme_data[i]['sequence_unavailable']).replace('False', '')

        if 'papers' not in enzyme_data[i]:
            enzyme_data[i]['papers'] = 0
        else:
            enzyme_data[i]['papers'] = len(enzyme_data[i]['papers'])

        if 'owner' in enzyme_data[i]:
            owner_id = str(enzyme_data[i]['owner'])
            if owner_id not in owners_dict:
                owner = User.objects(id=enzyme_data[i]['owner'])[0]
                owners_dict[owner_id] = f"{owner.first_name} {owner.last_name}"
            enzyme_data[i]['owner'] = owners_dict[owner_id]
        else:
            enzyme_data[i]['owner'] = ''

        if 'sequence' in enzyme_data[i]:
            if len(enzyme_data[i]['sequence']) > 50:
                enzyme_data[i]['sequence'] = enzyme_data[i]['sequence'][0:50] + "..."

    return enzyme_data