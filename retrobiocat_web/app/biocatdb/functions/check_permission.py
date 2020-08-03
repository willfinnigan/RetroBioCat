from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.app.app import user_datastore

def check_paper_permission(user_id, paper):
    user = User.objects(id=user_id).select_related()[0]

    if user.has_role('super_contributor'):
        return True
    if paper.owner == user:
        return True
    if user.has_role('enzyme_champion'):
        champ_types = [e.enzyme_type for e in user.enzyme_champion]
        if any(i in champ_types for i in paper.tags):
            return True
    return False

def check_seq_permissions(user_id, seq):
    user = User.objects(id=user_id).select_related()[0]

    if user.has_role('super_contributor'):
        return True
    if seq.owner == user:
        return True
    if user.has_role('enzyme_champion'):
        champ_types = [e.enzyme_type for e in user.enzyme_champion]
        if seq.enzyme_type in champ_types:
            return True
    return False

