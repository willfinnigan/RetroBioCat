import mongoengine as db
from retrobiocat_web.mongo.models.user_models import User
from datetime import datetime
import uuid
from bson import ObjectId
from retrobiocat_web.mongo.models.comments import Comment

class Reaction(db.Document):
    name = db.StringField(unique=True)
    smarts = db.ListField(db.StringField())
    enzyme_types = db.ListField(db.StringField())
    cofactors = db.DictField()
    positive_tests = db.ListField(db.StringField())
    negative_tests = db.ListField(db.StringField())
    type = db.StringField()
    experimental = db.BooleanField(default=False)
    two_step = db.BooleanField(default=False)
    requires_absence_of_water = db.BooleanField(default=False)

class Issue(db.Document):
    reaction = db.ReferenceField(Reaction)
    issue_reaction_smiles = db.StringField()
    issue_reaction_svg = db.StringField()
    raised_by = db.ReferenceField(User, reverse_delete_rule=2)
    status = db.StringField(default='Open')
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    public = db.BooleanField(default=False)
    date = db.DateTimeField(default=datetime.utcnow)

class ReactionSuggestion(db.Document):
    name = db.StringField()
    smarts = db.ListField(db.StringField())
    details = db.StringField()
    owner = db.ReferenceField(User)
    status = db.StringField(default='Open')
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    date = db.DateTimeField(default=datetime.utcnow)
