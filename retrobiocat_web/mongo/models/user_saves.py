import mongoengine as db
from retrobiocat_web.mongo.models.user_models import User

class Network(db.Document):
    uuid = db.StringField(primary_key=True)
    target_smiles = db.StringField()
    name = db.StringField(max_length=120, required=True)
    description = db.StringField(max_length=255)
    public = db.BooleanField()
    owner = db.ReferenceField(User, reverse_delete_rule=db.CASCADE)
    data = db.DictField()
    time = db.DateTimeField()

class MyMolecule(db.Document):
    owner = db.ReferenceField(User, reverse_delete_rule=db.CASCADE)
    smiles = db.StringField()
    name = db.StringField()
    svg = db.StringField()
