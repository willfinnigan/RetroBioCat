import mongoengine as db

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
