import mongoengine as db
from flask_security import UserMixin, RoleMixin

class Role(db.Document, RoleMixin):
    name = db.StringField(max_length=80, unique=True)
    description = db.StringField(max_length=255)

    def __unicode__(self):
        return self.name

    def __str__(self):
        return self.name

class User(db.Document, UserMixin):
    email = db.StringField(max_length=255, unique=True)
    password = db.StringField(max_length=255)
    active = db.BooleanField(default=True)
    fs_uniquifier = db.StringField(max_length=255)
    confirmed_at = db.DateTimeField()
    roles = db.ListField(db.ReferenceField(Role), default=[])
    first_name = db.StringField(max_length=255)
    last_name = db.StringField(max_length=255)
    affiliation = db.StringField(max_length=255)
    email_opt_in = db.BooleanField()
    current_login_at = db.DateTimeField()
    current_login_ip = db.StringField()
    last_login_at = db.DateTimeField()
    last_login_ip = db.StringField()
    login_count = db.IntField()

    def __unicode__(self):
        return self.email

    def __str__(self):
        return self.email