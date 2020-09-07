import mongoengine as db
from retrobiocat_web.mongo.models.user_models import User
from datetime import datetime
import uuid
from bson import ObjectId

class Comment(db.Document):
    owner = db.ReferenceField(User)
    text = db.StringField()
    date = db.DateTimeField(default=datetime.utcnow)
