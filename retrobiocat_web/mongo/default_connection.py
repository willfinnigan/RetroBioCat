import mongoengine as db
import os

def make_default_connection():
    MONGODB_HOST = os.environ.get('MONGO_HOST') or 'localhost'
    MONGODB_DB = os.environ.get('MONGODB_DB') or 'mydatabase'
    MONGODB_PORT = os.environ.get('MONGODB_PORT') or 27017
    MONGO_USERNAME = os.environ.get('MONGO_USERNAME') or ''
    MONGO_PASSWORD = os.environ.get('MONGO_PASSWORD') or ''

    db.connect(MONGODB_DB,
               host=MONGODB_HOST,
               port=MONGODB_PORT,
               username=MONGO_USERNAME,
               password=MONGO_PASSWORD)
