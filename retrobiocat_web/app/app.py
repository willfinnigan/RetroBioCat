from flask import Flask
from redis import Redis
import rq
from retrobiocat_web.config import Config
from flask_talisman import Talisman
from flask_wtf.csrf import CSRFProtect
from flask_jsglue import JSGlue
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_security import Security, MongoEngineUserDatastore, hash_password, current_user
from flask_mail import Mail
from flask_admin import Admin
from flask_session import Session
from flask_mongoengine import MongoEngine
import datetime
from mongoengine import disconnect
from mongoengine import Q

csrf = CSRFProtect()
jsglue = JSGlue()
limiter = Limiter(key_func=get_remote_address, default_limits=["100/minute", "1800/hour"])
talisman = Talisman(content_security_policy=False)
mail = Mail()
admin_ext = Admin()
db = MongoEngine()
session = Session()

from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.app.user_model_forms import ExtendedConfirmRegisterForm, ExtendedRegisterForm

from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Paper, Molecule, Activity
from retrobiocat_web.app.admin import MyAdminIndexView, MyModelView

user_datastore = MongoEngineUserDatastore(db, User, Role)

from retrobiocat_web.app import main_site, retrobiocat, biocatdb

def quick_db(config_class=Config):
    print("Make a quick app to init mongoengine..")
    app = Flask(__name__)
    app.config.from_object(config_class)
    disconnect()
    db.init_app(app)
    print('..done')
    return app

def create_app(config_class=Config, use_talisman=True):
    print("Create app...")
    app = Flask(__name__)
    app.config.from_object(config_class)

    print("Init task queues...")
    app.redis = Redis.from_url(app.config['REDIS_URL'])
    app.task_queue = rq.Queue('tasks', connection=app.redis, default_timeout=600)
    app.network_queue = rq.Queue('network', connection=app.redis, default_timeout=600)
    app.pathway_queue = rq.Queue('pathway', connection=app.redis, default_timeout=600)
    app.retrorules_queue = rq.Queue('retrorules', connection=app.redis, default_timeout=600)
    app.db_queue = rq.Queue('db', connection=app.redis, default_timeout=2400)

    print("Init addons...")
    if use_talisman == True:
        talisman.init_app(app, content_security_policy=False)

    csrf.init_app(app)
    disconnect()
    db.init_app(app)
    jsglue.init_app(app)
    session.init_app(app)
    limiter.init_app(app)
    mail.init_app(app)
    security = Security(app, user_datastore,
                        confirm_register_form=ExtendedConfirmRegisterForm,
                        register_form=ExtendedRegisterForm)

    print("Prepare admin views..")
    admin_ext.init_app(app, index_view=MyAdminIndexView())
    admin_ext.add_view(MyModelView(User))
    admin_ext.add_view(MyModelView(Role))
    admin_ext.add_view(MyModelView(EnzymeType))
    admin_ext.add_view(MyModelView(Sequence))
    admin_ext.add_view(MyModelView(Paper))
    admin_ext.add_view(MyModelView(Molecule))
    admin_ext.add_view(MyModelView(Activity))

    # Create a user to test with
    @app.before_first_request
    def create_user():
        admin = user_datastore.find_or_create_role('admin', description='admin role')
        user_datastore.find_or_create_role('enzyme_types_admin', description='enzyme_types_admin')
        user_datastore.find_or_create_role('contributor', description='contributor')
        user_datastore.find_or_create_role('super_contributor', description='contributor')
        user_datastore.find_or_create_role('rxn_rules_admin', description='rxn_rules_admin')
        user_datastore.find_or_create_role('paper_adder', description='paper_finder')
        user_datastore.find_or_create_role('experimental', description='experimental')
        user_datastore.find_or_create_role('enzyme_champion', description='enzyme_champion')
        user_datastore.find_or_create_role('enzyme_teams', description='enzyme_teams')

        if not user_datastore.get_user(app.config['ADMIN_EMAIL']):
            user = user_datastore.create_user(email=app.config['ADMIN_EMAIL'],
                                              password=hash_password(app.config['ADMIN_PASSWORD']),
                                              first_name='Admin',
                                              last_name='',
                                              affiliation='RetroBioCat',
                                              confirmed_at=datetime.datetime.now())
            user_datastore.add_role_to_user(user, admin)
        print("done")

    @app.context_processor
    def inject_login_mode():
        inject_dict = {}
        inject_dict['login_mode'] = app.config['USE_EMAIL_CONFIRMATION']

        if current_user.is_authenticated:
            user = User.objects(id=current_user.id).select_related()[0]
            if user.has_role('enzyme_teams') and user.enzyme_teams is not None:
                inject_dict['enzyme_teams'] = [enz_type_obj.enzyme_type for enz_type_obj in user.enzyme_teams]
            if user.has_role('enzyme_champion') and user.enzyme_champion is not None:
                inject_dict['enzyme_champion'] = [enz_type_obj.enzyme_type for enz_type_obj in user.enzyme_champion]
            if user.has_role('contributor'):
                inject_dict['user_papers_need_data'] = len(Paper.objects(Q(owner=user) & (Q(status='Data required') | Q(status='Enzymes need protein sequences'))))

            inject_dict['total_team_notifications'] = 0
            inject_dict['team_notifications'] = {}
            inject_dict['champ_notifications'] = {}

            for enz_type in inject_dict['enzyme_teams']:
                num_papers = len(Paper.objects(Q(tags=enz_type) & Q(owner=None) & (Q(status='Data required') | Q(status='Enzymes need protein sequences'))))
                inject_dict['team_notifications'][enz_type] = num_papers
                inject_dict['total_team_notifications'] += num_papers
            for enz_type in inject_dict['enzyme_champion']:
                num_papers = len(Paper.objects(Q(tags=enz_type) & Q(status='Complete - Awaiting review')))
                inject_dict['champ_notifications'][enz_type] = num_papers
                inject_dict['total_team_notifications'] += num_papers

        return inject_dict

    print("Register blueprints...")
    with app.app_context():
        app.register_blueprint(main_site.bp)
        app.register_blueprint(retrobiocat.bp)
        app.register_blueprint(biocatdb.bp)

        return app









