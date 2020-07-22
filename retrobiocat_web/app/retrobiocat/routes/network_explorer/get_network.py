from retrobiocat_web.app.retrobiocat import bp
from flask import render_template
import json
from flask import current_app, redirect
import uuid
from retrobiocat_web.mongo.models.user_saves import Network
from retrobiocat_web.app.app import user_datastore
from flask_security import current_user


# This function is not used and can be deleted
def get_users_saves(all_saves):
    user = user_datastore.get_user(current_user.id)
    user_saves = []
    for save_tuple in all_saves:
        if len(Network.objects(uuid=save_tuple[2])) > 0:
            mongo_network = Network.objects(uuid=save_tuple[2])[0]
            if mongo_network.owner == user:
                user_saves.append(save_tuple)
    return user_saves

def load_from_mongo(id):
    data = False
    if len(Network.objects(uuid=id)) > 0:

        if current_user.is_authenticated:
            user = user_datastore.get_user(current_user.id)
        else:
            user = False

        mongo_network = Network.objects(uuid=id)[0]
        if (mongo_network.public == True) or (mongo_network.owner == user):
            data = mongo_network.data
            if mongo_network.owner != user:
                data['save_id'] = str(uuid.uuid4())
                data['save_links'] = []
    return data

def load_from_redis(id):
    return json.loads(current_app.redis.get(id))

def save_new_redis(data):
    new_task_id = str(uuid.uuid4())
    current_app.redis.mset({new_task_id: json.dumps(data)})
    current_app.redis.expire(new_task_id, 5*60)
    return new_task_id

@bp.route("/network_explorer/<task_id>/", methods=["GET"])
def network_explorer(task_id):
    data = load_from_mongo(task_id)
    if data != False:
        new_task_id = save_new_redis(data)
        return redirect('/network_explorer/' + new_task_id + '/')
    else:
        # Should put an elif and have error on else
        result = load_from_redis(task_id)
        return render_template('network_explorer/network_explorer.html',
                               save_id=result['save_id'],
                               save_name=result['save_name'],
                               save_links=result['save_links'],
                               nodes=result['nodes'],
                               edges=result['edges'],
                               options=result['options'],
                               max_reactions=result['max_reactions'],
                               retrorules_diameters=[],
                               task_id=task_id)



