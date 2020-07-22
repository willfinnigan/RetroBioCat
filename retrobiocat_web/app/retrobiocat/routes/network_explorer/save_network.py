from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request
import json
from flask import current_app
from retrobiocat_web.mongo.models.user_saves import Network
from retrobiocat_web.app.app import user_datastore
from flask_security import current_user
import datetime

def add_save(save_links, unique_id, name):
    saved = False
    to_save = [name, str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")), unique_id]

    for i, save_tuple in enumerate(save_links):
        if save_tuple[2] == unique_id:
            save_links[i] = to_save
            saved = True
    if saved == False:
        save_links.append(to_save)
    return save_links

@bp.route('/_save_network', methods=['GET', 'POST'])
def save_network_func():
    if not current_user.is_authenticated:
        return jsonify(result={'error' : 'Please log in to save'})

    task_id = request.form['task_id']
    unique_id = request.form['unique_id']
    name = request.form['name']
    public = bool(request.form['public'])
    description = request.form['description']

    data = json.loads(current_app.redis.get(task_id))
    data['save_id'] = unique_id
    data['save_links'] = add_save(data['save_links'], unique_id, name)

    current_app.redis.mset({task_id: json.dumps(data)})

    user = user_datastore.get_user(current_user.id)

    if len(Network.objects(uuid=task_id)) > 0:
        old_save = Network.objects(uuid=task_id)[0]
        if old_save.owner != user:
            return jsonify(result={'error' : "Can not overwrite another user's save"})

    network = Network(uuid=unique_id,
                      name=name,
                      description=description,
                      public=public,
                      owner=user,
                      target_smiles=data["target_smiles"],
                      data=data,
                      time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    network.save()

    result = {'save_links' : data['save_links'][::-1]}
    return jsonify(result=result)