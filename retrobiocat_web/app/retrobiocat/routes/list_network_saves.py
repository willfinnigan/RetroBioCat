from retrobiocat_web.app.retrobiocat import bp
from flask import render_template, jsonify, request
from flask import redirect
import uuid
from retrobiocat_web.mongo.models.user_saves import Network
from retrobiocat_web.app.app import user_datastore
from flask_security import current_user
from flask import url_for
from rdkit import Chem

@bp.route('/my_network_saves')
def my_network_saves():
    if not current_user.is_authenticated:
        return redirect('/')

    user = user_datastore.get_user(current_user.id)

    saves = []
    for network in Network.objects(owner=user):
        try:
            img = str(Chem.MolFromSmiles(network.target_smiles))
            img = img[:-1] + 'width="100" height="100" >'
        except:
            img = ''

        network_dict = {'Name' : network.name,
                        'Target' : img,
                        'Decription' : network.description,
                        'Saved at' : str(network.time),
                        'Sharable' : str(network.public),
                        'UUID' : network.uuid,
                        'Link': url_for('retrobiocat.network_explorer', task_id=network.uuid, _external=True)}
        saves.append(network_dict)

    headings = ['Name', 'Target', 'Decription', 'Saved at', 'Sharable']


    return render_template('my_network_saves/my_network_saves.html',
                           headings=headings, row_saves=saves)


@bp.route('/delete_network_save', methods=['GET', 'POST'])
def delete_network_save():
    save_id = request.form['save_id']

    user = user_datastore.get_user(current_user.id)
    query = Network.objects(uuid=save_id)
    if len(query) > 0:
        to_delete = query[0]
        if to_delete.owner == user:
            to_delete.delete()

    return jsonify({})



