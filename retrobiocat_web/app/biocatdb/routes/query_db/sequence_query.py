from flask import render_template, jsonify, session, request
from retrobiocat_web.app.biocatdb import bp



# "localhost:5000/sequences?enzyme_type=CAR&enzyme_name=mpCAR"

@bp.route("/sequences", methods=["GET"])
def show_sequences():

    args = request.args.to_dict()

    if 'enzyme_type' in args:
        pass
    if 'enzyme_name' in args:
        pass

    return 'success', 200
