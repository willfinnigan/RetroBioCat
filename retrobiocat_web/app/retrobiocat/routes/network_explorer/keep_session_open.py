from retrobiocat_web.app.retrobiocat import bp, forms
from flask import render_template, jsonify, session, request, make_response
from flask import current_app


@bp.route('/_keep_session_open', methods=['GET', 'POST'])
def keep_session_open():
    task_id = request.form['task_id']
    current_app.redis.expire(task_id, 5*60)
    return jsonify({})