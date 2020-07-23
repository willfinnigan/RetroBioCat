from retrobiocat_web.app.contributions import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
import yaml
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.init_db import rxn_rules_to_db
from retrobiocat_web.mongo.init_db import biocatdb_excel_to_db, make_molecule_db
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Molecule
import tempfile
from werkzeug.utils import secure_filename
import os
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField
from flask_wtf.file import FileField
import pandas as pd
from retrobiocat_web.app.contributions.functions import paper_status, papers_functions, papers_crossref, create_building_block_db
import time
import datetime
from pathlib import Path
from retrobiocat_web.retro.evaluation.starting_material import StartingMaterialEvaluator

class InitDB(FlaskForm):
    rxns = FileField("Reactions")
    biocatdb = FileField("biocatdb")
    sequences = FileField("sequences")
    submit = SubmitField('Submit')

class Assign(FlaskForm):
    submit = SubmitField('Submit')


@bp.route('/init_db', methods=['GET', 'POST'])
@roles_required('admin')
def init_db():
    form = InitDB()

    if form.validate_on_submit() == True:
        if form.rxns.data != None:
            rxns_data = form.rxns.data
            filename = secure_filename(rxns_data.filename)
            try:
                rxns_data.save(filename)
                yaml_dict = rxn_rules_to_db.load_yaml_dict(filename)
                Reaction.drop_collection()
                current_app.db_queue.enqueue(rxn_rules_to_db.load_into_mongo, yaml_dict)
                flash("Yaml loaded ok, reactions_to_db in queue", 'success')
            except:
                flash("Failed to load reactions", 'fail')
            os.remove(filename)

        if form.biocatdb.data != None:
            filename = secure_filename(form.biocatdb.data.filename)
            form.biocatdb.data.save(filename)
            try:
                df = biocatdb_excel_to_db.load_df(filename)

                Paper.drop_collection()
                Activity.drop_collection()
                Molecule.drop_collection()
                Sequence.drop_collection()
                current_app.db_queue.enqueue(biocatdb_excel_to_db.df_to_db, df)
                current_app.db_queue.enqueue(make_molecule_db.make_fp_db)
                flash("Biocatdb excel loaded ok, added to queue", "success")
            except:
                flash("Problem loading biocatdb excel", "fail")
            os.remove(filename)

        if form.sequences.data != None:
            filename = secure_filename(form.sequences.data.filename)
            form.sequences.data.save(filename)
            try:

                if '.xlsx' in filename:
                    df = pd.read_excel(filename)
                elif '.csv' in filename:
                    df = pd.read_csv(filename)
                else:
                    df = None
                    Exception('File must be .csv or .xlsx')

                current_app.db_queue.enqueue(task_add_sequence_data, df)
                flash("Sequences excel loaded ok, added to queue", "success")
            except Exception as e:
                flash(f"Problem loading sequences excel - {e}", "fail")


            os.remove(filename)

    return render_template('init_db/init_db.html', form=form)

@bp.route('/other_admin_functions', methods=['GET', 'POST'])
@roles_required('admin')
def other_admin_functions():
    return render_template('init_db/other_admin_functions.html')

@bp.route('/_assign_papers', methods=['GET', 'POST'])
@roles_required('admin')
def secret_assign_papers():
    current_app.db_queue.enqueue(task_assign_papers)
    current_app.db_queue.enqueue(task_get_paper_metadata)
    current_app.db_queue.enqueue(biocatdb_init_complete)

    result = {'status': 'success',
              'msg': 'assigning papers',
              'issues': []}

    return jsonify(result=result)

def biocatdb_init_complete():
    print("SUBSTRATE SPECIFICITY DATABASE INITIALISATION COMLPETE")

def task_assign_papers():
    users = User.objects()
    papers = Paper.objects()

    for paper in papers:
        paper_status.update_status(paper)

    for user in users:
        usernames = get_usernames(user)
        for paper in papers:
            if paper.added_by is None or paper.added_by == '':
                activities = Activity.objects(paper=paper)
                for activity in activities:
                    if does_username_match(usernames, activity.added_by_string):
                        activity.added_by = user
                        activity.save()

                        if paper.added_by is None or paper.added_by == '':
                            paper.added_by = user
                            paper.owner = user
                            paper.save()

def get_usernames(user):
    usernames = []
    if len(User.objects(last_name=user.last_name)) == 1 and len(user.last_name) > 2:
        usernames.append(str(user.last_name).lower())
    if len(User.objects(first_name=user.first_name)) == 1 and len(user.first_name) > 2:
        usernames.append(str(user.first_name).lower())

    usernames.append(f"{user.first_name} {user.last_name}".lower())
    usernames.append(f"{user.first_name} {user.last_name}, {user.affiliation}".lower())

    for name in usernames:
        if len(name) > 3:
            usernames.remove(name)

    return usernames

def does_username_match(usernames, added_by_str):
    added_by_str = str(added_by_str).lower()
    if added_by_str != 'nan':
        for name in usernames:
            if name != 'tom':
                if name in added_by_str or name == added_by_str:
                    return True
    return False

def task_add_sequence_data(df):
    users = User.objects()

    for i, row in df.iterrows():
        seq_query = Sequence.objects(enzyme_name=row['Enzyme name'])
        if len(seq_query) != 0:
            seq = seq_query[0]
            seq.sequence = row['Protein Sequence']

            for user in users:
                usernames = get_usernames(user)
                if does_username_match(usernames, row['Added by']):
                    seq.added_by = user
                    seq.owner = user

            seq.save()

def task_get_paper_metadata():

    papers = Paper.objects()
    for paper in papers:
        papers_functions.tag_paper_with_enzyme_types(paper)

        if paper.authors is None or paper.authors == [''] or paper.authors == []:
            title, authors_list, journal, date, cite_mini = papers_crossref.get_metadata_from_crossref(paper.doi)
            if cite_mini == '':
                title, authors_list, journal, date, cite_mini = papers_functions.query_pubmed(paper.doi)
                time.sleep(3)

            if cite_mini != '':
                paper.title = title
                paper.authors = authors_list
                paper.journal = journal
                paper.date = date
                paper.short_citation = cite_mini
                paper.save()

                paper_status.update_status(paper)

def task_init_building_block_db(df):

    df_smi = df[['SMILES']]
    df_smi = df_smi.drop_duplicates()

    db_path = str(Path(__file__).parents[4]) + '/retro/data/buyability/building_blocks.db'

    create_building_block_db.make_db(db_path)
    create_building_block_db.df_to_sql(db_path, df_smi)
    create_building_block_db.create_index(db_path, 'SMILES')


def test_db(db_path):

    try:
        ev = StartingMaterialEvaluator(alternative_db_path=db_path)
        t0 = time.time()
        print(ev.eval('O=CC(=O)Cc1ccccc1'))
        t1 = time.time()
        print(f"Time  = {round(t1 - t0, 4)}")
    except:
        print('DB test failed')
        return False

    return True




if __name__ == '__main__':
    data_folder = str(Path(__file__).parents[4]) + '/retro/data/buyability'
    print(data_folder)