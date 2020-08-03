from retrobiocat_web.app.biocatdb import bp
from flask import render_template, flash, redirect, url_for, session, request, jsonify
from retrobiocat_web.app.biocatdb.model_forms import PaperInfo
from flask_security import roles_required, current_user
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Paper
from retrobiocat_web.app.biocatdb.functions import papers_functions, papers_crossref
from retrobiocat_web.app.app import user_datastore
from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField

class AddPaper(FlaskForm):
    doi = StringField("DOI")
    submit = SubmitField('Go')

@bp.route('/launch_add_paper', methods=['GET', 'POST'])
@roles_required('paper_adder')
def launch_add_paper():
    form = AddPaper()

    if form.validate_on_submit() == True:
        doi = form.doi.data.replace(' ','')
        list_html_to_remove = ['https://doi.org/', 'http://doi.org/', 'http://dx.doi.org/']
        for to_remove in list_html_to_remove:
            if to_remove in doi:
                doi = doi.replace(to_remove, '')
        session['doi'] = doi
        return redirect(url_for("biocatdb.create_paper"))

    else:
        return render_template('add_paper_workflow/add_paper_start.html', form=form)

@bp.route("/create_paper/", methods=["GET", "POST"])
@roles_required('paper_adder')
def create_paper():
    form = PaperInfo()

    if 'doi' in session:
        doi = session['doi']
    else:
        doi = ''

    if form.validate_on_submit() == False:
        user = User.objects(id=current_user.id)[0]
        papers = Paper.objects(doi__iexact=doi)
        if len(papers) != 0:

            paper = papers[0]
            if paper.owner == user or user.has_role('super_contributor'):
                flash("Paper already in the database", 'success')
                return redirect(url_for("biocatdb.submission_main_page", paper_id=paper.id))

            elif paper.owner == None:
                flash('Paper already in database, self-assign the paper to add data', 'warning')
                return redirect(url_for("biocatdb.papers_that_need_data"))

            elif paper.owner != user and not user.has_role('super_contributor'):
                flash("Paper already in the database and is assigned to another user", 'danger')

            else:
                flash("error", 'danger')
                return redirect(url_for("biocatdb.launch_add_paper"))

        else:
            title, authors_list, journal, date, cite_mini = papers_crossref.get_metadata_from_crossref(doi)
            if cite_mini == '':
                title, authors_list, journal, date, cite_mini = papers_functions.query_pubmed(doi)
            form.title.data = title
            form.authors.data = papers_functions.list_to_string(authors_list)
            form.journal.data = journal
            form.date.data = date
            form.short_cit.data = cite_mini
            form.doi.data = doi.replace(' ','').lower()

            can_self_assign = papers_functions.can_self_assign(user)

            if cite_mini == '':
                flash("Paper not found in crossref, pubmed or db, please add manually", 'fail')
            else:
                flash("Paper found, please check information", 'success')


            return render_template('add_paper_workflow/edit_paper_information.html',
                                   form=form, can_self_assign=can_self_assign)

    else:
        user = User.objects(id=current_user.id)[0]
        session.pop('doi')
        if form.self_assign.data == True:
            owner = user
        else:
            owner = None

        new_paper = Paper(doi=form.doi.data.replace(' ','').lower(),
                          short_citation=form.short_cit.data,
                          title=form.title.data,
                          html='https://doi.org/'+form.doi.data,
                          journal=form.journal.data,
                          date=form.date.data,
                          tags=form.tags.data.split(', '),
                          authors=form.authors.data.split(', '),
                          owner=owner,
                          added_by=user,
                          status='Data required')
        new_paper.save()
        flash("Paper saved", 'success')

        if owner == user:
            return redirect(url_for("biocatdb.submission_main_page", paper_id=new_paper.id))

        else:
            return redirect(url_for("biocatdb.launch_add_paper"))


