from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, TextAreaField, SubmitField, BooleanField, DateField
from wtforms.validators import DataRequired, NumberRange, ValidationError, Length
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Paper

def is_doi_taken(form, field):
    for obj in Paper.objects():
        if field.data == obj.doi:
            raise ValidationError(f'{field.data} is already in the database')

def is_type_taken(form, field):
    for obj in EnzymeType.objects():
        if field.data == obj.enzyme_type:
            raise ValidationError(f'{field.data} is already an enzyme type in the database')


class PaperInfo(FlaskForm):
    short_cit = StringField(validators=[DataRequired()])
    doi = StringField(validators=[DataRequired(), is_doi_taken])
    journal = StringField()
    date = DateField(validators=[DataRequired()])
    title = StringField()
    authors = StringField()
    self_assign = BooleanField(default=False)
    tags = StringField()
    submit = SubmitField('Save')

class EnzymeTypeForm(FlaskForm):
    enzyme_type = StringField(validators=[DataRequired(), Length(max=120), is_type_taken])
    full_name = StringField()
    other_abbreviations = StringField()
    description = TextAreaField(validators=[Length(max=255)])
    submit = SubmitField()