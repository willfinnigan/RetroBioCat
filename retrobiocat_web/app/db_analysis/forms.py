from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.mongo.models.biocatdb_models import SSN_record

class SSN_Form(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    alignment_score = IntegerField('Score', default=100, validators=[NumberRange(min=40, max=300)])
    hide_mutants = BooleanField('Hide mutants', default=False)
    only_biocatdb = BooleanField('Only RetroBioCat database sequences', default=False)
    submit = SubmitField('Submit')

    def set_choices(self):
        ssn_records = SSN_record.objects().select_related()
        enzyme_types = sorted([ssn.enzyme_type.enzyme_type for ssn in ssn_records])
        self.enzyme_type.choices = [(c, c) for c in enzyme_types]
