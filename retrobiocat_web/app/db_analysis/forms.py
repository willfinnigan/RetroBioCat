from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.mongo.models.biocatdb_models import SSN_record

class SSN_Form(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    alignment_score = IntegerField('Alignment score cut-off', default=300, validators=[NumberRange(min=1, max=500)])
    include_mutants = BooleanField('Include mutants', default=False)
    only_biocatdb = BooleanField('Only RetroBioCat database sequences', default=False)
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in (list(SSN_record.objects().distinct('enzyme_type')))]
