from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.mongo.models.biocatdb_models import SSN_record, EnzymeType

class SSN_Form(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    alignment_score = IntegerField('Score', validators=[NumberRange(min=40, max=300)])
    hide_mutants = BooleanField('Hide mutants', default=False)
    only_biocatdb = BooleanField('Only RetroBioCat database sequences', default=False)
    submit = SubmitField('Submit')

    def set_choices(self):
        ssn_records = SSN_record.objects().distinct('enzyme_type')
        enzyme_types = EnzymeType.objects()

        list_enzyme_types = []
        enzyme_descriptions = {}
        for enz_type in enzyme_types:
            if enz_type in ssn_records:
                enzyme_descriptions[enz_type.enzyme_type] = f"{enz_type.enzyme_type} - {enz_type.full_name}"
                list_enzyme_types.append(enz_type.enzyme_type)

        list_enzyme_types = sorted(list_enzyme_types)

        self.enzyme_type.choices = []
        for key in list_enzyme_types:
            self.enzyme_type.choices.append((key, enzyme_descriptions[key]))
