from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.mongo.models.biocatdb_models import SSN_record, EnzymeType

class SSN_Form(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    alignment_score = IntegerField('Score', default=100, validators=[NumberRange(min=40, max=300)])
    hide_mutants = BooleanField('Hide mutants', default=False)
    only_biocatdb = BooleanField('Only RetroBioCat database sequences', default=False)
    submit = SubmitField('Submit')

    def set_choices(self):
        ssn_records = SSN_record.objects().distinct('enzyme_type')
        enzyme_types = EnzymeType.objects()

        list_enzyme_types = []
        for enz_type in enzyme_types:
            if enz_type in ssn_records:
                list_enzyme_types.append(enz_type.enzyme_type)

        list_enzyme_types = sorted(list_enzyme_types)
        self.enzyme_type.choices = [(c, c) for c in list_enzyme_types]


if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    ssn_records = SSN_record.objects().distinct('enzyme_type')
    print(ssn_records)
    enzyme_types = EnzymeType.objects()

    list_enzyme_types = []
    for enz_type in enzyme_types:
        if enz_type in ssn_records:
            list_enzyme_types.append(enz_type.enzyme_type)

    list_enzyme_types = sorted(list_enzyme_types)
    print(list_enzyme_types)