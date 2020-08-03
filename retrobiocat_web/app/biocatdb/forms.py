from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.enzyme_identification import query_mongodb
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType

def is_accepted_by_rdkit(form, field):
    if node_analysis.rdkit_smile(field.data) == None:
        if field.data != '':
            raise ValidationError('SMILES not accepted by rdkit')

def is_reaction(form, field):
    reaction_names = list(field.data.split(", "))
    for reaction in reaction_names:
        if reaction not in (query_mongodb.get_reactions_in_db() + ['All']):
            if reaction != '':
                raise ValidationError('Reaction not defined in main_site')

def is_enzyme(form, field):
    enzyme_names = list(field.data.split(", "))
    for enzyme in enzyme_names:
        if enzyme not in (query_mongodb.get_enzymes_in_db() + ['All']):
            if enzyme != ['']:
                raise ValidationError('Enzyme not defined in main_site')

substrate_node_choices = [('Starting material', 'Starting material'),
                          ('Relative complexity', 'Relative complexity'),
                          ('Off', 'Off')]

enzyme_node_choices = [('Substrate specificity', 'Substrate specificity'),
                       ('Complexity change', 'Complexity change'),
                        ('Off', 'Off')]

edges_choices = [('Off', 'Off'), ('Complexity change', 'Complexity change')]

specificity_data_choices = [('All', 'All'),
                            ('Categorical', 'Categorical'),
                            ('Quantitative', 'Quantitative'),
                            ('Specific Activity', 'Specific Activity'),
                            ('Conversion', 'Conversion')]

class SubstrateForm(FlaskForm):
    enzymes = StringField('Enzymes', validators=[DataRequired(), is_enzyme])
    reactions = StringField('Reactions', validators=[is_reaction])
    data_level = SelectField('Data level', choices=specificity_data_choices)
    num_choices = IntegerField('Max enzymes per substrate', default=4, validators=[NumberRange(min=1, max=100)])
    max_hits = IntegerField('Max similar substrates', default=10, validators=[NumberRange(min=1, max=100)])
    target_smiles = StringField('Product SMILES', validators=[is_accepted_by_rdkit])
    similarity = DecimalField('Similarity cutoff', default=0.6, validators=[NumberRange(min=0.1, max=1)])
    auto_data = BooleanField('Include automatically generated data', default=False)
    submit = SubmitField('Submit')

class Network_Vis_Options(FlaskForm):
    colour_substrates = SelectField('Colour substrate nodes', choices=substrate_node_choices)
    colour_reactions = SelectField('Colour reaction nodes', choices=enzyme_node_choices)

class SequenceByName(FlaskForm):
    name = StringField('Enzyme name', validators=[DataRequired()])
    submit = SubmitField('Submit')

class SequenceByType(FlaskForm):
    def __init__(self, choices):
        super().__init__()
        self.type = SelectField('Enzyme type', choices=choices)
        self.submit = SubmitField('Submit')

class SequenceBySequence(FlaskForm):
    sequence = StringField('Protein sequence', validators=[DataRequired()])
    submit = SubmitField('Submit')

class SequenceSearch(FlaskForm):
    enzyme_type = SelectField('Enzyme type')
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in ['All'] + (list(EnzymeType.objects().distinct('enzyme_type')))]
        self.enzyme_type.choices = [('All', 'All'), ('AAD', 'AAD'), ('AADH', 'AADH'), ('AAO', 'AAO'), ('ADC', 'ADC'), ('ADH', 'ADH'), ('ADH(FAD)', 'ADH(FAD)'), ('AHR', 'AHR'), ('AKR', 'AKR'), ('ALR', 'ALR'), ('Adenylating amidase', 'Adenylating amidase'), ('AlDH', 'AlDH'), ('AlOx', 'AlOx'), ('AlaDH', 'AlaDH'), ('Aldolase', 'Aldolase'), ('AmDH', 'AmDH'), ('AmOx', 'AmOx'), ('Amidase', 'Amidase'), ('BBE', 'BBE'), ('BVMO', 'BVMO'), ('CAR', 'CAR'), ('CMT', 'CMT'), ('Chemical', 'Chemical'), ('DC', 'DC'), ('DERA Aldolase', 'DERA Aldolase'), ('EDDS lyase', 'EDDS lyase'), ('EH', 'EH'), ('ERED', 'ERED'), ('Esterase', 'Esterase'), ('FMO', 'FMO'), ('Hydratase', 'Hydratase'), ('Hydroxynitrile Lyase', 'Hydroxynitrile Lyase'), ('IRED', 'IRED'), ('KRED', 'KRED'), ('Kinase', 'Kinase'), ('Limonene Hydratase', 'Limonene Hydratase'), ('Lipase', 'Lipase'), ('NHase', 'NHase'), ('NIR', 'NIR'), ('NMT', 'NMT'), ('NTR', 'NTR'), ('Nitrilase', 'Nitrilase'), ('Oxd', 'Oxd'), ('P450', 'P450'), ('PAL', 'PAL'), ('PAM', 'PAM'), ('PNP', 'PNP'), ('PSase', 'PSase'), ('Penicillin Acylase', 'Penicillin Acylase'), ('Phosphopentomutase', 'Phosphopentomutase'), ('SMO', 'SMO'), ('SOI', 'SOI'), ('TA', 'TA'), ('TAL', 'TAL'), ('TAM', 'TAM'), ('TDL', 'TDL'), ('TPL', 'TPL'), ('Threonine Aldolase', 'Threonine Aldolase'), ('TrpS', 'TrpS'), ('XOR', 'XOR')]

        print(self.enzyme_type.choices)



