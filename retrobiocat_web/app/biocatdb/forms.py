from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, SubmitField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired, NumberRange, ValidationError
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.enzyme_identification import query_mongodb
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence, Activity, SeqSimNet
from retrobiocat_web.mongo.models.reaction_models import Reaction

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
    enzymes = SelectField('Enzyme type', validators=[DataRequired(), is_enzyme])
    reactions = SelectField('Reaction', validators=[is_reaction])
    data_level = SelectField('Data level', choices=specificity_data_choices)
    num_choices = IntegerField('Max enzymes per substrate', default=4, validators=[NumberRange(min=1, max=100)])
    max_hits = IntegerField('Max similar substrates', default=10, validators=[NumberRange(min=1, max=100)])
    target_smiles = StringField('Product SMILES', validators=[is_accepted_by_rdkit])
    similarity = DecimalField('Similarity cutoff', default=0.6, validators=[NumberRange(min=0.1, max=1)])
    auto_data = BooleanField('Include automatically generated data', default=False)
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzymes.choices = [(c, c) for c in ['All'] + (list(Activity.objects().distinct('enzyme_type')))]
        self.reactions.choices = [(c, c) for c in ['All'] + (list(Activity.objects().distinct('reaction')))]

    def validate(self):
        if not FlaskForm.validate(self):
            return False
        result = True
        seen = set()

        if self.enzymes.data == 'All' and self.reactions.data == 'All':
            self.enzymes.errors.append("Can't both be All")
            self.reactions.errors.append("Can't both be All")
            result = False
        else:
            seen.add(self.enzymes.data)
            seen.add(self.reactions.data)
        return result

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
        self.enzyme_type.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_type')))]

class PapersSearch(FlaskForm):
    enzyme_type = SelectField('Enzyme type')
    enzyme_name = SelectField('Enzyme name')
    reaction = SelectField('Reaction')
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_type')))]
        self.enzyme_name.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_name')))]
        self.reaction.choices = [(c, c) for c in ['All'] + (list(Reaction.objects().distinct('name')))]

class SubstrateScopeForm(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    enzyme_name = SelectField('Enzyme name')
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_type')))]
        self.enzyme_name.choices = [(c, c) for c in ['All'] + (list(Sequence.objects().distinct('enzyme_name')))]

    def validate(self):
        if not FlaskForm.validate(self):
            return False
        result = True
        seen = set()

        if self.enzyme_type.data == 'All' and self.enzyme_name.data == 'All':
            self.enzyme_type.errors.append("Can't both be All")
            self.enzyme_name.errors.append("Can't both be All")
            result = False
        else:
            seen.add(self.enzyme_type.data)
            seen.add(self.enzyme_name.data)
        return result

class SSN_Form(FlaskForm):
    enzyme_type = SelectField('Enzyme name')
    alignment_score = IntegerField('Alignment score cut-off', default=300, validators=[NumberRange(min=1, max=500)])
    combine_mutants = BooleanField('Combine mutants with parents', default=True)
    only_biocatdb = BooleanField('Only BioCatDB Sequences', default=False)
    max_nodes = IntegerField('Max nodes', default=1000, validators=[NumberRange(min=10, max=2000)])
    edges_to_cluster_on = IntegerField('Edges for clusters', default=10, validators=[NumberRange(min=5, max=30)])
    submit = SubmitField('Submit')

    def set_choices(self):
        self.enzyme_type.choices = [(c, c) for c in (list(SeqSimNet.objects().distinct('enzyme_type')))]

