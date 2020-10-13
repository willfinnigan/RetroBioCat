import mongoengine as db
from rdkit import Chem
from retrobiocat_web.mongo.functions import sequence_functions

class EnzymeType(db.Document):
    enzyme_type = db.StringField(max_length=120, unique=True, required=True)
    full_name = db.StringField(default='')
    description = db.StringField(default='')
    other_abbreviations = db.ListField(db.StringField())
    bioinformatics_status = db.StringField(default='Idle')

    def __unicode__(self):
        return self.enzyme_type

    def __str__(self):
        return self.enzyme_type

    meta = {'indexes': ['enzyme_type']}

from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.comments import Comment

class Paper(db.Document):
    doi = db.StringField(unique=True)
    short_citation = db.StringField(default='')
    html = db.StringField(default='')
    owner = db.ReferenceField(User)
    added_by = db.ReferenceField(User)
    edits_by = db.ListField(db.ReferenceField(User))

    title = db.StringField(default='')
    authors = db.ListField(db.StringField())
    journal = db.StringField(default='')
    date = db.DateField()

    status = db.StringField(default='')
    tags = db.ListField(db.StringField())
    reviewed = db.BooleanField(default=False)
    reviewed_by = db.ReferenceField(User)
    has_issues = db.BooleanField(default=False)
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))

    def __unicode__(self):
        return self.short_citation

    def __str__(self):
        return self.short_citation

class Sequence(db.Document):
    enzyme_type = db.StringField(max_length=120, required=True)
    enzyme_name = db.StringField(max_length=120, unique=True, required=True)
    other_names = db.ListField(db.StringField())
    n_tag = db.StringField()
    sequence = db.StringField(default='')
    c_tag = db.StringField()
    sequence_unavailable = db.BooleanField(default=False)

    accession = db.StringField(max_length=20, default='')
    other_identifiers = db.ListField(db.StringField(max_length=20))

    structure = db.BooleanField(default=False)
    pdb = db.StringField(default='')
    mutant_of = db.StringField(default='')
    notes = db.StringField(default='')
    bioinformatics_ignore = db.BooleanField(default=False)

    owner = db.ReferenceField(User)
    added_by = db.ReferenceField(User)
    edits_by = db.ListField(db.ReferenceField(User))
    papers = db.ListField(db.ReferenceField(Paper))

    blast = db.DateTimeField(default=None)
    alignments_made = db.DateTimeField()

    objects_to_update = []

    def update_name(self, new_name):
        new_name = sequence_functions.sanitise_string(new_name)
        q = Sequence.objects(enzyme_name=new_name)
        if len(q) != 0:
            return False, 'Name already exists'
        else:
            print(f"Updating sequence name: {new_name}")

        # update mutants_of
        mutants = Sequence.objects(mutant_of=self.enzyme_name)
        for mut in mutants:
            mut.enzyme_name = new_name
            mut.save()

        # update activity
        acts = Activity.objects(enzyme_name=self.enzyme_name)
        for act in acts:
            act.enzyme_name = new_name
            act.save()

        self.enzyme_name = new_name
        self.save()

        return True, 'Update successful'

    def update_type(self, new_type):
        q = EnzymeType(enzyme_type=new_type)
        if len(q) == 0:
            return False, 'Enzyme type not found'

        # update activity
        acts = Activity.objects(enzyme_name=self.enzyme_name)
        for act in acts:
            act.enzyme_type = new_type
            act.save()

        self.enzyme_type = new_type
        self.blast = None
        self.save()

        return True, 'Update successful'

    def update_sequence(self, seq_string):
        seq_string = sequence_functions.sanitise_sequence(seq_string)

        bad_chars = sequence_functions.sequence_check(seq_string)
        if len(bad_chars) != 0:
            return False, f'Invalid sequence chars: {bad_chars}'

        self.sequence = seq_string
        self.blast = None
        self.save()

    def __unicode__(self):
        return self.enzyme_name

    def __str__(self):
        return self.enzyme_name

    meta = {'indexes': ['enzyme_name', 'enzyme_type']}

class Molecule(db.DynamicDocument):
    smiles = db.StringField()
    mol = db.BinaryField()
    #mfp2_2048 = db.StringField()
    #rdfp2_2048 = db.StringField()

    def get_mol(self):
        return Chem.Mol(self.mol)

    meta = {'indexes': ['smiles']}

class Activity(db.Document):
    enzyme_type = db.StringField()
    enzyme_name = db.StringField()
    reaction = db.StringField()
    short_citation = db.StringField()
    html_doi = db.StringField()

    added_by = db.ReferenceField(User)
    added_by_string = db.StringField()
    paper = db.ReferenceField(Paper)
    edits_by = db.ListField(db.ReferenceField(User))

    substrate_1_smiles = db.StringField()
    substrate_2_smiles = db.StringField()
    product_1_smiles = db.StringField()
    cascade_num = db.StringField()

    temperature = db.StringField()
    ph = db.StringField()
    solvent = db.StringField()
    other_conditions = db.StringField()
    notes = db.StringField()
    reaction_vol = db.StringField()
    formulation = db.StringField()
    biocat_conc = db.StringField()

    kcat = db.FloatField()
    km = db.FloatField()
    mw = db.FloatField()

    substrate_1_conc = db.StringField()
    substrate_2_conc = db.StringField()
    specific_activity = db.FloatField()
    conversion = db.FloatField()
    conversion_time = db.FloatField()

    categorical = db.StringField()
    binary = db.BooleanField()

    selectivity = db.StringField()
    auto_generated = db.BooleanField()

    meta = {'indexes': ['enzyme_name', 'enzyme_type', 'reaction']}

class Tag(db.Document):
    seq = db.StringField()
    n_term = db.BooleanField(default=False)
    c_term = db.BooleanField(default=False)

    def __unicode__(self):
        return self.seq

    def __str__(self):
        return self.seq

class UniRef50(db.Document):
    # These fields are shared with Sequence
    enzyme_name = db.StringField()
    sequence = db.StringField()

    enzyme_type = db.ReferenceField(EnzymeType)
    result_of_blasts_for = db.ListField(db.ReferenceField(Sequence))
    seq_match = db.ReferenceField(Sequence)
    blast_round = db.IntField()
    alignments_made = db.DateTimeField()

    protein_name = db.StringField()
    tax = db.StringField()
    tax_id = db.StringField()

    meta = {'indexes': ['enzyme_type', 'enzyme_name']}

    def __unicode__(self):
        return self.enzyme_name

    def __str__(self):
        return self.enzyme_name

class UniRef90(db.DynamicDocument):
    pass

class Alignment(db.DynamicDocument):
    pass

class SeqSimNet(db.DynamicDocument):
    pass

class SSN_record(db.Document):
    enzyme_type = db.ReferenceField(EnzymeType)
    status = db.StringField(default='Newly created')











