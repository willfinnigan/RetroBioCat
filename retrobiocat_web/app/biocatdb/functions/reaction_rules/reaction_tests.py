import yaml
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.retro.rdchiral.main import rdchiralReaction
from retrobiocat_web.retro.generation.network_generation.network import Network
from retrobiocat_web.retro.generation.retrosynthesis_engine.retrosynthesis_engine import RuleApplicator
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile

class ReactionTester(object):

    def __init__(self):
        self.state = 'danger'
        self.msgs = {'danger': 'Tests failed',
                     'success': 'Tests passed',
                     'warning': 'Tests passed with warnings'}
        self.issues = []

    def get_msg(self):
        return self.msgs[self.state]

    def run(self, selection, name, smarts, cofactors, positive_tests, negative_tests,
            type, experimental, two_step, requires_absence_of_water):
        self.state = 'success'

        self._test_name(name, selection)
        self._test_cofactors_dict(cofactors)
        list_rxns = self._test_smarts(smarts)

        if list_rxns is not None:
            self._positive_tests(positive_tests, list_rxns)
            self._negative_tests(negative_tests, list_rxns)


    def _test_name(self, rxn_name, selection):
        if rxn_name == '' or rxn_name == ' ' or rxn_name == 'Empty template':
            self.issues.append('Reaction must have a name')
            self.state = 'danger'

        if rxn_name != selection:
            if len(Reaction.objects(name=rxn_name)) != 0:
                self.issues.append('Reaction name already exists')
                self.state = 'danger'

    def _test_cofactors_dict(self, cofactors):

        try:
            cofactors_dict = yaml.load(cofactors, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load cofactors yaml')
            return

        if cofactors_dict is None:
            self.state = 'danger'
            self.issues.append('At least one entry must be made for enzymes/cofactors')
            return

        for enz in list(cofactors_dict.keys()):
            if enz not in list(EnzymeType.objects().distinct('enzyme_type')):
                self.state = 'danger'
                self.issues.append(f'Enzyme type {enz} is not defined')

            elif type(cofactors_dict[enz]) != dict:
                self.state = 'danger'
                self.issues.append(f'Cofactors for {enz} are not structured as a dictionary')

            elif 'cofactors_plus' not in cofactors_dict[enz] or 'cofactors_minus' not in cofactors_dict[enz]:
                self.state = 'danger'
                self.issues.append(f'Cofactors must be defined')

            elif type(cofactors_dict[enz]['cofactors_plus']) != list or type(cofactors_dict[enz]['cofactors_minus']) != list:
                self.state = 'danger'
                self.issues.append(f'Cofactors must be given as lists')

    def _test_smarts(self, smarts):
        try:
            smarts_list = yaml.load(smarts, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load smarts yaml')
            return

        try:
            rxn_list = []
            for sma in smarts_list:
                rxn_list.append(rdchiralReaction(sma))
        except:
            self.state = 'danger'
            self.issues.append('Could not load reactions')
            return

        return rxn_list

    def _positive_tests(self, positive_tests, list_rxns):
        empty_network = Network()
        rule_applicator = RuleApplicator(empty_network)
        rxns = {'tests': list_rxns}

        try:
            positive_tests = yaml.load(positive_tests, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load positive tests yaml')
            return

        for test_product in positive_tests:
            try:
                rdkit_smile(test_product)
            except:
                self.state = 'danger'
                self.issues.append(f'Positive test SMILE: {test_product} not accepted by rdkit')
                return

        try:
            for test_product in positive_tests:
                reaction_outcomes = self._apply_reactions(empty_network, rule_applicator, test_product, rxns)
                if len(reaction_outcomes) == 0:
                    self.state = 'danger'
                    self.issues.append(f'Reaction not in outcomes for tested positive product: {test_product}')
        except:
            self.state = 'danger'
            self.issues.append('Problem running positive tests')
            return

        return True

    def _negative_tests(self, negative_tests, list_rxns):
        empty_network = Network()
        rule_applicator = RuleApplicator(empty_network)
        rxns = {'tests': list_rxns}

        try:
            negative_tests = yaml.load(negative_tests, Loader=yaml.FullLoader)
        except:
            self.state = 'danger'
            self.issues.append('Could not load negative tests yaml')
            return

        for test_product in negative_tests:
            try:
                rdkit_smile(test_product)
            except:
                self.state = 'danger'
                self.issues.append(f'Negative test SMILE: {test_product} not accepted by rdkit')
                return

        for test_product in negative_tests:
            reaction_outcomes = self._apply_reactions(empty_network, rule_applicator, test_product, rxns)
            if len(reaction_outcomes) != 0:
                self.state = 'danger'
                self.issues.append(f'Reaction should not be outcomes for tested negative product: {test_product}')

        try:
            for test_product in negative_tests:
                reaction_outcomes = self._apply_reactions(empty_network, rule_applicator, test_product, rxns)
                if len(reaction_outcomes) != 0:
                    self.state = 'danger'
                    self.issues.append(f'Reaction should not be outcomes for tested negative product: {test_product}')
        except:
            self.state = 'danger'
            self.issues.append('Problem running negative tests')
            return

        return True

    def _apply_reactions(self, empty_network, rule_applicator, test_product, rxns):
        reaction_outcomes = []
        empty_network.settings["combine_enantiomers"] = True
        reaction_outcomes += rule_applicator.apply_rules(test_product, rxns)
        empty_network.settings["combine_enantiomers"] = False
        reaction_outcomes += rule_applicator.apply_rules(test_product, rxns)
        return reaction_outcomes




