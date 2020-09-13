
from retrobiocat_web.retro.generation.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine, RuleApplicator, ReactionSelector, GraphManipulator
from pathlib import Path
import pandas as pd
import pickle
import sqlite3
import uuid

from rdkit.Chem import AllChem, rdmolops
from retrobiocat_web.retro.generation.node_analysis import rdkit_smile
from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.retro.generation.retrosynthesis_engine import aizynthfinder_actions


class AIZynth_RuleApplicator(RuleApplicator):

    def __init__(self, network):
        super().__init__(network)
        self.action_applier = aizynthfinder_actions.ActionApplier()

    def run(self, smile, rxns, graph, explicit_hydrogens=False):
        precursor_dict = self.apply_rules(smile, rxns)
        precursor_dict = self._remove_precursors_already_in_graph(graph, smile, precursor_dict)

        if self.network.settings['remove_simple'] == True:
            for name in precursor_dict:
                precursor_dict[name] = self._remove_simple_precursors(precursor_dict[name])

        return precursor_dict

    def get_rxns(self, smi):
        if self.action_applier.policy_model == None:
            self.action_applier.load_model()
        reactions = self.action_applier.get_actions(smi)

class AIZynth_ReactionSelector(ReactionSelector):

    def __init__(self, network, graphPruner):
        super().__init__(network, graphPruner)

class AIZynth_GraphManipulator(GraphManipulator):

    def __init__(self, network):
        super().__init__(network)

class AIZynth_RetrosynthesisEngine(RetrosynthesisEngine):

    def __init__(self, network):
        super().__init__(network)
        self.reactionSelector = AIZynth_ReactionSelector(network, self.graphPruner)
        self.ruleApplication = AIZynth_RuleApplicator(network)
        self.graphManipulator = AIZynth_GraphManipulator(network)

    def single_step(self, smile, rxns, graph, disallowedProducts=None):
        if self._should_rules_be_applied(smile, graph) == False:
            return [],[]

        rxnsSubstratesDict = self.ruleApplication.run(smile, rxns, graph)
        listProducts, listReactions = self.graphManipulator.add_nodes_to_graph(rxnsSubstratesDict, smile, graph)
        listProducts, listReactions = self.reactionSelector.remove_disallowed_products(disallowedProducts, listProducts, listReactions)

        return listProducts, listReactions