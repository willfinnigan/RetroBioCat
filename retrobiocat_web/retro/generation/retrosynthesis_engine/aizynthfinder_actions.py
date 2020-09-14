""" This code has been modified from the AIZynthfinder repository at https://github.com/MolecularAI/aizynthfinder """

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import functools
from tensorflow.keras.metrics import top_k_categorical_accuracy
from tensorflow.keras.models import load_model
import tensorflow
import logging
import os
from pathlib import Path
from retrobiocat_web.retro.rdchiral.main import rdchiralReaction

os.environ['KMP_DUPLICATE_LIB_OK']='True'

top10_acc = functools.partial(top_k_categorical_accuracy, k=10)
top10_acc.__name__ = "top10_acc"

top50_acc = functools.partial(top_k_categorical_accuracy, k=50)
top50_acc.__name__ = "top50_acc"

CUSTOM_OBJECTS = {"top10_acc": top10_acc, "top50_acc": top50_acc}

tf_logger = tensorflow.get_logger()
tf_logger.setLevel(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

data_folder = str(Path(__file__).parents[3]) + '/retro/data/aizynthfinder'

class LocalKerasModel:
    def __init__(self, filename):
        self.model = load_model(filename, custom_objects=CUSTOM_OBJECTS)
        try:
            self._model_dimensions = int(self.model.input.shape[1])
        except AttributeError:
            self._model_dimensions = int(self.model.input[0].shape[1])

    def __len__(self):
        return self._model_dimensions

    def predict(self, input_):
        return self.model.predict(input_)

class ActionApplier():

    def __init__(self):
        self.policy_model = None
        self.templates = None
        self.cutoff_cumulative = 0.995
        self.cutoff_number = 50
        self.template_column = 'retro_template'

    def load_model(self):
        if self.policy_model == None:
            policy_path = data_folder + '/uspto_model.hdf5'
            self.policy_model = LocalKerasModel(policy_path)
        if self.templates == None:
            templates_path = data_folder + '/uspto_templates.hdf5'
            self.templates = pd.read_hdf(templates_path, "table")

    def get_actions(self, smi):
        reactions = []
        priors = []

        mol = Chem.MolFromSmiles(smi)

        all_transforms_prop = self._predict(mol)
        probable_transforms_idx = self._cutoff_predictions(all_transforms_prop)

        possible_moves = self.templates.iloc[probable_transforms_idx]
        probs = all_transforms_prop[probable_transforms_idx]

        priors.extend(probs)
        for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
            reaction = {}
            metadata = dict(move)
            del metadata[self.template_column]
            metadata["policy_probability"] = float(probs[idx])
            metadata["template_code"] = move_index

            reaction = {'smarts': move[self.template_column],
                        'metadata': metadata,
                        'prior': priors[idx]}

            reactions.append(reaction)

        return reactions

    def get_rxns(self, smile):
        if self.policy_model == None:
            self.load_model()

        reactions = self.get_actions(smile)
        rxns = {}
        for reaction in reactions:
            name = f"Chem_{reaction['metadata']['classification']}"
            if name not in rxns:
                rxns[name] = []
            rxns[name].append(rdchiralReaction(reaction['smarts']))
        return rxns

    def _predict(self, mol):
        fingerprint = self._get_fingerprint(mol, 2, nbits=len(self.policy_model))
        fp_arr = fingerprint.reshape([1, len(self.policy_model)])
        return np.array(self.policy_model.predict(fp_arr)).flatten()

    def _get_fingerprint(self, rd_mol, radius, nbits=None):
        """
        Returns the Morgan fingerprint of the molecule
        """
        if nbits:
            key = (radius, nbits)
        else:
            key = (radius,)

        bitvect = AllChem.GetMorganFingerprintAsBitVect(rd_mol, *key)
        array = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(bitvect, array)
        return array

    def _cutoff_predictions(self, predictions):
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """
        sortidx = np.argsort(predictions)[::-1]
        cumsum = np.cumsum(predictions[sortidx])
        if any(cumsum >= self.cutoff_cumulative):
            maxidx = np.argmin(cumsum < self.cutoff_cumulative)
        else:
            maxidx = len(cumsum)

        maxidx = min(maxidx, self.cutoff_number) or 1
        return sortidx[:maxidx]

aizynth_action_applier = ActionApplier()

if __name__ == '__main__':
    reactions = aizynth_action_applier.get_actions('CCCCCO')
    print(reactions)

