from rdkit.Chem import rdDepictor
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions


def apply_smiles_to_filename_check(smile):
    new_string = ''
    smile = str(smile)
    for letter in smile:
        if letter == '#':
            new_string += '=-'
        elif letter == '/':
            new_string += 'fs'
        else:
            new_string += letter
    return new_string

def save_smiles_image(smile, directory=''):
    mol = AllChem.MolFromSmiles(smile)
    filename_smile = apply_smiles_to_filename_check(smile)
    filename = directory + str(filename_smile) + '.png'
    Draw.MolToFile(mol, filename)

def get_images_of_substrates(list_smiles, graph_dir='retrobiocat/retrobiocat_web/app', img_dir='static/mol_images/'):
    directory = graph_dir + '/' + img_dir
    for smile in list_smiles:
        save_smiles_image(smile, directory=directory)

def moltosvg(mol,molSize=(150,150),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg

def rxntosvg(list_rxns, rxnSize=(900,300)):
    list_svgs = []
    for i, rxn in enumerate(list_rxns):
        rdkit_reaction = rxn.rxn
        d = Draw.MolDraw2DSVG(rxnSize[0], rxnSize[1])
        d.DrawReaction(rdkit_reaction)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        list_svgs.append(svg)

    return list_svgs

def smiles_rxn_to_svg(smiles_rxn, rxnSize=(900,300)):
    rxn = rdChemReactions.ReactionFromSmarts(smiles_rxn, useSmiles=True)
    d = Draw.MolDraw2DSVG(rxnSize[0], rxnSize[1])
    d.DrawReaction(rxn)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    return svg