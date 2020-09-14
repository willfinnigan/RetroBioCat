from rdkit.Chem import AllChem, Draw
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from requests.utils import quote

from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs

import time


def moltosvg(mol,molSize=(300,300)):
    """
    this func ra ra ra

    Args:
        mol:
        molSize:

    Returns:

    """

    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    opts = drawer.drawOptions()
    opts.padding = 0.1
    opts.addStereoAnnotation = True
    drawer.SetFontSize(1.0)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    return svg

def smile_to_svg_url(smiles, size=(300,300)):
    mol = AllChem.MolFromSmiles(smiles)
    svg = moltosvg(mol, molSize=size)
    url = "data:image/svg+xml;charset=utf-8," + quote(svg)
    return url

def morgan_fingerprint_vis(product_smiles, sim_smiles):
    mol = Chem.MolFromSmiles(sim_smiles)
    refmol = Chem.MolFromSmiles(product_smiles)

    fp1 = SimilarityMaps.GetMorganFingerprint(mol)
    fp2 = SimilarityMaps.GetMorganFingerprint(refmol)

    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, SimilarityMaps.GetMorganFingerprint)
    return DataStructs.FingerprintSimilarity(fp1, fp2)
