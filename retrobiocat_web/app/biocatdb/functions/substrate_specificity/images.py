from requests.utils import quote
from rdkit.Chem import rdDepictor
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def smitosvg_url(smi, molSize=(150,150), kekulize=True):
    try:
        mol = Chem.MolFromSmiles(smi)
        url = moltosvg_url(mol, molSize=molSize, kekulize=kekulize)
        return url
    except:
        return ''

def moltosvg_url(mol,molSize=(150,150),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    url = "data:image/svg+xml;charset=utf-8," + quote(svg)

    return url