import pubchempy

def get_pubchem_compound_from_smiles(smiles):
    try:
        list_compounds = pubchempy.get_compounds(smiles, namespace='smiles')
    except:
        print('Pubchempy: smiles not found')
        return None

    if len(list_compounds) > 0:
        c = list_compounds[0]
        return c

    else:
        print('Pubchempy: list cids was zero')
        return None