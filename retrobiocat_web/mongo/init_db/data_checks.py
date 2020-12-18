import numpy as np

def check_is_float(to_test):
    if type(to_test) != float:
        try:
            to_test = float(to_test)
        except:
            return None
    if np.isnan(to_test) == True:
        return None

    return to_test

def check_is_nan(to_test):
    if str(to_test) == 'NaN' or str(to_test) == '' or str(to_test) == ' ':
        return None
    if type(to_test) != str:
        if isinstance(to_test, (int, float, complex)):
            if np.isnan(to_test) == True:
                return None
            else:
                to_test = str(to_test)
    return to_test

def get_mol(smile):
    try:
        mol = Chem.MolFromSmiles(smile)
    except:
        mol = None
    return mol

def get_smile(smile):
    smile = check_is_nan(smile)
    if smile == None:
        return ''
    return smile