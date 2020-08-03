import pandas as pd
import numpy as np


def get_level(amount, levels):
    cat = 'None'
    current_cat_amount = 0

    for category, minimum in levels.items():
        if amount >= minimum:
            if minimum > current_cat_amount:
                cat = category
                current_cat_amount = minimum
    return cat

def binary_from_category(data_dict):
    if data_dict.get('binary', '') == '':
        if data_dict.get('categorical', '') != '':
            if data_dict['categorical'] == 'None':
                data_dict['binary'] = 0
            else:
                data_dict['binary'] = 1
    return data_dict

def category_from_conversion(data_dict, levels):
    if type(data_dict.get('conversion','')) != str:
        if data_dict.get('categorical', '') == '':
            data_dict['categorical'] = get_level(data_dict['conversion'], levels)
    return data_dict

def category_from_sa(data_dict, levels):

    if type(data_dict.get('specific_activity','')) != str:
        if data_dict.get('categorical', '') == '':
            data_dict['categorical'] = get_level(data_dict['specific_activity'], levels)
    return data_dict

def sa_from_kinetics(data_dict, S):

    if type(data_dict.get('km', '')) != str and type(data_dict.get('kcat', '')) != str and type(data_dict.get('mw', '')) != str:
        if type(data_dict.get('specific_activity','')) == str:
            data_dict['specific_activity'] = calc_sa(data_dict['kcat'], data_dict['km'], data_dict['mw'], DEFAULT_CONC)
            data_dict['substrate_1_conc'] = f"{S} mM - calculated from kinetics"
    return data_dict

def calc_sa(kcat, km, mw, s):
    if kcat == 0 or s == 0 or km == 0:
        umol_min_mg = 0
    else:
        umols_enz = 1000/mw
        vmax = kcat*umols_enz
        umol_min_mg = vmax * (s/(s+km))
        umol_min_mg = round(umol_min_mg, 2)
    return umol_min_mg


CATS_CONVERSION = {'High': 60,
                   'Medium': 10,
                   'Low': 0.01,
                   'None': 0}

CATS_SA = {'High': 1,
           'Medium': 0.25,
           'Low': 0.05,
           'None': 0}

DEFAULT_CONC = 10

def remove_empty_columns(data_dict):
    for key in list(data_dict.keys()):
        if data_dict[key] == '' and key != '_id':
            data_dict.pop(key)
    return data_dict

def make_sure_binary_is_int(data_dict):
    try:
        if 'binary' in data_dict:
            data_dict['binary'] = int(data_dict['binary'])
    except:
        if 'binary' in data_dict:
            data_dict.pop('binary')
    return data_dict

def set_activities(data_dict):
    data_dict = sa_from_kinetics(data_dict, DEFAULT_CONC)
    data_dict = category_from_sa(data_dict, CATS_SA)
    data_dict = category_from_conversion(data_dict, CATS_CONVERSION)
    data_dict = binary_from_category(data_dict)
    data_dict = make_sure_binary_is_int(data_dict)
    return data_dict

