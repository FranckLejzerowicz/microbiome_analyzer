# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import pandas as pd
from os.path import basename


def check_metadata_cases_dict(meta: str, meta_pd: pd.DataFrame,
                              cases_dict: dict) -> dict:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param cases_dict: cases to check for possible use (need metadata presence)
    :return: checked cases.
    """
    meta_pd_vars = set(meta_pd.columns.tolist())
    for variable, factors_lists in cases_dict.items():
        if variable not in meta_pd_vars:
            print('  ! [not analyzed] variable %s not in %s' % (variable, basename(meta)))
            cases_dict.pop(variable)
        else:
            factors = set(meta_pd[variable].unique().tolist())
            for factors_list in factors_lists:
                factors_common = set(factors_list) & factors
                if sorted(factors_common) != sorted(factors):
                    print('  ! [not analyzed] factors of variable %s not in %s:' % (variable, basename(meta)))
                    for factor in factors:
                        print('     - %s' % factor)
                    cases_dict.pop(variable)
    return cases_dict


def check_metadata_formulas(meta: str, meta_pd: pd.DataFrame,
                            formulas: dict) -> dict:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param formulas: formulas to check for possible use (need metadata presence)
    :return: checked formulas.
    """
    meta_pd_vars = set(meta_pd.columns.tolist())
    for formula_name, formula in formulas.items():
        terms = set(re.split('[*+-/]+', formula.strip('"').strip("'")))
        for variable in terms:
            if variable not in meta_pd_vars:
                print('  ! [not analyzed] variable %s not in %s' % (variable, basename(meta)))
                formulas.pop(variable)
    return formulas
