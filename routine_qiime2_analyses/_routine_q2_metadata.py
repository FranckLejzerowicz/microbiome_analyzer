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
                              cases_dict: dict, analysis: str) -> dict:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param cases_dict: cases to check for possible use (need metadata presence).
    :param analysis: current qiime2 analysis.
    :return: checked cases.
    """
    to_pop = set()
    meta_pd_vars = set(meta_pd.columns.tolist())
    for variable, factors_lists in cases_dict.items():
        if variable == 'ALL':
            continue
        if variable not in meta_pd_vars:
            print('  [%s] variable %s not in %s' % (analysis, variable, basename(meta)))
            to_pop.add(variable)
        else:
            factors = set(meta_pd[variable].unique().astype(str).tolist())
            for factors_list in factors_lists:
                if factors_list[0][0] in ['>', '<']:
                    continue
                factors_common = set(factors_list) & factors
                if sorted(factors_common) != sorted(factors_list):
                    factors_print = ', '.join([factor for factor in factors_list[:5]])
                    if len(factors_list) > 5:
                        factors_print = '%s, ...' % factors_print
                    print('  [%s] factors of variable %s not in %s (%s)' % (
                        analysis, variable, basename(meta), factors_print))
                    to_pop.add(variable)
    if to_pop:
        for pop in to_pop:
            cases_dict.pop(pop)
    return cases_dict


def check_metadata_formulas(meta: str, meta_pd: pd.DataFrame,
                            formulas: dict, analysis: str) -> dict:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param formulas: formulas to check for possible use (need metadata presence).
    :return: checked formulas.
    """
    to_pop = []
    meta_pd_vars = set(meta_pd.columns.tolist())
    for formula_name, formula in formulas.items():
        terms = [x.lower() for x in set(re.split('[*+-/]+', formula.strip('"').strip("'")))]
        for variable in terms:
            if variable not in meta_pd_vars:
                print('  [%s] variable %s not in %s' % (analysis, variable, basename(meta)))
                to_pop.append(variable)
    if to_pop:
        for pop in to_pop:
            formulas.pop(pop)
    return formulas


def check_metadata_models(meta: str, meta_pd: pd.DataFrame,
                          songbird_models: dict, analysis: str) -> dict:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param songbird_models: groups to factor in the songbird model (need metadata presence).
    :param analysis: current qiime2 analysis.
    :return: checked songbird models.
    """
    meta_pd_vars = [x.lower() for x in set(meta_pd.columns.tolist())]
    models = {}
    for model, formula_ in songbird_models.items():
        formula = formula_.strip('"').strip("'")
        if formula.startswith('C('):
            formula_split = [formula.split('C(')[-1].split(',')[0].strip().strip()]
            formula = formula.replace('C(%s' % formula_split[0], 'C(%s' % formula_split[0].lower())
            if 'Diff' in formula:
                levels = formula.split("levels=['")[-1].split("']")[0].split("','")
            elif 'Treatment(' in formula:
                levels = formula.split("Treatment('")[-1].split("')")[0]
        else:
            formula = formula.lower()
            formula_split = [x.lower() for x in re.split('[+/:*]', formula)]
            levels = []
        common_with_md = set(meta_pd_vars) & set(formula_split)
        if sorted(set(formula_split)) != sorted(common_with_md):
            only_formula = sorted(set(formula_split) ^ common_with_md)
            print('Songbird formula term(s) missing in metadata:\n  %s\n  [not used]: %s=%s' % (
                ', '.join(sorted(only_formula)), model, formula))
            continue
        if len(levels):
            levels_set = [str(x) for x in meta_pd[formula_split[0]].unique()]
            if 'Diff' in formula:
                if not set(levels).issubset(levels_set):
                    not_in_levels = set(levels_set).difference(levels)
                    print('Songbird formula factors(s) missing in metadata "%s":\n  %s' % (
                        formula_split[0], not_in_levels))
                    continue
            elif 'Treatment(' in formula:
                levels = formula.split("Treatment('")[-1].split("')")[0]
                if levels not in levels_set:
                    print('Songbird formula factors(s) missing in metadata "%s":\n  %s' % (
                        formula_split[0], levels))
                    continue
        models[model] = formula
    return models


def check_metadata_testing_groups(meta: str, meta_pd: pd.DataFrame,
                                  main_testing_groups: tuple, analysis: str) -> list:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param main_testing_groups: groups to test for (need metadata presence).
    :param analysis: current qiime2 analysis.
    :return: checked testing groups.
    """
    meta_pd_vars = set(meta_pd.columns.tolist())
    main_testing = [variable for variable in main_testing_groups if variable in meta_pd_vars]
    for variable in main_testing_groups:
        if variable not in main_testing:
            print('  [%s] variable %s not in %s' % (analysis, variable, basename(meta)))
    return main_testing