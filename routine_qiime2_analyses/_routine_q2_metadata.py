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
                    print('  [%s] factors (%s) of variable %s not in metadata:\n\t%s ' % (
                        analysis, factors_print, variable, basename(meta)))
                    if len([x for x in factors_print if len(x) == 1]) == len(factors_print):
                        print('    -> this is certainly due to non-nested list in the yaml file.')
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
        terms = [x.lower() for x in set(re.split('[*+-/)(]+', formula.strip('"').strip("'")))]
        for variable in terms:
            if variable not in meta_pd_vars:
                print('  [%s] variable %s not in %s' % (analysis, variable, basename(meta)))
                to_pop.append(variable)
    if to_pop:
        for pop in to_pop:
            formulas.pop(pop)
    return formulas


def check_metadata_models(meta: str, meta_pd: pd.DataFrame,
                          songbird_models: dict) -> dict:
    """
    :param meta: metadata file name.
    :param meta_pd: metadata table.
    :param songbird_models: groups to factor in the songbird model (need metadata presence).
    :param analysis: current qiime2 analysis.
    :return: checked songbird models.
    """
    models = {}
    for model, formula_ in songbird_models.items():
        vars = set()
        drop = {}
        formula = formula_.strip('"').strip("'")
        if formula.startswith('C('):
            formula_split = formula.split('C(')[-1].split(')')
            formula_split_c = formula_split[0].split(',')[0].strip().strip()
            formula = 'C(%s)' % formula_split[0].replace(formula_split_c, formula_split_c.lower())
            formula += formula_split[1].lower()
            if 'Diff' in formula:
                levels = {formula_split_c.lower(): [
                    x.strip().strip('"').strip("'")
                    for x in formula.split("levels=['")[-1].split("']")[0].split(",")
                ]}
            elif 'Treatment(' in formula:
                levels = {formula_split_c.lower(): [
                    formula.split("Treatment('")[-1].split("')")[0]
                ]}
            vars.add(formula_split_c.lower())
            vars.update(set([x.lower() for x in re.split('[+/:*]', formula_split[1]) if x]))
        else:
            formula_split = re.split('[+/:*]', formula)
            vars.update(set([x.lower() for x in formula_split]))
            levels = {}

        common_with_md = set(meta_pd.columns.str.lower().values) & vars
        if sorted(vars) != sorted(common_with_md):
            only_formula = sorted(vars ^ common_with_md)
            print('\t\tSongbird formula term(s) missing in metadata:\n\t\t  %s\n\t\t  [not used]: %s=%s\n\t\t%s' % (
                ', '.join(only_formula), model, formula, meta))
            continue

        if levels:
            # print('levels', levels)
            # print(meta_pd['sex'].unique())
            levels_set = sorted([x for x in meta_pd[formula_split_c.lower()].unique() if str(x) != 'nan'])
            # print('levels_set', levels_set)
            if 'Diff' in formula:
                cur_levels = levels[formula_split_c.lower()]
                common_levels = set(levels_set) & set(cur_levels)
                only_meta = set(levels_set) ^ common_levels
                only_model = set(cur_levels) ^ common_levels
                if len(only_model):
                    print('\t\tSongbird formula "Diff" factors(s) missing in metadata "%s": %s\n\t\t%s' % (
                        formula_split_c, list(only_model), meta))
                    continue
                if len(only_meta):
                    drop[formula_split_c.lower()] = list(only_meta)
                    print('\t\tSongbird formula "Diff" factors(s) incomplete for metadata "%s":\n\t\t'
                          '  -> skipping samples with %s\n\t\t%s' % (formula_split_c, list(only_meta), meta))
            elif 'Treatment(' in formula:
                levels = {formula_split_c.lower(): formula.split("Treatment('")[-1].split("')")[0]}
                if levels[formula_split_c.lower()] not in levels_set:
                    print('\t\tSongbird formula "Treatment" factors(s) missing in metadata "%s":\n\t\t  %s\n\t\t%s' % (
                        formula_split_c, levels, meta))
                    continue
        models[model] = [formula, vars, drop]
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
    main_testing = []
    for variable in main_testing_groups:
        if variable not in meta_pd_vars:
            print('  [%s] variable %s not in %s' % (analysis, variable, basename(meta)))
            continue
        if meta_pd[variable].unique().size > (meta_pd.shape[0] * .8):
            print('  [%s] variable %s from %s not suitable for permanova' % (analysis, variable, basename(meta)))
            continue
        main_testing.append(variable)
    return main_testing


def rename_duplicate_columns(meta_subset):
    meta_subset_cols = []
    meta_subset_copy = meta_subset.copy()
    for col in meta_subset.columns:
        if col in meta_subset_cols:
            meta_subset_cols.append('%s.%s' % (col, meta_subset_cols.count(col)))
        else:
            meta_subset_cols.append(col)
    meta_subset_copy.columns = meta_subset_cols
    return meta_subset_copy


