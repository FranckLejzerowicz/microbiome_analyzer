# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import random
random.seed(10)
import numpy as np
import pandas as pd
from os.path import basename
from sklearn.model_selection import train_test_split


def get_train_column(meta, new_meta_pd, meta_vars, train, new_meta, new_meta_ct):
    if train.isdigit() or train.replace('.', '').isdigit():
        train_column = 'TrainTest'
        if train.isdigit():
            train_int = int(train)
            if train_int < (0.1 * new_meta_pd.shape[0]):
                train_perc = 0.1
            else:
                train_perc = train_int / new_meta_pd.shape[0]
        else:
            train_float = float(train)
            if 0 < train_float < 1:
                train_perc = train_float
            else:
                train_column = ''
                print('\t\t\t[SONGBIRD] Float passed as percent of samples for'
                              ' training not valid (must be in range 0-1)')
                return None

        new_meta_vars_pd = new_meta_pd[meta_vars].copy()
        cat_vars = [x for x in meta_vars if str(new_meta_vars_pd[x].dtype) == 'object']
        # print("cat_vars")
        # print(cat_vars)
        if cat_vars:
            new_meta_cat_pd = new_meta_vars_pd[cat_vars].copy()
            new_meta_cat_pd['concat_cols'] = new_meta_cat_pd.apply(
                func=lambda x: '_'.join([str(y) for y in x]), axis=1)
            rep_d = dict(('_'.join([str(i) for i in r]), list(r)) for r in new_meta_cat_pd[cat_vars].values)
            vc = new_meta_cat_pd['concat_cols'].value_counts()
            # print("vc")
            # print(vc)
        if cat_vars and vc.size < new_meta_cat_pd.shape[0] * 0.5:
            if 1 in vc.values:
                vc_in = vc[vc > 1].index.tolist()
                new_meta_cat_pd_in = new_meta_cat_pd.loc[new_meta_cat_pd['concat_cols'].isin(vc_in)]
            else:
                new_meta_cat_pd_in = new_meta_cat_pd.copy()
            X = np.array(new_meta_cat_pd_in.values)
            y = new_meta_cat_pd_in.index.tolist()
            # print("new_meta_cat_pd_in['concat_cols'].unique()")
            # print(new_meta_cat_pd_in['concat_cols'].unique())
            # print("train_perc")
            # print(train_perc)
            if new_meta_cat_pd_in['concat_cols'].unique().size < 2:
            # if train_perc < new_meta_cat_pd_in['concat_cols'].unique().size:
                return None

            _, __, ___, train_samples = train_test_split(
                X, y, test_size=train_perc,
                stratify=new_meta_cat_pd_in['concat_cols'].tolist()
            )
            new_meta_cat_pd[train_column] = [
                'Train' if x in train_samples else
                'Test' for x in new_meta_cat_pd.index
            ]
            if new_meta_ct:
                ct = pd.crosstab(new_meta_cat_pd[train_column],
                                 new_meta_cat_pd['concat_cols']).T.reset_index()
                ct = pd.concat([
                    pd.DataFrame(
                        [rep_d[x] for x in ct['concat_cols']],
                        columns=cat_vars, index=ct.index
                    ),
                    ct[['Train', 'Test']]
                ], axis=1)
                ct.to_csv(new_meta_ct, sep='\t')
        else:
            train_samples = random.sample(
                new_meta_pd.index.tolist(),
                k=int(train_perc * new_meta_pd.shape[0])
            )
        new_meta_vars_pd[train_column] = [
            'Train' if x in train_samples else
            'Test' for x in new_meta_pd.index
        ]
    else:
        train_samples = []
        if train in new_meta_pd.columns:
            if {'Train', 'Test'}.issubset(new_meta_pd[train]):
                train_column = train
                new_meta_vars_pd = new_meta_pd.loc[
                    new_meta_pd[train].isin(['Train', 'Test'])
                ]
            else:
                train_column = ''
                print('\t\t\t[SONGBIRD] Columns passed for training do '
                              'not have "Train" and "Test" factors')
                print('\t\t\t\->', meta)
                return None
        else:
            train_column = ''
            print('\t\t\t[SONGBIRD] Columns passed for training not exists')
            print('\t\t\t\->', meta)
            return None
    if new_meta:
        new_meta_vars_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
    return train_column, train_samples


def get_metadata_train_test(meta_fp, meta_pd, meta_vars, new_meta, train, drop, new_meta_ct):

    if train in meta_pd.columns:
        meta_vars.append(train)

    new_meta_pd = meta_pd[meta_vars]
    new_meta_pd = new_meta_pd.loc[~new_meta_pd.isna().any(1)]
    new_meta_pd = rename_duplicate_columns(new_meta_pd)

    if drop:
        to_remove = pd.concat([
            new_meta_pd[meta_var.lower()].isin(var_drop) for meta_var, var_drop in drop.items()],
            axis=1
        ).any(axis=1)
        new_meta_pd = new_meta_pd.loc[~to_remove]

    train_column, train_samples = get_train_column(
        meta_fp, new_meta_pd, meta_vars, str(train), new_meta, new_meta_ct)
    return train_column, train_samples


def make_train_test_column(meta_fp: str, train_test_dict: dict,
                           meta_pd: pd.DataFrame, dat1: str,
                           dat2: str = '') -> (pd.DataFrame, set):
    train_cols = set()
    meta_tt_pd = meta_pd.copy()
    for dat in [dat1, dat2]:
        if dat in train_test_dict['datasets']:
            for tt, tt_vars in train_test_dict['datasets'][dat].items():
                train_column, train_samples = get_metadata_train_test(
                    meta_fp, meta_pd.set_index('sample_name'), tt_vars, '', train_test_dict['train'], {}, '')
                train_cols.add(tt)
                meta_tt_pd[tt] = [
                    'Train' if x in train_samples else
                    'Test' for x in meta_tt_pd.sample_name.tolist()
                ]
            break
    return meta_tt_pd, train_cols


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
            formula_split = formula.split('C(')[-1].rsplit(')', 1)
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
            formula = formula.lower()
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
                                  main_testing_groups: tuple,
                                  p_perm_tests_min: int,
                                  analysis: str) -> list:
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
        meta_var = meta_pd[variable]
        if variable not in meta_pd_vars:
            print('  [%s] variable %s not in %s' % (analysis, variable, basename(meta)))
            continue
        if meta_var.unique().size > (meta_pd.shape[0] * .8):
            print('  [%s] variable %s from %s not suitable for permanova' % (analysis, variable, basename(meta)))
            continue
        meta_var_vc = meta_var.value_counts()
        meta_var_vc = meta_var_vc[meta_var_vc >= p_perm_tests_min]
        if not meta_var_vc.size >= 2:
            print('  [%s] less than 2 factors in variable %s have >= %s samples' % (analysis, variable, p_perm_tests_min))
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


