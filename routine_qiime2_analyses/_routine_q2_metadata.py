# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import sys
import random
random.seed(10)
import numpy as np
import pandas as pd
from os.path import basename, splitext
from sklearn.model_selection import train_test_split


def get_common_columns(
        meta1_pd: pd.DataFrame,
        meta2_pd: pd.DataFrame) -> list:
    """Get these columns names that are identical in two metadata

    Parameters
    ----------
    meta1_pd : pd.DataFrame
        A metadata table
    meta2_pd : pd.DataFrame
        Another metadata table

    Returns
    -------
    common_cols : list
        Columns names in common between
        the two input metadata
    """
    common_cols = set(meta1_pd.columns) & set(meta2_pd.columns)
    common_cols = [x for x in common_cols if x != 'sample_name']
    return common_cols


def get_different_columns(
        meta_subset1: pd.DataFrame,
        meta_subset2: pd.DataFrame,
        common_cols: list) -> list:
    """Find which metadata columns have the
    same name but their content differ.

    Parameters
    ----------
    meta_subset1 : pd.DataFrame
        A metadata table
    meta_subset2 : pd.DataFrame
        Another metadata table
    common_cols : list
        Metadata columns that are in common
        between the two metadata tables

    Returns
    -------
    diff_cols : list
        Metadata columns that are
        different in contents.
    """
    diff_cols = []
    for c in common_cols:
        try:
            meta_col1 = meta_subset1[c].tolist()
            meta_col2 = meta_subset2[c].tolist()
        except:
            print(meta_subset1[c])
            sys.exit(1)
        if meta_col1 != meta_col2:
            diff_cols.append(c)
    return diff_cols


def get_meta_subset(
        meta1_pd: pd.DataFrame,
        meta2_pd: pd.DataFrame,
        sams: set):
    """Merges two metadata tables. This function pay attention to not
    reproduce columns in common when these columns have the same
    contents. Otherwise, these columns are kept in duplicates.

    Parameters
    ----------
    meta1_pd : pd.DataFrame
        A metadata table
    meta2_pd : pd.DataFrame
        Another metadata table
    sams : set
        Samples names in common between
        the two input metadata

    Returns
    -------
    meta_subset : pd.DataFrame
        Metadata formed by the merging
        of the two input metadata
    """
    meta_subset1 = meta1_pd.loc[meta1_pd.sample_name.isin(sams), :].copy()
    meta_subset2 = meta2_pd.loc[meta2_pd.sample_name.isin(sams), :].copy()
    meta_subset1 = rename_duplicate_columns(meta_subset1)
    meta_subset2 = rename_duplicate_columns(meta_subset2)
    comm_cols = get_common_columns(meta_subset1, meta_subset2)
    diff_cols = get_different_columns(meta_subset1, meta_subset2, comm_cols)
    if len(diff_cols):
        meta_subset2.rename(
            columns=dict((c, '%s.copy' % c) for c in diff_cols), inplace=True
        )
    sn = 'sample_name'
    meta_subset = meta_subset1.merge(
        meta_subset2, on=([sn] + [c for c in comm_cols if c not in diff_cols])
    )
    sorting_col = [sn] + [x for x in meta_subset.columns.tolist() if x != sn]
    meta_subset = meta_subset[sorting_col]
    return meta_subset


def get_train_perc_from_numeric(train, new_meta_pd):
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
            print(
                '\t\t\t[SONGBIRD] Float passed as percent of samples for'
                ' training not valid (must be in range 0-1)')
            return None
    return train_perc


def get_cat_vars_and_vc(
        vars: list,
        vars_pd: pd.DataFrame) -> tuple:
    """

    Parameters
    ----------
    vars : list
    vars_pd : pd.DataFrame

    Returns
    -------
    cat_vars : list
    cat_pd : pd.DataFrame
    vc : pd.Series
    rep_d : dict
    """
    cat_vars = [x for x in vars if str(vars_pd[x].dtype) == 'object']
    rep_d = {}
    cat_pd = None
    vc = None
    if cat_vars:
        cat_pd = vars_pd[cat_vars].copy()
        cat_pd['concat_cols'] = cat_pd.apply(
            func=lambda x: '_'.join([str(y) for y in x]), axis=1)
        rep_d = dict(('_'.join([str(i) for i in r]), list(r))
                     for r in cat_pd[cat_vars].values)
        vc = cat_pd['concat_cols'].value_counts()
    return cat_vars, cat_pd, vc, rep_d


def write_cross_tab(
        meta_fp: str,
        cat_pd: pd.DataFrame,
        cat_vars: list,
        train_samples: pd.Series,
        train_col: str,
        rep_d: dict):
    """

    Parameters
    ----------
    meta_fp : str
    cat_pd : pd.DataFrame
    cat_vars : list
    train_samples : pd.Series
    train_col : str
    rep_d : dict
    """
    cat_pd[train_col] = ['Train' if x in train_samples
                         else 'Test' for x in cat_pd.index]
    new_meta_ct = '%s_cv.txt' % splitext(meta_fp)[0]
    ct = pd.crosstab(
        cat_pd[train_col],
        cat_pd['concat_cols']
    ).T.reset_index()
    ct = pd.concat([pd.DataFrame(
        [rep_d[x] for x in ct['concat_cols']],
        columns=cat_vars, index=ct.index),
        ct[['Train', 'Test']]], axis=1)
    ct.to_csv(new_meta_ct, index=False, sep='\t')


def make_train_test_from_cat(
        cat_pd: pd.DataFrame,
        vc: pd.Series,
        train_perc: float,
        meta_fp: str,
        cat_vars: list,
        train_col: str,
        rep_d: dict):
    """

    Parameters
    ----------
    cat_pd : pd.DataFrame
    vc : pd.Series
    train_perc : float
    meta_fp : str
    cat_vars : list
    train_col : str
    rep_d : dict

    Returns
    -------
    train_samples
    """
    if 1 in vc.values:
        vc_in = vc[vc > 1].index.tolist()
        cat_pd_in = cat_pd.loc[cat_pd['concat_cols'].isin(vc_in)]
    else:
        cat_pd_in = cat_pd.copy()
    X = np.array(cat_pd_in.values)
    y = cat_pd_in.index.tolist()
    if cat_pd_in['concat_cols'].unique().size >= 2:
        _, __, ___, train_samples = train_test_split(
            X, y, test_size=train_perc,
            stratify=cat_pd_in['concat_cols'].tolist()
        )
        write_cross_tab(meta_fp, cat_pd, cat_vars,
                        train_samples, train_col, rep_d)
        return train_samples
    return None


def get_train_column(meta_fp, new_meta_pd, meta_vars, train, new_meta):
    if train.isdigit() or train.replace('.', '').isdigit():
        train_column = 'TrainTest'
        print(train)
        print(train)
        print(train)
        print(train)
        train_perc = get_train_perc_from_numeric(train, new_meta_pd)
        print(train_perc)
        print(train_perc)
        print(train_perc)
        print(train_perc)
        vars_pd = new_meta_pd[meta_vars].copy()
        cat_vars, cat_pd, vc, rep_d = get_cat_vars_and_vc(meta_vars, vars_pd)
        if cat_vars and vc.size < cat_pd.shape[0] * 0.5:
            train_samples = make_train_test_from_cat(
                cat_pd, vc, train_perc, meta_fp, cat_vars, train_column, rep_d
            )
        else:
            train_samples = random.sample(
                new_meta_pd.index.tolist(),
                k=int(train_perc * new_meta_pd.shape[0])
            )
        vars_pd[train_column] = [
            'Train' if x in set(train_samples) else
            'Test' for x in new_meta_pd.index
        ]
    else:
        train_samples = []
        if train in new_meta_pd.columns:
            if {'Train', 'Test'}.issubset(new_meta_pd[train]):
                train_column = train
                vars_pd = new_meta_pd.loc[
                    new_meta_pd[train].isin(['Train', 'Test'])
                ]
            else:
                train_column = ''
                print('\t\t\t[SONGBIRD] Columns passed for training do '
                              'not have "Train" and "Test" factors')
                print('\t\t\t\->', meta_fp)
                return None, None
        else:
            train_column = ''
            print('\t\t\t[SONGBIRD] Columns passed for training not exists')
            print('\t\t\t\->', meta_fp)
            return None, None
    if new_meta:
        print('WRITE:', new_meta)
        vars_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
    return train_column, train_samples


def get_metadata_train_test(
        meta_fp, meta_pd, meta_vars, new_meta, train, drop):

    if train in meta_pd.columns:
        meta_vars.append(train)

    new_meta_pd = meta_pd[meta_vars]
    new_meta_pd = new_meta_pd.loc[~new_meta_pd.isna().any(1)]
    new_meta_pd = rename_duplicate_columns(new_meta_pd)

    if drop:
        to_remove = pd.concat([
            new_meta_pd[meta_var.lower()].isin(var_drop)
            for meta_var, var_drop in drop.items()],
            axis=1
        ).any(axis=1)
        new_meta_pd = new_meta_pd.loc[~to_remove]

    train_column, train_samples = get_train_column(
        meta_fp, new_meta_pd, meta_vars, str(train), new_meta)
    return train_column, train_samples


def make_train_test_column(meta_fp: str, train_test_dict: dict,
                           meta_pd: pd.DataFrame, dat1: str,
                           dat2: str = '') -> (pd.DataFrame, set):
    train_cols = set()
    meta_tt_pd = meta_pd.copy()
    for dat in [dat1, dat2]:
        if 'datasets' in train_test_dict and dat in train_test_dict['datasets']:
            for tt, tt_vars in train_test_dict['datasets'][dat].items():
                train_column, train_samples = get_metadata_train_test(
                    meta_fp, meta_pd.set_index('sample_name'), tt_vars,
                    '', train_test_dict['train'], {})
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
            meta_subset_cols.append(
                '%s.%s' % (col, meta_subset_cols.count(col))
            )
        else:
            meta_subset_cols.append(col)
    meta_subset_copy.columns = meta_subset_cols
    return meta_subset_copy


