# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import sys
import random
import numpy as np
import pandas as pd
from os.path import splitext
from sklearn.model_selection import train_test_split
from microbiome_analyzer._io_utils import check_vals
np.set_printoptions(precision=2, suppress=True)


# def get_subset(
#         case_vals: list,
#         case_var: str,
#         form: str = None
# ) -> str:
#     """
#     Get the current subset, which is the concatenation of:
#      - diversity metric.
#      - metadata variable.
#      - metadata variable's values.
#      - formula.
#
#     Parameters
#     ----------
#     case_vals
#     case_var
#     form
#
#     Returns
#     -------
#
#     """
#     if len(case_vals):
#         case = '%s_%s' % (case_var, '-'.join(
#             [x.replace('<', 'below').replace('>', 'above') for x in case_vals]
#         ))
#     else:
#         case = case_var
#     if form:
#         case = '%s_%s' % (case, form)
#     case = re.sub('[ /()]', '', case.replace('__', '_'))
#     return case


def get_subset(
        case_var: str,
        form: str = None
) -> str:
    """
    Get the current subset, which is the concatenation of:
     - diversity metric.
     - metadata variable.
     - metadata variable's values.
     - formula.

    Parameters
    ----------
    case_var
    form

    Returns
    -------

    """
    case = case_var
    if form:
        case = '%s_%s' % (case, form)
    case = re.sub('[ /()]', '', case.replace('__', '_'))
    return case


# def get_subset_pd(
#         meta_pd: pd.DataFrame,
#         subset: str,
#         variable: str,
#         factors: list
# ) -> pd.DataFrame:
#     """
#     Perform subset.
#
#     Parameters
#     ----------
#     meta_pd
#     subset
#     variable
#     factors
#
#     Returns
#     -------
#
#     """
#     if 'ALL' in subset:
#         new_meta_pd = meta_pd.copy()
#     elif len([x for x in factors if x[0] == '>' or x[0] == '<']):
#         new_meta_pd = meta_pd.copy()
#         for case_val in factors:
#             if case_val[0] == '>':
#                 new_meta_pd = new_meta_pd[
#                     new_meta_pd[variable].astype(float) >= float(case_val[1:])
#                 ].copy()
#             elif case_val[0] == '<':
#                 new_meta_pd = new_meta_pd[
#                     new_meta_pd[variable].astype(float) <= float(case_val[1:])
#                 ].copy()
#     else:
#         new_meta_pd = meta_pd[meta_pd[variable].isin(factors)].copy()
#     return new_meta_pd


def get_subset_pd(
        meta_pd: pd.DataFrame,
        dat: str,
        subset: str,
        subsets: dict
) -> pd.DataFrame:
    """
    Perform subset.

    Parameters
    ----------
    meta_pd : pd.DataFrame
        Metadata table in full
    dat : str
        Name of the dataset
    subset : str
        Name of the samples subset
    subsets : dict
        Conditions for the samples subset

    Returns
    -------
    new_meta_pd : pd.DataFrame
        Metadata table reduced to the samples
    """
    if 'ALL' in subset:
        new_meta_pd = meta_pd.copy()
    else:
        vars_sams = get_sample_subset(meta_pd, dat, subsets)
        if not vars_sams:
            return None
        sams = list(set.intersection(*vars_sams))
        new_meta_pd = meta_pd[meta_pd['sample_name'].isin(sams)].copy()
    return new_meta_pd


def meta_subset(
        meta_pd: pd.DataFrame,
        var: str,
        vals: list
) -> set:
    """Collect the set of samples matching the given variable's factors.

    Parameters
    ----------
    meta_pd : pd.DataFrame
        Metadata table
    var : str
        Name of a metadata variable
    vals : tuple
        Metadata factors for the current variable

    Returns
    -------
    sams : str
        Set of samples for the current metadata variable and factors
    """
    new_meta_pd = meta_pd.copy()
    if 'ALL' in var:
        pass
    elif len([x for x in vals if x[0] == '>' or x[0] == '<']):
        for case_val in vals:
            if case_val[0] == '>':
                new_meta_pd = new_meta_pd[
                    new_meta_pd[var].astype(float) >= float(case_val[1:])
                ].copy()
            elif case_val[0] == '<':
                new_meta_pd = new_meta_pd[
                    new_meta_pd[var].astype(float) <= float(case_val[1:])
                ].copy()
    else:
        new_meta_pd = meta_pd[meta_pd[var].isin(vals)].copy()
    sams = set(new_meta_pd.sample_name)
    return sams


def get_sample_subset(
        meta_pd: pd.DataFrame,
        dat: str,
        vars_vals: dict,
) -> list:
    """

    Parameters
    ----------
    meta_pd : pd.DataFrame
        Matadata table
    dat : str
        Name of the dataset
    vars_vals : dict
        Factors per variable to subset for
        {'var1': ('factor1', 'factor2'), 'var2': ...}

    Returns
    -------
    vars_sams : list
        Samples of the subset
    """
    vars_sams = []
    for var, vals in vars_vals.items():
        if check_vals(meta_pd, dat, 'subset', var, vals):
            continue
        vars_sams.append(meta_subset(meta_pd, var, vals))
    return vars_sams


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
    meta_subset_pd = meta_subset1.merge(
        meta_subset2, on=([sn] + [c for c in comm_cols if c not in diff_cols])
    )
    sorting_col = [sn] + [x for x in meta_subset_pd.columns.tolist() if x != sn]
    meta_subset_pd = meta_subset_pd[sorting_col]
    return meta_subset_pd


def get_cat_vars_and_vc(
        vars_: list,
        vars_pd: pd.DataFrame) -> tuple:
    """

    Parameters
    ----------
    vars_ : list
    vars_pd : pd.DataFrame

    Returns
    -------
    cat_vars : list
    cat_pd : pd.DataFrame
    vc : pd.Series
    rep_d : dict
    """
    cat_vars = [x for x in vars_ if str(vars_pd[x].dtype) == 'object']
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


def make_random_train_test(
        new_meta_pd: pd.DataFrame,
        train_perc: float,
):
    """

    Parameters
    ----------
    new_meta_pd : pd.DataFrame
    train_perc : float

    Returns
    -------
    train_samples
    """
    train_samples = random.sample(
        new_meta_pd.index.tolist(),
        k=int(train_perc * new_meta_pd.shape[0]))
    return train_samples


def make_train_test_from_cat(
        new_meta_pd: pd.DataFrame,
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
    new_meta_pd : pd.DataFrame
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
    train_size = np.floor(cat_pd_in.shape[0] * (1 - train_perc))
    n_classes = cat_pd_in['concat_cols'].unique().size
    if (n_classes >= 2) and (train_size >= n_classes):
        _, __, ___, train_samples = train_test_split(
            X, y, test_size=train_perc,
            stratify=cat_pd_in['concat_cols'].tolist()
        )
        write_cross_tab(
            meta_fp, cat_pd, cat_vars, train_samples, train_col, rep_d)
    else:
        train_samples = make_random_train_test(new_meta_pd, train_perc)
    return train_samples


def rename_duplicate_columns(meta_subset):
    cols = []
    meta_subset_copy = meta_subset.copy()
    for col in meta_subset.columns:
        if col in cols:
            cols.append('%s.%s' % (col, cols.count(col)))
        else:
            cols.append(col)
    meta_subset_copy.columns = cols
    return meta_subset_copy


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
