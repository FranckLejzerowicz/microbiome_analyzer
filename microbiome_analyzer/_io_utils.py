# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import sys

import biom
import yaml
import pkg_resources
import pandas as pd

from os.path import basename, dirname, isfile, isdir, splitext
from microbiome_analyzer.core.commands import run_export
from microbiome_analyzer._scratch import to_do, rep

RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources/wol")


def convert_to_biom(data_pd: pd.DataFrame) -> biom.Table:
    data = data_pd.values
    sample_ids = data_pd.columns.tolist()
    observ_ids = data_pd.index.tolist()
    return biom.Table(data, observ_ids, sample_ids, table_id='')


def check_features(data_pd: biom.Table) -> tuple:

    """Parse the feature names anb check whether these are
    all DNA sequences, whether they contain genome IDs,
    or if they have been corrected for " " ot ";"

    Parameters
    ----------
    data_pd : biom.Table
        BIOM table

    Returns
    -------
    genome_ids : dict
        Mapping Genome ID -> corrected ID (taxonomy without ";" or " ")
    dna : bool
        Whether the features are all DNA sequences or not
    correct : bool
        Whether the feature names have been corrected
    """
    # features = data_pd.index.tolist()
    features = data_pd.ids(axis='observation')
    # regex to find the first non-DNA character
    not_dna = re.compile('[^ACGTN].*?')
    # init genome IDs mapping to corrected ID (when also containing taxonomy)
    dna = True
    correct = False
    genome_ids = {}
    for feature in features:
        dna = check_if_not_dna(feature, dna, not_dna)
        if re.match('G\d{9}', feature):
            if ';' in feature:
                correct = True
                feature_corr = feature.replace(';', '|').replace(' ', '')
            else:
                feature_corr = feature
            genome_ids[re.search('G\d{9}', feature).group(0)] = feature_corr
    return genome_ids, dna, correct


def check_if_not_dna(
        feature: str,
        dna: bool,
        not_dna) -> bool:
    """
    Checks whether a feature name do not contain
    any DNA character, as for a SEPP-based sequence
    placement need all features to be DNA sequences.

    Parameters
    ----------
    feature : str
        Feature name
    dna : bool
        Whether the previous feature name is not DNA
    not_dna
        Regex pattern to find the first non-DNA character

    Returns
    -------
    dna : bool
        Whether the current feature name is not DNA
    """
    if dna and bool(not_dna.search(feature)):
        dna = False
    return dna


def check_vals(
        meta_pd: pd.DataFrame,
        dat: str,
        analysis: str,
        var: str,
        vals: list,
        test: str=None
):
    meta_pd_vars = set(meta_pd.columns.tolist())
    if var not in meta_pd_vars:
        print('[%s][%s] no variable %s in "%s"' % (dat, analysis, var, data.dat))
        return True
    else:
        factors = meta_pd[var].unique()
        if vals[0][0] in ['>', '<'] and str(meta_pd[var].dtype) == 'object':
            return True
        factors_common = set(vals) & set(factors.astype(str).tolist())
        if var == test and factors_common == 1:
            print('[%s][%s] subset to %s==["%s"] leave only one category' % (
                dat, analysis, var, '", "'.join(vals)))
            return True
        # if sorted(factors_common) != sorted(vals):
        if not set(factors_common) & set(vals):
            vals_print = ', '.join([val for val in vals[:5]])
            if len(vals) > 5:
                vals_print = '%s, ...' % vals_print
            print('[%s][%s] factors of %s not found (%s)' % (
                dat, analysis, var, vals_print))
            if len([x for x in vals if len(x) == 1]) == len(vals):
                print('[%s] -> check nested list in yaml file.' % dat)
            return True
    return False


def add_q2_type(meta_pd: pd.DataFrame, meta: str, cv: str, tests: list,
                add_q2: bool = False, n=4) -> bool:
    cvs = []
    for test in tests:
        meta_pd = meta_pd.replace({test: dict(
            (x, x.replace('(', '').replace(')', '').replace('/', ''))
            for x in meta_pd[test].astype(str).unique() if str(x) != 'nan'
            and x != x.replace('(', '').replace(')', '').replace('/', ''))})
        cv_pd = meta_pd[test].fillna('NA').value_counts()
        cv_pd = cv_pd[cv_pd >= n]
        if cv_pd.size == 1:
            return True
        if sum(cv_pd) < (n*2):
            return True
        meta_pd = meta_pd.loc[meta_pd[test].isin(cv_pd.index)]
        if add_q2:
            q2types = pd.DataFrame(
                [(['#q2:types'] + ['categorical'] * (meta_pd.shape[1] - 1))],
                columns=meta_pd.columns.tolist())
            meta_pd = pd.concat([q2types, meta_pd])
        cv_pd = pd.DataFrame(cv_pd).reset_index()
        cv_pd = cv_pd.rename(columns={test: 'value'})
        cv_pd['variable'] = test
        cvs.append(cv_pd)
    cv_pd = pd.concat(cvs)
    cv_pd.to_csv(rep(cv), index=False, sep='\t', header=False)
    meta_pd.fillna('NA').to_csv(rep(meta), index=False, sep='\t')
    return False


def subset_dm(
        meta_pd: pd.DataFrame,
        sams: list,
        t1: str,
        t2: str,
        col: str='sample_name'
) -> pd.DataFrame:
    """
    Subset the metadata to the samples in two input distance matrices.

    Parameters
    ----------
    meta_pd : pd.DataFrame
        Metadata table
    sams : list
        All possible samples
    t1 : str
        First distance matrix
    t2 : str
        Second distance matrix
    col : str
        Name of the sample name column (default to 'sample_name')

    Returns
    -------
    new_meta_pd : pd.DataFrame
        Metadata table subset
    """
    with open(rep(t1)) as f:
        for line in f:
            s1 = set(line.strip().split('\t'))
            break
    with open(rep(t2)) as f:
        for line in f:
            s2 = set(line.strip().split('\t'))
            break
    sams = list(s1 & s2 & set(sams))
    new_meta_pd = meta_pd.loc[meta_pd[col].isin(sams)].copy()
    return new_meta_pd


def subset_meta(
        meta_pd: pd.DataFrame,
        sams: list,
        variables: list,
        test: str = '',
        terms: list = [],
        col: str='sample_name'
) -> pd.DataFrame:
    """
    Perform metadata subset.

    Parameters
    ----------
    meta_pd : pd.DataFrame
        Metadata table
    sams : list
    variables : list
    test : str = ''
    terms : list = []
    col : str

    Returns
    -------
    new_meta_pd : pd.DataFrame
        Metadata table subset
    """
    # init new metadata as being a copy of the input, and a columns names set
    new_meta_pd, cols = meta_pd.copy(), set()
    # first subset to the selected samples for the column with sample names
    if variables:
        new_meta_pd = new_meta_pd.loc[new_meta_pd[col].isin(sams)]
        # if variables are given for subset, also keep these in the metadata
        cols.update(variables)
    if test:
        # if a stat test is to be done, also keep the values in the metadata
        cols.add(test)
    if terms:
        # if varibales are terms in models, also keep them in the metadata
        cols.update(terms)
    # subset the metadata variables to all useful columns
    new_meta_pd = new_meta_pd[(['sample_name'] + list(cols))]
    return new_meta_pd


def get_wol_tree(wol_tree: str) -> tuple:
    """
    """
    if wol_tree == 'resources/wol_tree.nwk':
        wol_tree_fp = '%s/wol_tree.nwk' % RESOURCES
        return wol_tree_fp, ''
    if wol_tree:
        sys.exit('"%s" does not exist\nExiting...' % wol_tree)
    if not wol_tree.endswith('.nwk'):
        if wol_tree.endswith('qza'):
            wol_tree_nwk = '%s.nwk' % splitext(wol_tree)[0]
            if isfile(wol_tree_nwk):
                print('"%s" found (export of "%s"?)' % (wol_tree_nwk, wol_tree))
            cmd = run_export(wol_tree, wol_tree_nwk, 'phylogeny')
            return wol_tree_nwk, cmd
        else:
            sys.exit('"%s" is not a .nwk (tree) file or not a'
                     ' qiime2 Phylogeny artefact\nExiting...' % wol_tree)
    else:
        return wol_tree, ''


def get_sepp_tree(sepp_tree: str) -> str:
    """
    Get the full path of the reference database for SEPP.

    :param sepp_tree: database to use.
    :return: path of the reference database for SEPP.
    """
    if not sepp_tree or not isfile(sepp_tree):
        sys.exit('"%s" does not exist\nExiting...' % sepp_tree)
    if not sepp_tree.endswith('qza'):
        sys.exit('"%s" is not a qiime2 Phylogeny\nExiting...' % sepp_tree)
    if basename(sepp_tree) in ['sepp-refs-silva-128.qza',
                               'sepp-refs-gg-13-8.qza']:
        return sepp_tree
    else:
        message = (
                '%s is not:\n'
                '- "sepp-refs-silva-128.qza"\n'
                '- "sepp-refs-gg-13-8.qza"\n'
                'Download: https://docs.qiime2.org/<YYYY.MM>/data-resources/\n' 
                'Exiting...' % sepp_tree)
        sys.exit(mes)


def get_taxonomy_classifier(classifier: str) -> str:
    """
    Get the full path of the reference taxonomic classifier.

    Parameters
    ----------
    classifier : str
        Database to use.

    Returns
    -------
    classifier : str
        Database to use and that exists.

    """
    if not isfile(classifier):
        print('%s does not exist\nExiting...' % classifier)
        sys.exit(0)
    if not classifier.endswith('qza'):
        print('%s is not a qiime2 artefact\nExiting...' % classifier)
        sys.exit(0)

    classifiers = [
        'silva-132-99-nb-classifier.qza',
        'silva-132-99-515-806-nb-classifier.qza',
        'gg-13-8-99-nb-classifier.qza',
        'gg-13-8-99-515-806-nb-classifier.qza'
    ]
    if basename(classifier) in classifiers:
        return classifier
    else:
        print('%s is not:\n%s'
              'Download: https://docs.qiime2.org/2020.2/data-resources/'
              '#taxonomy-classifiers-for-use-with-q2-feature-classifier)\n'
              'Exiting...' % (classifier, '\n- %s' % '\n- '.join(classifiers)))
        sys.exit(0)


def parse_g2lineage() -> dict:
    """

    Returns
    -------

    """
    g2lineage_fp = '%s/g2lineage.txt' % RESOURCES
    g2lineage = {}
    for line in open(g2lineage_fp).readlines():
        line_split = line.strip().split('\t')
        g2lineage[line_split[0]] = line_split[1]
    return g2lineage


def write_filtered_tsv(tsv_out: str, tsv_pd: pd.DataFrame) -> None:
    """

    Parameters
    ----------
    tsv_out
    tsv_pd

    Returns
    -------

    """
    tsv_sams_col = tsv_pd.reset_index().columns[0]
    tsv_pd = tsv_pd.reset_index().rename(
        columns={tsv_sams_col: 'Feature ID'}).set_index('Feature ID')
    tsv_pd.reset_index().to_csv(tsv_out, index=False, sep='\t')


def make_fp_dir(fp):
    if not isdir(rep(dirname(fp))):
        os.makedirs(rep(dirname(fp)))
