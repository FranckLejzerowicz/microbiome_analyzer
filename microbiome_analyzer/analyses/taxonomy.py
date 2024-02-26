# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import biom
import pandas as pd
from os.path import isdir

from microbiome_analyzer._inputs import read_meta_pd
from microbiome_analyzer._io_utils import (
    get_taxonomy_classifier, parse_g2lineage)
from microbiome_analyzer._scratch import io_update, to_do, rep
from microbiome_analyzer.core.commands import (
    write_fasta, write_taxonomy_sklearn, run_export, run_import, run_summary)


def get_tax_tables(tax_fp: str) -> tuple:
    """

    Parameters
    ----------
    tax_fp : str

    Returns
    -------
    tax_pd : pd.DataFrame
    split_taxa_pd : pd.DataFrame
    """
    # read taxonomy with features as index, format and collect features IDs list
    tax_pd = pd.read_csv(tax_fp, header=0, sep='\t', dtype=str)
    tax_pd.rename(columns={tax_pd.columns[0]: 'Feature ID'}, inplace=True)
    # perform taxonomic split on the Taxon list and give the Feature as index
    split_taxa_pd = get_split_taxonomy(tax_pd.Taxon.tolist())
    split_taxa_pd.index = tax_pd['Feature ID'].tolist()
    split_taxa_pd.index.name = 'Feature ID'
    return tax_pd, split_taxa_pd


def get_split_taxonomy(taxa, extended=False, taxo_sep=';'):
    """

    Parameters
    ----------
    taxa
    extended
    taxo_sep

    Returns
    -------

    """
    split_lens = set()
    split_taxa = []
    for tdx, taxon in enumerate(taxa):
        if str(taxon) == 'nan':
            split_lens.add(1)
            split_taxa.append(pd.Series(['Unassigned']))
        else:
            taxon_split = [x.strip() for x in str(taxon).split(taxo_sep)
                           if len(x.strip()) and not x.startswith('x__')]
            split_lens.add(len(taxon_split))
            split_taxa.append(pd.Series(taxon_split))
    # if the parsed and split taxonomies have  very variable number of fields
    # or very long split results  it is terminated here as a taxonomy that
    # does not make sense
    if len(split_lens) > 15 or max(split_lens) > 15:
        return pd.DataFrame([[x] for x in taxa], columns=['not_really_taxon'])

    # build a dataframe from the spit taxonomy
    split_taxa_pd = pd.DataFrame(split_taxa)
    # add a column label
    ALPHA = 'ABCDEFGHIJKLMNOPQRST'
    split_taxa_pd.columns = ['Taxolevel_%s' % (ALPHA[idx])
                             for idx in range(split_taxa_pd.shape[1])]
    # add dummie 0-1 encoding of columns items that are less les 50 per column
    if extended:
        padded_extended_pd = extend_split_taxonomy(split_taxa_pd)
        if padded_extended_pd.shape[0]:
            split_taxa_pd = pd.concat([
                split_taxa_pd, padded_extended_pd], axis=1)
    return split_taxa_pd


def extend_split_taxonomy(split_taxa_pd: pd.DataFrame):
    """

    Parameters
    ----------
    split_taxa_pd

    Returns
    -------

    """
    to_concat = []
    for col in split_taxa_pd.columns.tolist():
        if split_taxa_pd[col].unique().size > 50:
            continue
        split_taxa_dummy = split_taxa_pd[col].str.get_dummies()
        split_taxa_dummy.columns = ['%s__%s' % (x, col) for x
                                    in split_taxa_dummy.columns]
        to_concat.append(split_taxa_dummy)
    if len(to_concat):
        return pd.concat(to_concat, axis=1)
    else:
        return pd.DataFrame()


def parse_split_taxonomy(
        split_taxa_pd: pd.DataFrame) -> dict:
    """

    Parameters
    ----------
    split_taxa_pd : pd.DataFrame

    Returns
    -------
    ranks : dict
    """
    torm = []
    ranks = {}
    # not_collapsable = False
    # parse each column/levels to determine which is
    # (i) to be removed or (ii) labelled with consistent rank (e.g. "k__"))
    for col in split_taxa_pd.columns:
        # if some levels of the split taxonomy are the feature IDs themselves:
        # remove these levels (as it would be a copy of the non-collapsed table)
        col_features = split_taxa_pd[col].tolist()
        if split_taxa_pd.index.tolist() == col_features:
            torm.append(col)
        else:
            rank = get_ranks_from_split_taxonomy(split_taxa_pd, col)
            if rank and rank != col:
                ranks[col] = rank
    # remove columns to be removed
    if torm:
        split_taxa_pd.drop(columns=torm, inplace=True)
    return ranks


def get_ranks_from_split_taxonomy(split_taxa_pd, col):
    """

    Parameters
    ----------
    split_taxa_pd
    col

    Returns
    -------

    """
    rank = [x.split('__')[0] for x in split_taxa_pd[col]
            if str(x) not in ['nan', 'None', 'Unassigned']]
    if len(set(rank)) == 1:
        return rank[0]
    else:
        rank = [x.split('_')[0] for x in split_taxa_pd[col]
                if str(x) not in ['nan', 'None', 'Unassigned']]
        if len(set(rank)) == 1:
            return rank[0]
    return ''


def edit_split_taxonomy(
        ranks: dict,
        split_taxa_pd: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters
    ----------
    ranks
    split_taxa_pd
    """
    if len(ranks) == split_taxa_pd.shape[1]:
        split_taxa_pd = split_taxa_pd.rename(columns=ranks)
    else:
        alpha = 'ABCDEFGHIJKLMNOPQRST'
        cols = [alpha[x] for x in range(split_taxa_pd.shape[1])]
        split_taxa_pd = pd.DataFrame(
            [pd.Series([
                '%s__%s' % (cols[idx], str(x).replace(' ', '_'))
                for idx, x in enumerate(row) if str(x) != 'nan'
            ]) for row in split_taxa_pd.values],
            columns=cols, index=split_taxa_pd.index.tolist()
        )
    return split_taxa_pd


def get_taxonomy_command(self, dat, data):
    """

    Parameters
    ----------
    self
    dat
    data

    Returns
    -------

    """
    cmd = ''
    if data.tax[0] == 'wol':
        cmd = run_taxonomy_wol(self, dat, data)
    elif data.tax[0] in ['amplicon', 'sklearn']:
        if self.config.classifier:
            cmd = run_taxonomy_amplicon(self, dat, data)
        else:
            print('No classifier passed for 16S data\nExiting...')
    else:
        cmd = run_taxonomy_others(self, dat, data)
    return cmd


def run_taxonomy_wol(self, dat, data) -> str:
    """

    Parameters
    ----------
    self
    dat
    data

    Returns
    -------

    """
    cmd = ''
    out_qza, out_tsv, feats = data.tax[1], data.tax[2], data.features
    if self.config.force or to_do(out_qza):
        g2lineage = parse_g2lineage()
        feats_rev = dict((y, x) for x, y in feats.items())
        if to_do(out_tsv):
            readme = out_tsv.replace('taxonomy.tsv', 'readme.txt')
            with open(rep(readme), 'w') as o:
                o.write('Taxonomy file based on WOL reference taxonomy\n')
            with open(rep(out_tsv), 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in data.data[''].ids(axis='observation'):
                    if feats_rev[feat] in g2lineage:
                        o.write('%s\t%s\n' % (feat, g2lineage[feats_rev[feat]]))
                    else:
                        o.write('%s\t%s\n' % (feat, feat.replace('|', '; ')))
        io_update(self, i_f=out_tsv, o_f=out_qza, key=dat)
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_amplicon(self, dat: str, data) -> str:
    """

    Parameters
    ----------
    self
    dat
    data

    Returns
    -------

    """
    cmd = ''
    out_qza, out_tsv = data.tax[1], data.tax[2]
    classifier_qza = get_taxonomy_classifier(self.config.classifier)
    readme = out_tsv.replace('taxonomy.tsv', 'readme.txt')
    with open(rep(readme), 'w') as o:
        o.write('Taxonomy created using "%s"\n' % classifier_qza)
    if not to_do(out_tsv) and to_do(out_qza):
        cmd += run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
        io_update(self, i_f=out_tsv, o_f=out_qza, key=dat)
    else:
        odir_seqs = '%s/sequences/%s' % (self.dir, dat)
        seqs_fasta = '%s/sequences.fasta' % odir_seqs
        seqs_qza = '%s/sequences.qza' % odir_seqs
        if not isdir(rep(odir_seqs)):
            os.makedirs(rep(odir_seqs))
        if self.config.force or to_do(seqs_qza):
            cmd += write_fasta(seqs_fasta, seqs_qza, data.data[''])
            io_update(self, i_f=seqs_fasta, o_f=seqs_qza, key=dat)
        if self.config.force or to_do(out_qza):
            cmd += write_taxonomy_sklearn(classifier_qza, seqs_qza, out_qza)
            cmd += run_export(out_qza, out_tsv, '')
            i_f = [classifier_qza]
            if not to_do(seqs_qza):
                i_f.append(seqs_fasta)
            io_update(self, i_f=i_f, o_f=[out_qza, out_tsv], key=dat)
    return cmd


def run_taxonomy_others(
        self,
        dat,
        data
) -> str:
    """

    Parameters
    ----------
    self
    dat
    data

    Returns
    -------

    """
    cmd = ''
    out_qza, out_tsv = data.tax[1], data.tax[2]
    if self.config.force or to_do(out_qza):
        if to_do(out_tsv):
            readme = out_tsv.replace('taxonomy.tsv', 'readme.txt')
            with open(rep(readme), 'w') as o:
                o.write('Taxonomy created by copying features names\n')
            with open(rep(out_tsv), 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in data.data[''].ids(axis='observation'):
                    o.write('%s\t%s\n' % (feat, feat))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
        io_update(self, i_f=out_tsv, o_f=out_qza, key=dat)
    return cmd


def get_edit_taxonomy_command(self, dat, data):
    """

    Parameters
    ----------
    self
    dat
    data

    Returns
    -------

    """
    cmd = ''
    tax_qza, tax_tsv = data.tax[1], data.tax[2]
    tax_pd, split_taxa_pd = get_tax_tables(rep(tax_tsv))
    taxa = tax_pd['Taxon'].tolist()
    taxa_edit = get_taxa_edit(taxa)
    if split_taxa_pd.shape[1] > 1:
        ranks = parse_split_taxonomy(split_taxa_pd)
        split_taxa_pd = edit_split_taxonomy(ranks, split_taxa_pd)
        split_taxa_pd = split_taxa_pd.fillna('NONE').replace({
            col: {'NONE': '%s__' % col} for col in split_taxa_pd.columns})
        taxa_edit = split_taxa_pd.apply(
            lambda x: ';'.join(x).replace(',', '_'), axis=1).values.tolist()
    if taxa != taxa_edit:
        tax_pd['Taxon'] = taxa_edit
        tax_pd.to_csv(rep(tax_tsv), index=False, sep='\t')
        cmd = run_import(tax_tsv, tax_qza, 'FeatureData[Taxonomy]')
        io_update(self, i_f=tax_qza, o_f=tax_tsv, key=dat)
    return cmd


def get_taxa_edit(taxo):
    taxo_edit = []
    for tax in taxo:
        if str(tax) == 'nan':
            taxo_edit.append(tax)
        elif not tax.strip('_'):
            taxo_edit.append(tax)
        else:
            taxo_edit.append(tax.replace(',', '_'))
    return taxo_edit


def get_gids(features: dict, biom_table: biom.Table):
    """

    Parameters
    ----------
    features : dict
        Mapping genome -> feature
    biom_table : biom.Table
        BIOM table

    Returns
    -------
    gid_features : dict
        Mapping genome -> feature
    """
    data_feats = set(biom_table.ids(axis='observation'))
    gid_features = dict(gid for gid in features.items() if gid[1] in data_feats)
    return gid_features


def get_split_levels(
        collapse_levels: dict, split_taxa_pd: pd.DataFrame) -> tuple:
    """

    Parameters
    ----------
    collapse_levels
    split_taxa_pd

    Returns
    -------

    """
    split_levels = {}
    empties = {}
    for taxo_name, header_index in collapse_levels.items():
        empties[taxo_name] = set()
        split_taxa_index = get_split_taxa_index(split_taxa_pd, header_index)
        if not split_taxa_index:
            continue
        split_levels[str(taxo_name)] = split_taxa_index
        if split_taxa_index > split_taxa_pd.shape[1]:
            continue
        for tdx, tax in enumerate(
                split_taxa_pd.iloc[:, (split_taxa_index-1)].tolist()):
            tax_edit = tax.replace(str(taxo_name), '').strip('__')
            if str(tax) != 'nan' and len(tax_edit) < 3:
                empties[taxo_name].add(
                    ';'.join(split_taxa_pd.iloc[tdx, :split_taxa_index]))
    return split_levels, empties


def get_split_taxa_index(split_taxa_pd: pd.DataFrame, header_index):
    """
    Parameters
    ----------
    split_taxa_pd : pd.DataFrame
    header_index

    Returns
    -------
    split_taxa_index : int
    """
    columns = split_taxa_pd.columns.tolist()
    if isinstance(header_index, int):
        if header_index > len(columns):
            return ''
        split_taxa_index = header_index
    else:
        if header_index not in columns:
            return ''
        split_taxa_index = columns.index(header_index) + 1
    return split_taxa_index


def fix_collapsed_data(
        self,
        dat: str,
        remove_empty: set,
        coll_biom: biom.Table,
        coll_tsv: str,
        coll_qza: str,
        coll_qzv: str,
        coll_meta: str,
):
    """
    Parameters
    ----------
    self
    dat : str
    remove_empty : set
    coll_biom : biom.Table
    coll_tsv : str
    coll_qza : str
    coll_qzv : str
    coll_meta : str

    Returns
    -------
    cmd : str
    """
    cmd = ''
    ids = set(coll_biom.ids(axis='observation'))
    common = remove_empty & ids
    if len(common):
        cmd += '\n# Features removed: %s\n'
        for c in common:
            cmd += '#  - %s\n' % c
        coll_biom.filter(ids_to_keep=list(ids - remove_empty))
        coll_biom.remove_empty()
        coll_pd = coll_biom.to_dataframe(dense=True)
        coll_pd.index.name = 'FeatureID'
        coll_pd.to_csv(rep(coll_tsv), index=True, sep='\t')

        coll_meta_pd = read_meta_pd(rep(coll_meta))
        if coll_meta_pd.index.size != coll_pd.columns.size:
            cmd += '# Feature metadata edited too: %s\n' % coll_meta
            coll_meta_pd = coll_meta_pd.loc[
                coll_meta_pd.sample_name.isin(coll_pd.columns.tolist())]
            coll_meta_pd.to_csv(rep(coll_meta), index=False, sep='\t')
    if to_do(coll_qza):
        cmd += run_import(coll_tsv, coll_qza, 'FeatureTable[Frequency]')
        cmd += run_summary(coll_qza, coll_qzv, coll_meta)
        io_update(self, i_f=coll_tsv, o_f=coll_qza, key=dat)
    return cmd


def find_matching_features(data, subset_regex) -> list:
    """
    Make a feature metadata from the regex to get the names
    of the features to keep.
    """
    feats = data.data[''].ids(axis='observation')
    tax_pd = data.taxa[0]
    to_keep_feats = {}
    for regex in subset_regex:
        to_keep_feats['%s_1' % regex] = pd.Series(
            feats, index=feats).astype(str).str.contains(str(regex))
        if tax_pd.shape[0]:
            tax_pd_regex = tax_pd.copy()
            tax_pd_regex = tax_pd_regex.set_index('Feature ID')
            tax_pd_regex = tax_pd_regex['Taxon'].str.contains(str(regex))
            to_keep_feats['%s_2' % regex] = tax_pd_regex
    to_keep_feats_pd = pd.DataFrame(to_keep_feats)
    to_keep_feats = to_keep_feats_pd.loc[to_keep_feats_pd.any(axis=1)]
    return to_keep_feats.index.tolist()
