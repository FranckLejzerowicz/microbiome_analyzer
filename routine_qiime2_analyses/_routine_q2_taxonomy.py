# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import pandas as pd
import numpy as np
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    read_yaml_file,
    read_meta_pd,
    get_taxonomy_classifier,
    get_job_folder,
    get_analysis_folder,
    parse_g2lineage,
    get_raref_tab_meta_pds,
    simple_chunks
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_barplots,
    write_seqs_fasta,
    write_taxonomy_sklearn,
    write_collapse_taxo,
    run_export,
    run_import
)
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def get_padded_new_rows_list(new_rows, max_new_rows):
    # create a new 'Not available' matrix of shape (n_features x n_fields)
    padded_new_rows_list = []
    for padded_row in new_rows:
        # update row if not of max length
        n_pad = max_new_rows - len(padded_row)
        if n_pad:
            to_pad = ['Not_available']*n_pad
            padded_row = padded_row + to_pad
        padded_new_rows_list.append(padded_row)
    return padded_new_rows_list


def add_alpha_level_label(taxa, padded_new_rows_list, max_new_rows):
    # add a alpha label for rank-identification prupose
    ALPHA = 'ABCDEFGHIJKL'
    cols = ['Node_level_%s' % (ALPHA[idx]) for idx in range(1, max_new_rows + 1)]
    padded_new_rows_pd = pd.DataFrame(
        padded_new_rows_list,
        index = taxa,
        columns = cols
    )
    return padded_new_rows_pd


def extend_split_taxonomy(split_taxa_pd: pd.DataFrame):
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


def get_split_taxonomy(taxa, extended=False, taxo_sep=';'):

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


def get_ranks_from_split_taxonomy(split_taxa_pd, col):
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
    return tax_pd, split_taxa_pd


def edit_split_taxonomy(
        ranks: dict,
        split_taxa_pd: pd.DataFrame) -> pd.DataFrame:
    """

    Parameters
    ----------
    ranks
    split_taxa_pd

    Returns
    -------

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


def get_taxo_levels(taxonomies: dict) -> dict:

    split_taxa_pds = {}
    # for each dataset and its taxonomic file
    for dat, tax_fp in taxonomies.items():
        # default as not necessary to rewrite
        rewrite = False
        # skip if the taxonomic file does not exist
        if not isfile(tax_fp[-1]):
            continue
        tax_fpo = '%s_splitaxa.tsv' % splitext(tax_fp[-1])[0]
        tax_pd, split_taxa_pd = get_tax_tables(tax_fp[-1])
        # write and collect (table, filename) if get_split_taxonomy()
        # found nothing to split (i.e. no taxonomic path)
        if split_taxa_pd.shape[1] == 1:
            split_taxa_pds[dat] = (split_taxa_pd, tax_fpo)
            split_taxa_pd.to_csv(tax_fpo, index=True, sep='\t')
            continue

        ranks = parse_split_taxonomy(split_taxa_pd)
        split_taxa_pd = edit_split_taxonomy(ranks, split_taxa_pd)
        split_taxa_pd.to_csv(tax_fpo, index=True, sep='\t')
        split_taxa_pds[dat] = (split_taxa_pd, tax_fpo)

        if rewrite:
            split_taxa_pd = pd.DataFrame({
                'Feature ID': split_taxa_pd.index.tolist(),
                'Taxon_edit': [';'.join([x for x in row if str(x) != 'nan']) for row in split_taxa_pd.values]
            })
            split_taxa_fpo = '%s_taxSplit.tsv' % splitext(tax_fp[-1])[0]
            tax_extended_pd = tax_pd.merge(split_taxa_pd, on='Feature ID', how='left')
            tax_extended_pd.to_csv(split_taxa_fpo, index=False, sep='\t')
    return split_taxa_pds


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
    if isinstance(header_index, int):
        split_taxa_index = header_index
    else:
        split_taxa_index = split_taxa_pd.columns.tolist().index(
            header_index) + 1
    return split_taxa_index


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
        split_levels[str(taxo_name)] = split_taxa_index
        if split_taxa_index > split_taxa_pd.shape[1]:
            continue
        for tdx, tax in enumerate(
                split_taxa_pd.iloc[:, (split_taxa_index-1)].tolist()):
            tax_edit = tax.replace(str(taxo_name), '').strip('__')
            if str(tax) != 'nan' and len(tax_edit) < 3:
                empties[taxo_name].add(
                    ';'.join(split_taxa_pd.iloc[tdx, :split_taxa_index])
                )
    return split_levels, empties


def make_pies(i_datasets_folder: str, split_taxa_pds: dict,
              datasets_rarefs: dict, datasets_read: dict) -> dict:

    pies_data = {}
    for dat in split_taxa_pds:
        pies_data[dat] = []
        odir = get_analysis_folder(i_datasets_folder, 'nestedness/%s' % dat)
        out_pdf = '%s/pies_%s.pdf' % (odir, dat)
        with PdfPages(out_pdf) as pdf:
            split_taxa, split_taxa_fp = split_taxa_pds[dat]
            ranks = split_taxa.columns.tolist()
            ranks_col = []
            for row in split_taxa.values:
                ranks_col.append(ranks[(len([x for idx, x in enumerate(row)
                                             if str(x).lstrip('%s_' % ranks[idx]).strip('_')])-1)])
            split_taxa['rank'] = ranks_col
            for idx, (tab_, meta) in enumerate(datasets_read[dat]):
                nsams = tab_.shape[1]
                cur_raref = datasets_rarefs[dat][idx]
                tab_sum = tab_.sum(1)
                tab_bool = tab_.astype(bool).sum(1)
                tab = pd.concat([
                    tab_sum,
                    tab_sum / tab_sum.sum(),
                    tab_bool,
                    (tab_bool / nsams) * 100,
                    split_taxa[['rank']]
                ], axis=1)
                tab.columns = ['abundance', 'abundance_percent',
                               'prevalence', 'prevalence_percent', 'rank']
                pies_data_raref = {}
                for min_abundance in [1, 2, 5, 10, 100]:
                    tab_min = tab.loc[tab.abundance >= min_abundance].copy()
                    abundance_bins = [int(x) for x in np.logspace(0, np.log10(tab_min['abundance'].max()+1), num=16)]
                    tab_min['abundance_bin'] = [
                        '%s-%s' % (abundance_bins[x], abundance_bins[x+1]) for x in
                        np.digitize(tab_min['abundance'], bins=abundance_bins[1:], right=True)
                    ]

                    prevalence_bins = [1, 2, 5] + list(range(10, 101, 10))
                    tab_min['prevalence_percent_bin'] = [
                        '%s-%s' % (prevalence_bins[x], prevalence_bins[x+1]) for x in
                        np.digitize(tab_min['prevalence_percent'], bins=prevalence_bins[1:], right=True)
                    ]

                    abundances = tab_min[
                              ['rank', 'abundance_bin', 'abundance']
                          ].groupby(['rank', 'abundance_bin']).count().unstack().fillna(0)
                    abundances.columns = abundances.columns.droplevel()
                    pies_data_raref['abundances'] = abundances

                    prevalences = tab_min[
                        ['rank', 'prevalence_percent_bin', 'prevalence']
                    ].groupby(['rank', 'prevalence_percent_bin']).count().unstack().fillna(0)
                    prevalences.columns = prevalences.columns.droplevel()
                    pies_data_raref['prevalences'] = prevalences

                    tab_gb = pd.concat([
                        tab_min.groupby('rank').count().iloc[:, 0],
                        tab_min[['rank', 'abundance']].groupby('rank').sum(),
                        tab_min[['rank', 'abundance_percent']].groupby('rank').sum(),
                    ], axis=1)
                    tab_gb.columns = ['count', 'abundance_sum', 'abundance_percent_sum']
                    pies_data_raref['tab_gb'] = tab_gb

                    f = plt.figure(figsize=(6, 6))
                    plt.pie([x[0] for x in tab_gb.values],
                            labels=['%s (%s)\n%s reads' % (r, row.iloc[0], row.iloc[1])
                                    for r, row in tab_gb.iterrows()], autopct='%1.2f',
                            startangle=90)
                    plt.title("Number of features assigned per taxon level: %s%s\nmin %s reads, %s features" % (
                        dat, cur_raref, min_abundance, tab_min.shape[0]), size=12)
                    plt.close()
                    pdf.savefig(f, bbox_inches='tight')

                    if abundances.shape[0]:
                        cols = sorted(abundances.columns.tolist(), key=lambda x: int(x.split('-')[0]))
                        ax = abundances.plot(kind='bar', stacked=True)
                        plt.ylabel('Number of features')
                        plt.title("Features per abundance group: %s%s\nmin %s reads, %s features" % (
                            dat, cur_raref, min_abundance, tab_min.shape[0]), size=12)
                        handles, _ = ax.get_legend_handles_labels()
                        plt.legend(handles, cols)
                        pdf.savefig(bbox_inches='tight')
                        plt.close()

                    if prevalences.shape[0]:
                        cols = sorted(prevalences.columns.tolist(), key=lambda x: int(x.split('-')[0]))
                        ax = prevalences.plot(kind='bar', stacked=True)
                        plt.ylabel('Number of features')
                        plt.title("Features per prevalence group: %s%s\nmin %s reads, %s features" % (
                            dat, cur_raref, min_abundance, tab_min.shape[0]), size=12)
                        handles, _ = ax.get_legend_handles_labels()
                        plt.legend(handles, cols)
                        pdf.savefig(bbox_inches='tight')
                        plt.close()
                pies_data[dat].append(pies_data_raref)
    return pies_data


def get_collapse_taxo(
        p_collapse_taxo: str,
        datasets_filt: dict) -> dict:
    """
    Parameters
    ----------
    p_collapse_taxo : str
    datasets_filt : dict

    Returns
    -------
    collapse_taxo : dict
    """
    collapse_taxo = read_yaml_file(p_collapse_taxo)
    collapse_taxo.update(
        dict((datasets_filt[dat], x) for dat, x in
             collapse_taxo.items() if dat in datasets_filt))
    return collapse_taxo


def collapse_paths(dat, tax, tab_fp, meta_fp):
    dat_tax = '%s_tx-%s' % (dat, tax)
    dat_coll = '%s_tx-%s' % (splitext(tab_fp)[0].split('/tab_')[-1], tax)
    coll_tsv = '%s_tx-%s.tsv' % (splitext(tab_fp)[0], tax)
    coll_qza = coll_tsv.replace('.tsv', '.qza')
    coll_meta = '%s_tx-%s.tsv' % (splitext(meta_fp)[0], tax)
    return dat_tax, dat_coll, coll_tsv, coll_qza, coll_meta


def fix_collapsed_data(
        remove_empty: set,
        coll_pd: pd.DataFrame,
        coll_tsv: str,
        coll_qza: str,
        coll_meta: str):
    """
    Parameters
    ----------
    remove_empty : set
    coll_pd : pd.DataFrame
    coll_tsv : str
    coll_qza : str
    coll_meta : str

    Returns
    -------
    cmd : str
    """
    cmd = ''
    if len(remove_empty & set(coll_pd.index)):
        coll_pd = coll_pd.drop(index=list(remove_empty & set(coll_pd.index)))
        coll_pd = coll_pd.loc[:, coll_pd.sum() > 0]
        coll_pd.to_csv(coll_tsv, index=True, sep='\t')

        coll_meta_pd = read_meta_pd(coll_meta)
        if coll_meta_pd.index.size != coll_pd.columns.size:
            coll_meta_pd = coll_meta_pd.loc[
                coll_meta_pd.sample_name.isin(coll_pd.columns.tolist())]
            coll_meta_pd.to_csv(coll_meta, index=False, sep='\t')

        cmd = run_import(coll_tsv, coll_qza, 'FeatureTable[Frequency]')
    return cmd


def run_collapse(i_datasets_folder: str, datasets: dict, datasets_filt: dict,
                 datasets_read: dict, datasets_features: dict,
                 datasets_phylo: dict, split_taxa_pds: dict,
                 taxonomies: dict, p_collapse_taxo: str, datasets_rarefs: dict,
                 datasets_collapsed: dict, datasets_collapsed_map: dict,
                 force: bool, prjct_nm: str, qiime_env: str, chmod: str,
                 noloc: bool, run_params: dict, filt_raref: str,
                 jobs: bool) -> dict:

    collapse_taxo = get_collapse_taxo(p_collapse_taxo, datasets_filt)
    main_written = 0
    collapsed = {}
    datasets_update = {}
    datasets_read_update = {}
    datasets_features_update = {}
    datasets_phylo_update = {}
    stop_for_collapse = False
    job_folder = get_job_folder(i_datasets_folder, 'collapsed_taxo')
    job_folder2 = get_job_folder(i_datasets_folder, 'collapsed_taxo/chunks')
    run_pbs = '%s/3_run_collapsed_taxo_%s%s.sh' % (
        job_folder, prjct_nm, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tab_meta_fps in datasets.items():
            if dat not in collapse_taxo:
                continue
            # get the taxonomic levels
            collapse_levels = collapse_taxo[dat]
            split_taxa_pd, split_taxa_fp = split_taxa_pds[dat]
            split_levels, remove_empties = get_split_levels(collapse_levels,
                                                            split_taxa_pd)
            collapsed[dat] = split_levels

            # files that will be collapsed using qiime2
            tax_qza, tax_fp = taxonomies[dat][1:]

            written = 0
            out_sh = '%s/run_collapsed_taxo_%s_%s%s.sh' % (
                job_folder2, prjct_nm, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, tab_meta_fp in enumerate(tab_meta_fps):
                    tab_fp, meta_fp = tab_meta_fp
                    tab_qza = '%s.qza' % splitext(tab_fp)[0]
                    for tax, level in split_levels.items():
                        coll_paths = collapse_paths(dat, tax, tab_fp, meta_fp)
                        dat_tax = coll_paths[0]
                        dat_coll = coll_paths[1]
                        coll_tsv = coll_paths[2]
                        coll_qza = coll_paths[3]
                        coll_meta = coll_paths[4]
                        if isfile(coll_tsv) and isfile(coll_meta):
                            coll_pd = pd.read_csv(
                                coll_tsv, index_col=0, header=0, sep='\t')
                            if coll_pd.shape[0] < 5:
                                continue
                            cmd = fix_collapsed_data(
                                remove_empties[tax], coll_pd, coll_tsv,
                                coll_qza, coll_meta)
                            if cmd:
                                cur_sh.write('echo "%s"\n' % cmd)
                                cur_sh.write('%s\n\n' % cmd)
                                main_written += 1
                                written += 1
                            datasets_read_update.setdefault(dat_tax, []).append([coll_pd, coll_meta_pd])
                            datasets_collapsed.setdefault(dat, []).append(dat_coll)
                            datasets_collapsed_map[dat_coll] = dat
                            datasets_update.setdefault(dat_tax, []).append([coll_tsv, coll_meta])
                            datasets_rarefs.setdefault(dat_tax, []).append(datasets_rarefs[dat][idx])
                            datasets_phylo_update[dat_tax] = ('', 0)
                        else:
                            written += 1
                            main_written += 1
                            stop_for_collapse = True
                            cmd = write_collapse_taxo(
                                tab_qza, tax_qza, coll_qza, coll_tsv,
                                meta_fp, coll_meta, level, remove_empties[tax])
                            cur_sh.write('echo "%s"\n' % cmd)
                            cur_sh.write('%s\n\n' % cmd)

            if written:
                run_xpbs(
                    out_sh, out_pbs, '%s.cllps.%s%s' % (
                        prjct_nm, dat, filt_raref),
                    qiime_env, run_params["time"], run_params["n_nodes"],
                    run_params["n_procs"], run_params["mem_num"],
                    run_params["mem_dim"], chmod, written, 'single', o,
                    noloc, jobs)

    if main_written:
        print_message('# Collapse features for taxo levels defined in %s' %
                      p_collapse_taxo, 'sh', run_pbs, jobs)

    if stop_for_collapse:
        print('Stopping here as this collapse must be run first for other '
              'analyses to work')
        sys.exit(0)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_phylo.update(datasets_phylo_update)

    return collapsed


def run_taxonomy_others(force: bool, tsv_pd: pd.DataFrame,
                        out_qza: str, out_tsv: str) -> str:
    """
    :param force: Force the re-writing of scripts for all commands.
    :param tsv_pd: Current features table for the current dataset.
    :param out_qza: Taxonomy classification output to generate.
    :param out_tsv: Taxonomy classification output exported.
    """
    cmd = ''
    if force or not isfile(out_qza):
        if not isfile(out_tsv):
            with open(out_tsv, 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in tsv_pd.index:
                    o.write('%s\t%s\n' % (feat, feat))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_wol(force: bool, tsv_pd: pd.DataFrame, out_qza: str,
                     out_tsv: str, cur_datasets_features: dict) -> str:
    """
    :param force: Force the re-writing of scripts for all commands.
    :param tsv_pd: Current features table for the current dataset.
    :param out_qza: Taxonomy classification output to generate.
    :param out_tsv: Taxonomy classification output exported.
    :param cur_datasets_features: gotu -> features name containing gotu
    """
    cmd = ''
    if force or not isfile(out_qza):
        g2lineage = parse_g2lineage()
        rev_cur_datasets_features = dict(
            (y, x) for x, y in cur_datasets_features.items())
        if not isfile(out_tsv):
            with open(out_tsv, 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in tsv_pd.index:
                    if rev_cur_datasets_features[feat] in g2lineage:
                        o.write('%s\t%s\n' % (
                            feat, g2lineage[rev_cur_datasets_features[feat]]))
                    else:
                        o.write('%s\t%s\n' % (feat, feat.replace('|', '; ')))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_amplicon(
        dat: str, i_datasets_folder: str, force: bool, tsv_pd: pd.DataFrame,
        out_qza: str, out_tsv: str, i_classifier: str) -> str:
    """
    :param dat: Current dataset.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param force: Force the re-writing of scripts for all commands.
    :param tsv_pd: Current features table for the current dataset.
    :param out_qza: Taxonomy classification output to generate.
    :param out_tsv: Taxonomy classification output exported.
    :param i_classifier: Path to the taxonomic classifier.
    """
    cmd = ''
    if isfile(out_tsv) and not isfile(out_qza):
        cmd += run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    else:
        ref_classifier_qza = get_taxonomy_classifier(i_classifier)
        odir_seqs = get_analysis_folder(i_datasets_folder, 'seqs/%s' % dat)
        out_fp_seqs_rad = '%s/seq_%s' % (odir_seqs, dat)
        out_fp_seqs_fasta = '%s.fasta' % out_fp_seqs_rad
        out_fp_seqs_qza = '%s.qza' % out_fp_seqs_rad
        if force or not isfile(out_fp_seqs_qza):
            cmd += write_seqs_fasta(out_fp_seqs_fasta, out_fp_seqs_qza, tsv_pd)
        if force or not isfile(out_qza):
            cmd += write_taxonomy_sklearn(out_qza, out_fp_seqs_qza,
                                          ref_classifier_qza)
            cmd += run_export(out_qza, out_tsv, '')
    return cmd


def get_gid_features(features: dict, data_pd: pd.DataFrame):
    """

    Parameters
    ----------
    features : dict
        Mapping genome -> feature
    data_pd : pd.DataFrame
        Feature table

    Returns
    -------
    gid_features : dict
        Mapping genome -> feature
    """
    data_feats = set(data_pd.index)
    gid_features = dict(gid for gid in features.items() if gid[1] in data_feats)
    return gid_features


def get_edit_taxonomy_command(data):
    cmd = ''
    out_pd = pd.read_csv(data.tax[2], dtype=str, sep='\t')
    taxo = out_pd['Taxon'].tolist()
    taxo_edit = get_taxa_edit(taxo)
    if taxo != taxo_edit:
        out_pd['Taxon'] = taxo_edit
        out_pd.to_csv(data.tax[2], index=False, sep='\t')
        cmd = run_import(
            data.tax[2], data.tax[1], 'FeatureData[Taxonomy]')
    return cmd


def get_taxonomy_command(dat, config, data):
    cmd = ''
    if data.tax:
        if data.tax[0] == 'wol':
            cmd = run_taxonomy_wol(
                config.force, data.data[0], data.tax[1],
                data.tax[2], data.features)
        elif data.tax[0] == 'amplicon':
            if config.i_classifier:
                cmd = run_taxonomy_amplicon(
                    dat, config.i_datasets_folder, config.force, data.data[0],
                    data.tax[1], data.tax[2], config.i_classifier)
            else:
                print('No classifier passed for 16S data\nExiting...')
    # elif not data.data[0].index.is_numeric():
    # and also add classyfire
    else:
        cmd = run_taxonomy_others(
            config.force, data.data[0], data.tax[1], data.tax[2])
    return cmd


def run_taxonomy(
        method: str, i_datasets_folder: str, datasets: dict,
        datasets_read: dict, datasets_phylo: dict, datasets_features: dict,
        datasets_filt_map: dict, i_classifier: str, taxonomies: dict,
        force: bool, prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
        run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
    """

    Parameters
    ----------
    method
    i_datasets_folder : str
        Path to the folder containing the data/metadata subfolders.
    datasets : dict
        Mappring dataset name -> [data file path, metadata file path].
    datasets_read : dict
        Mapping dataset name -> [data table, metadata table]
    datasets_phylo : dict
        To be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    datasets_features : dict
        Mapping dataset name -> list of features names in
                                the dataset tsv / biom file.
    datasets_filt_map : dict
    i_classifier : str
        Path to the taxonomic classifier.
    taxonomies : dict
        Mapping Dataset name -> [method, assignment qza]
    force : bool
        Force the re-writing of scripts for all commands.
    prjct_nm : str
        Short nick name for your project.
    qiime_env : str
        Name of your qiime2 conda environment (e.g. qiime2-2019.10).
    chmod : str
        Whether to change permission of output files (default: 744).
    noloc : str
    run_params : dict
    filt_raref : str
    jobs : bool
    chunkit : int

    Returns
    -------

    """
    job_folder = get_job_folder(i_datasets_folder, 'taxonomy')
    job_folder2 = get_job_folder(i_datasets_folder, 'taxonomy/chunks')
    amplicon_datasets = [dat for dat, (tree, correction) in
                         datasets_phylo.items() if tree == 'amplicon']
    wol_datasets = [dat for dat, (tree, correction)
                    in datasets_phylo.items() if tree == 'wol']

    main_written = 0
    to_chunk = []
    run_pbs = '%s/1_run_taxonomy_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds_ in datasets_read.items():
            out_sh = '%s/run_taxonomy_%s_%s%s.sh' % (
                job_folder2, prjct_nm, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            if dat in datasets_filt_map:
                taxonomies[dat] = taxonomies[datasets_filt_map[dat]]
                continue
            written = 0
            with open(out_sh, 'w') as cur_sh:
                for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
                    if idx:
                        continue
                    tsv, meta = datasets[dat][idx]
                    if not isinstance(tsv_meta_pds[0], pd.DataFrame) and \
                            tsv_meta_pds[0] == 'raref':
                        if not isfile(tsv):
                            print('Must have run rarefaction to use it '
                                  'further...\nExiting')
                            sys.exit(0)
                        tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                        datasets_read[dat][idx] = [tsv_pd, meta_pd]
                    else:
                        tsv_pd, meta_pd = tsv_meta_pds

                    odir = get_analysis_folder(i_datasets_folder,
                                               'taxonomy/%s' % dat)
                    out_rad = '%s/tax_%s' % (odir, dat)

                    if dat in amplicon_datasets:
                        out_qza = '%s_%s.qza' % (out_rad, method)
                        out_tsv = '%s.tsv' % splitext(out_qza)[0]
                        taxonomies[dat] = [method, out_qza, out_tsv]
                        if not i_classifier:
                            print('No classifier passed for 16S '
                                  'data\nExiting...')
                            continue
                        cmd = run_taxonomy_amplicon(
                            dat, i_datasets_folder, force, tsv_pd,
                            out_qza, out_tsv, i_classifier)
                    else:
                        out_qza = '%s.qza' % out_rad
                        out_tsv = '%s.tsv' % out_rad
                        if dat in wol_datasets:
                            cur_datasets_features = datasets_features[dat]
                            taxonomies[dat] = ['wol', out_qza, out_tsv]
                            cmd = run_taxonomy_wol(
                                force, tsv_pd, out_qza,
                                out_tsv, cur_datasets_features)
                        else:
                            if len([x for x in tsv_pd.index if str(
                                    x).isdigit()]) == tsv_pd.shape[0]:
                                continue
                            taxonomies[dat] = ['feat', out_qza, out_tsv]
                            cmd = run_taxonomy_others(
                                force, tsv_pd, out_qza, out_tsv)
                    if cmd:
                        cur_sh.write('echo "%s"\n' % cmd)
                        cur_sh.write('%s\n\n' % cmd)
                        main_written += 1
                        written += 1
            if written:
                to_chunk.append(out_sh)
                if not chunkit:
                    run_xpbs(
                        out_sh, out_pbs,
                        '%s.tx.sklrn.%s%s' % (prjct_nm, dat, filt_raref),
                        qiime_env, run_params["time"], run_params["n_nodes"],
                        run_params["n_procs"], run_params["mem_num"],
                        run_params["mem_dim"], chmod, written, 'single', o,
                        noloc, jobs)

    if to_chunk and chunkit:
        simple_chunks(run_pbs, job_folder2, to_chunk, 'taxonomy',
                      prjct_nm, run_params["time"], run_params["n_nodes"],
                      run_params["n_procs"], run_params["mem_num"],
                      run_params["mem_dim"], qiime_env, chmod, noloc, jobs,
                      chunkit, None)

    if main_written:
        print_message('# Classify features using classify-sklearn', 'sh',
                      run_pbs, jobs)


def run_barplot(i_datasets_folder: str, datasets: dict, taxonomies: dict,
                force: bool, prjct_nm: str, qiime_env: str,
                chmod: str, noloc: bool, run_params: dict,
                filt_raref: str, jobs: bool, chunkit: int) -> None:
    """Visualize taxonomy with an interactive bar plot.

    Parameters
    ----------
    i_datasets_folder : str
        Path to the folder containing the data/metadata subfolders
    datasets : dict
        Mappig dataset name -> [tsv file path, metadata file path]
    taxonomies : dict
        Mappig dataset name -> [classification_method, tax_qza]
    force : bool
        Force the re-writing of scripts for all commands
    prjct_nm : str
        Short nick name for your project
    qiime_env : str
        Mame of a qiime2 conda environment
    chmod : str
        Whether to change permission of output files (defalt: 744)
    noloc : bool
    run_params : dict
    filt_raref : str
    jobs : bool
    chunkit : int

    Returns
    -------

    """
    job_folder = get_job_folder(i_datasets_folder, 'barplot')
    job_folder2 = get_job_folder(i_datasets_folder, 'barplot/chunks')

    written = 0
    to_chunk = []
    run_pbs = '%s/1_run_barplot_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds_ in datasets.items():
            out_sh = '%s/run_barplot_%s_%s%s.sh' % (job_folder2, prjct_nm,
                                                    dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for tsv_meta_pds in tsv_meta_pds_:
                    tsv, meta = tsv_meta_pds
                    if dat not in taxonomies:
                        continue
                    method, tax_qza, tax_tsv = taxonomies[dat]
                    if not method:
                        method = 'taxofromfile'
                    qza = '%s.qza' % splitext(tsv)[0]
                    odir = get_analysis_folder(i_datasets_folder,
                                               'barplot/%s' % dat)
                    out_qzv = '%s/bar_%s_%s.qzv' % (odir, dat, method)
                    if force or not isfile(out_qzv):
                        write_barplots(out_qzv, qza, meta, tax_qza, cur_sh)
                        written += 1
            to_chunk.append(out_sh)
            if not chunkit:
                run_xpbs(
                    out_sh, out_pbs, '%s.brplt.%s%s' % (
                        prjct_nm, dat, filt_raref),
                    qiime_env, run_params["time"], run_params["n_nodes"],
                    run_params["n_procs"], run_params["mem_num"],
                    run_params["mem_dim"], chmod, written, 'single', o,
                    noloc, jobs)

    if to_chunk and chunkit:
        simple_chunks(run_pbs, job_folder2, to_chunk, 'barplot',
                      prjct_nm, run_params["time"], run_params["n_nodes"],
                      run_params["n_procs"], run_params["mem_num"],
                      run_params["mem_dim"], qiime_env, chmod, noloc, jobs,
                      chunkit, None)

    if written:
        print_message('# Make sample compositions barplots',
                      'sh', run_pbs, jobs)


def get_raw_of_filtered_dataset(
        dat: str,
        datasets_filt_map: dict) -> str:
    """
    Parameters
    ----------
    dat : str
        Dataset name
    datasets_filt_map : dict
        Mapping filtered dataset name -> raw dataset name

    Returns
    -------
    dat_tax : str
        Raw dataset name
    """
    if dat in datasets_filt_map:
        dat_raw = datasets_filt_map[dat]
    else:
        dat_raw = dat
    return dat_raw


def get_precomputed_taxonomies(
        i_datasets_folder: str,
        datasets: dict,
        datasets_filt_map: dict,
        taxonomies: dict,
        method: str) -> None:
    """Update taxonomies dict with file found existing

    Parameters
    ----------
    i_datasets_folder
        Path to the folder containing the data/metadata subfolders
    datasets : dict
        Mapping dataset name -> [tsv file path, metadata file path]
    datasets_filt_map : dict
    taxonomies : dict
        Mapping dataset name -> [classification_method, tax_qza]
    method : str
    """
    for dat in datasets:
        dat_raw = get_raw_of_filtered_dataset(dat, datasets_filt_map)
        analysis_folder = get_analysis_folder(
            i_datasets_folder, 'taxonomy/%s' % dat_raw)

        tax_qza = '%s/tax_%s_%s.qza' % (analysis_folder, dat_raw, method)
        tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
        if isfile(tax_tsv):
            taxonomies[dat] = ['', tax_qza, tax_tsv]

        tax_qza = '%s/tax_%s.qza' % (analysis_folder, dat_raw)
        tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
        if isfile(tax_tsv):
            taxonomies[dat] = ['', tax_qza, tax_tsv]


def create_songbird_feature_metadata(i_datasets_folder: str,
                                     taxonomies: dict,
                                     q2_pd: pd.DataFrame):

    q2_pd = q2_pd.loc[(q2_pd.pair == 'no_pair') & (q2_pd.Pseudo_Q_squared > 0)]
    for dat in taxonomies.keys():
        dat_q2_pd = q2_pd.loc[q2_pd.dat.str.contains(dat)]
        dat_sbs = []
        for (pair, dat, dataset_filter, subset, model, songbird_filter,
             parameters, baseline, differentials, Pseudo_Q_squared) in dat_q2_pd.values:
            sb_pd = pd.read_table(differentials, index_col=0).iloc[1:]
            sb_pd.columns = ['%s__%s__%s__%s__%s__%s__%s (Q2=%s): %s' % (
                dat, dataset_filter, subset, model, songbird_filter,
                parameters, baseline, Pseudo_Q_squared, x
            ) for x in sb_pd.columns]
            dat_sbs.append(sb_pd)
        if len(dat_sbs):
            dat_sbs_pd = pd.concat(dat_sbs, axis=1, sort=False)
            odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % dat)
            fpo_tsv = '%s/sb_%s.tsv' % (odir, dat)
            fpo_qza = '%s/sb_%s.qza' % (odir, dat)
            dat_sbs_pd.reset_index().rename(
                columns={dat_sbs_pd.reset_index().columns.tolist()[0]: 'Feature ID'}
            ).to_csv(fpo_tsv, index=True, sep='\t')
            run_import(fpo_tsv, fpo_qza, 'FeatureData[Differential]')


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


def edit_taxonomies(
        i_datasets_folder: str, taxonomies: dict, force: bool, prjct_nm: str,
        qiime_env: str, chmod: str, noloc: bool, run_params: dict,
        filt_raref: str, jobs: bool, chunkit: int):

    job_folder = get_job_folder(i_datasets_folder, 'taxonomy')
    job_folder2 = get_job_folder(i_datasets_folder, 'taxonomy/chunks')

    main_written = 0
    to_chunk = []
    run_pbs = '%s/1_run_taxonomy_edit_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, (_, qza, tsv) in taxonomies.items():
            if not isfile(tsv):
                continue
            written = 0
            out_pd = pd.read_csv(tsv, dtype=str, sep='\t')
            taxo = out_pd['Taxon'].tolist()
            taxo_edit = get_taxa_edit(taxo)
            if taxo != taxo_edit:
                out_pd['Taxon'] = taxo_edit
                out_pd.to_csv(tsv, index=False, sep='\t')
                cmd = run_import(tsv, qza, 'FeatureData[Taxonomy]')

                out_sh = '%s/run_taxonomy_edit_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
                out_pbs = '%s.pbs' % splitext(out_sh)[0]
                with open(out_sh, 'w') as cur_sh:
                    cur_sh.write('echo "%s"\n' % cmd)
                    cur_sh.write('%s\n\n' % cmd)
                    main_written += 1
                    written += 1
                if written:
                    to_chunk.append(out_sh)
                if not chunkit:
                    run_xpbs(out_sh, out_pbs, '%s.tx.dt.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                             run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                             run_params["mem_num"], run_params["mem_dim"],
                             chmod, written, 'single', o, noloc, jobs)

    if to_chunk and chunkit:
        simple_chunks(run_pbs, job_folder2, to_chunk, 'taxonomy_edit',
                      prjct_nm, run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                      run_params["mem_num"], run_params["mem_dim"],
                      qiime_env, chmod, noloc, jobs, chunkit, None)

    if main_written:
        print_message('# Edit features taxonomy to not contain "," characters', 'sh', run_pbs, jobs)

