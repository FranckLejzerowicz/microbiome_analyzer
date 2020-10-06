# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import re
import pandas as pd
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_taxonomy_classifier,
    get_job_folder,
    get_analysis_folder,
    parse_g2lineage,
    get_raref_tab_meta_pds,
    get_collapse_taxo,
    get_songbird_outputs
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_barplots,
    write_seqs_fasta,
    write_taxonomy_sklearn,
    write_collapse_taxo,
    run_export,
    run_import
)


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
        split_taxa_dummy.columns = ['%s__%s' % (x, col) for x in split_taxa_dummy.columns]
        to_concat.append(split_taxa_dummy)
    if len(to_concat):
        return pd.concat(to_concat, axis=1)
    else:
        return pd.DataFrame()


def get_split_taxonomy(taxa, extended=False, taxo_sep=';'):

    # get the taxon name split per "taxon" level
    # split_chars = taxo_sep
    # if len([1 for x in taxa if '|' in x]) > (0.5 * len(taxa)):
    #     split_chars += "|\|"
    # if len([1 for x in taxa if '.' in x]) > (0.5 * len(taxa)):
    #     split_chars += "|\."

    split_lens = set()
    split_taxa = []
    for taxon in taxa:
        taxon_split = [x.strip() for x in str(taxon).split(taxo_sep) if len(x.strip()) and not x.startswith('x__')]
        # taxon_split = [x.strip() for x in re.split(split_chars, str(taxon)) if len(x.strip()) and not x.startswith('x__')]
        split_lens.add(len(taxon_split))
        split_taxa.append(taxon_split)

    if len(split_lens) > 15 or max(split_lens) > 15:
        return pd.DataFrame([[x] for x in taxa], columns=['not_really_taxon'])

    split_taxa_pd = pd.DataFrame(split_taxa)
    ALPHA = 'ABCDEFGHIJKLMNOPQRST'
    split_taxa_pd_cols = ['Taxolevel_%s' % (ALPHA[idx]) for idx in range(split_taxa_pd.shape[1])]
    split_taxa_pd.columns = split_taxa_pd_cols
    # get the max number of fields
    # max_new_rows = max([len(new_row) for new_row in new_rows])
    # padded_new_rows_list = get_padded_new_rows_list(new_rows, max_new_rows)
    # padded_new_rows_pd = add_alpha_level_label(taxa, padded_new_rows_list, max_new_rows)
    if extended:
        padded_extended_pd = extend_split_taxonomy(split_taxa_pd)
        if padded_extended_pd.shape[0]:
            split_taxa_pd = pd.concat([split_taxa_pd, padded_extended_pd], axis=1)
        # split_taxa_pd = split_taxa_pd.reset_index().rename(columns={'index': 'Taxon'})
    return split_taxa_pd


def get_taxo_levels(taxonomies: dict) -> dict:

    split_taxa_pds = {}
    for dat, tax_fp in taxonomies.items():
        rewrite = False
        tax_pd = pd.read_csv(tax_fp[-1], header=0, sep='\t', dtype=str)
        tax_pd.rename(columns={tax_pd.columns[0]: 'Feature ID'}, inplace=True)
        features = tax_pd['Feature ID'].tolist()
        split_taxa_pd = get_split_taxonomy(tax_pd.Taxon.tolist())
        if split_taxa_pd.shape[1] == 1:
            split_taxa_pds[dat] = split_taxa_pd
            continue

        torm = []
        for col in split_taxa_pd.columns:
            col_features = split_taxa_pd[col].tolist()
            if features == col_features:
                rewrite = True
                torm.append(col)

        if rewrite:
            split_taxa_pd = split_taxa_pd.drop(columns=torm)

        ranks = {}
        not_collapsable = False
        for col in split_taxa_pd.columns:
            rank = [x.split('_')[0] for x in split_taxa_pd[col] if str(x) not in ['nan', 'None']]
            if len(rank) == split_taxa_pd.shape[0]:
                if len(set(rank)) == 1:
                    ranks[col] = list(rank)[0]
            else:
                not_collapsable = True

        if not_collapsable:
            split_taxa_pds[dat] = split_taxa_pd
            continue

        if len(ranks) == split_taxa_pd.shape[1]:
            split_taxa_pd = split_taxa_pd.rename(columns=ranks)
        else:
            rewrite = True
            alpha = 'ABCDEFGHIJKLMNOPQRST'
            cols = [alpha[x] for x in range(split_taxa_pd.shape[1])]
            split_taxa_pd = pd.DataFrame(
                [['%s__%s' % (cols[idx], str(x).replace(' ', '_')) for idx, x in enumerate(row)]
                  for row in split_taxa_pd.values],
                columns=cols
            )
        split_taxa_pds[dat] = split_taxa_pd
        if rewrite:
            split_taxa_pd = pd.DataFrame({
                'Feature ID': features,
                'Taxon_edit': [';'.join([x for x in row if str(x)]) for row in split_taxa_pd.values]
            })
            split_taxa_fpo = '%s_taxSplit.tsv' % splitext(tax_fp[-1])[0]
            tax_extended_pd = tax_pd.merge(split_taxa_pd, on='Feature ID', how='left')
            tax_extended_pd.to_csv(split_taxa_fpo, index=False, sep='\t')

    return split_taxa_pds


def get_split_levels(dat, collapse_taxo: dict, split_taxa_pds: dict):
    split_levels = {}
    taxo = collapse_taxo[dat]
    split_taxa_pd = split_taxa_pds[dat]
    # taxo levels are the header of split_taxa_pd
    for taxo_name, taxo_header in taxo.items():
        if isinstance(taxo_header, int):
            split_levels[taxo_name] = taxo_header
        else:
            split_levels[taxo_name] = split_taxa_pd.columns.tolist().index(taxo_header)+1
    return split_levels


def run_collapse(i_datasets_folder: str, datasets: dict, datasets_read: dict,
                 datasets_features: dict, datasets_phylo: dict, split_taxa_pds: dict,
                 taxonomies: dict, p_collapse_taxo: str, datasets_rarefs: dict,
                 datasets_collapsed: dict, datasets_collapsed_map: dict, force: bool,
                 prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
                 run_params: dict, filt_raref: str) -> None:

    collapse_taxo = get_collapse_taxo(p_collapse_taxo)

    main_written = 0
    datasets_update = {}
    datasets_read_update = {}
    datasets_features_update = {}
    datasets_phylo_update = {}
    job_folder = get_job_folder(i_datasets_folder, 'collapsed_taxo')
    job_folder2 = get_job_folder(i_datasets_folder, 'collapsed_taxo/chunks')
    run_pbs = '%s/3_run_collapsed_taxo%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tab_meta_fps in datasets.items():
            if dat not in collapse_taxo:
                continue
            written = 0
            out_sh = '%s/run_collapsed_taxo_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                split_levels = get_split_levels(dat, collapse_taxo, split_taxa_pds)
                tax_qza, tax_fp = taxonomies[dat][1:]
                for idx, tab_meta_fp in enumerate(tab_meta_fps):
                    tab_fp, meta_fp = tab_meta_fp
                    tab_qza = '%s.qza' % splitext(tab_fp)[0]
                    for tax, level in split_levels.items():
                        dat_tax = '%s_tx-%s' % (dat, tax)
                        dat_collapsed = '%s_tx-%s' % (splitext(tab_fp)[0].split('/tab_')[-1], tax)
                        datasets_collapsed.setdefault(dat, []).append(dat_collapsed)
                        datasets_collapsed_map[dat_collapsed] = dat

                        collapsed_tsv = '%s_tx-%s.tsv' % (splitext(tab_fp)[0], tax)
                        collapsed_qza = collapsed_tsv.replace('.tsv', '.qza')
                        collapsed_meta = '%s_tx-%s.tsv' % (splitext(meta_fp)[0], tax)

                        datasets_update.setdefault(dat_tax, []).append([collapsed_tsv, collapsed_meta])
                        datasets_rarefs.setdefault(dat_tax, []).append(datasets_rarefs[dat][idx])
                        datasets_phylo_update[dat_tax] = ('', 0)

                        if isfile(collapsed_tsv) and isfile(collapsed_meta):
                            collapsed_pd = pd.read_csv(collapsed_tsv, index_col=0, header=0, sep='\t')
                            with open(collapsed_tsv) as f:
                                for line in f:
                                    break
                            collapsed_meta_pd = pd.read_csv(
                                collapsed_meta, header=0, sep='\t',
                                dtype={line.split('\t')[0]: str}
                            )
                            datasets_read_update.setdefault(dat_tax, []).append([collapsed_pd, collapsed_meta_pd])
                            continue
                        else:
                            written += 1
                            main_written += 1
                            write_collapse_taxo(tab_qza, tax_qza, collapsed_qza, collapsed_tsv,
                                                meta_fp, collapsed_meta, level, cur_sh)
                            # meta_pd = meta_pd.set_index('sample_name')
                            # collapsed_meta_pd = meta_pd.loc[collapsed_pd.columns.tolist()].copy()
                            # collapsed_pd.reset_index().to_csv(collapsed_tsv, index=False, sep='\t')
                            # collapsed_meta_pd.reset_index().to_csv(collapsed_meta, index=False, sep='\t')
                            #
                            # datasets_update[dat_collapsed] = [[collapsed_tsv, collapsed_meta]]
                            #
                            # # EXPORT AND READ TO ADD HERE:
                            # datasets_read_update[dat_collapsed] = [[collapsed_pd, collapsed_meta_pd.reset_index()]]
            if written:
                run_xpbs(out_sh, out_pbs, '%s.cllps.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                         run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                         run_params["mem_num"], run_params["mem_dim"],
                         chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Collapse features for taxo levels defined in %s' % p_collapse_taxo, 'sh', run_pbs)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_features.update(datasets_features_update)
    datasets_phylo.update(datasets_phylo_update)


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
        rev_cur_datasets_features = dict((y, x) for x, y in cur_datasets_features.items())
        if not isfile(out_tsv):
            with open(out_tsv, 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in tsv_pd.index:
                    if rev_cur_datasets_features[feat] in g2lineage:
                        o.write('%s\t%s\n' % (feat, g2lineage[rev_cur_datasets_features[feat]]))
                    else:
                        o.write('%s\t%s\n' % (feat, feat.replace('|', '; ')))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_amplicon(dat: str, i_datasets_folder: str, force: bool, tsv_pd: pd.DataFrame,
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
            cmd += write_taxonomy_sklearn(out_qza, out_fp_seqs_qza, ref_classifier_qza)
            cmd += run_export(out_qza, out_tsv, '')
    return cmd


def run_taxonomy(method: str, i_datasets_folder: str, datasets: dict, datasets_read: dict,
                 datasets_phylo: dict, datasets_features: dict, datasets_filt_map: dict,
                 i_classifier: str, taxonomies: dict, force: bool, prjct_nm: str, qiime_env: str,
                 chmod: str, noloc: bool, run_params: dict, filt_raref: str) -> None:
    """
    classify-sklearn: Pre-fitted sklearn-based taxonomy classifier

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param i_classifier: Path to the taxonomic classifier.
    :param taxonomies: dataset -> [method, assignment qza]
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'taxonomy')
    job_folder2 = get_job_folder(i_datasets_folder, 'taxonomy/chunks')
    amplicon_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'amplicon']
    wol_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'wol']

    main_written = 0
    run_pbs = '%s/1_run_taxonomy%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds_ in datasets_read.items():
            out_sh = '%s/run_taxonomy_%s%s.sh' % (job_folder2, dat, filt_raref)
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
                    if not isinstance(tsv_meta_pds[0], pd.DataFrame) and tsv_meta_pds[0] == 'raref':
                        if not isfile(tsv):
                            print('Must have run rarefaction to use it further...\nExiting')
                            sys.exit(0)
                        tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                        datasets_read[dat][idx] = [tsv_pd, meta_pd]
                    else:
                        tsv_pd, meta_pd = tsv_meta_pds

                    odir = get_analysis_folder(i_datasets_folder, 'taxonomy/%s' % dat)
                    out_rad = '%s/tax_%s' % (odir, dat)

                    if dat in amplicon_datasets:
                        out_qza = '%s_%s.qza' % (out_rad, method)
                        out_tsv = '%s.tsv' % splitext(out_qza)[0]
                        # if dat in taxonomies:
                        #     if out_tsv == [dat][-1]:
                        #         continue
                        #     tax_tsv_ref = taxonomies[dat][-1]
                        #     tax_tsv_ref_pd = pd.read_csv(tax_tsv_ref, header=0, index_col=0, sep='\t')
                        #     print(tax_tsv_ref_pd.shape)
                        #     print(tax_tsv_ref_pd.iloc[:3,:3])
                        #     print(tsv_pd.shape)
                        #     print(tsv_pd.iloc[:3,:3])
                        #     tax_pd = tax_tsv_ref_pd.loc[tsv_pd.index]
                        #     print(tax_pd.shape)
                        #     print(tax_pd.iloc[:3,:3])
                        #     print(tax_pvdv)
                        # else:
                        taxonomies[dat] = [method, out_qza, out_tsv]
                        if not i_classifier:
                            print('No classifier passed for 16S data\nExiting...')
                            continue
                        cmd = run_taxonomy_amplicon(dat, i_datasets_folder, force, tsv_pd,
                                                    out_qza, out_tsv, i_classifier)
                    else:
                        out_qza = '%s.qza' % out_rad
                        out_tsv = '%s.tsv' % out_rad
                        if dat in wol_datasets:
                            cur_datasets_features = datasets_features[dat]
                            taxonomies[dat] = ['wol', out_qza, out_tsv]
                            cmd = run_taxonomy_wol(force, tsv_pd, out_qza, out_tsv,
                                                   cur_datasets_features)
                        else:
                            if len([x for x in tsv_pd.index if str(x).isdigit()]) == tsv_pd.shape[0]:
                                continue
                            taxonomies[dat] = ['feat', out_qza, out_tsv]
                            cmd = run_taxonomy_others(force, tsv_pd, out_qza, out_tsv)
                    if cmd:
                        cur_sh.write('echo "%s"\n' % cmd)
                        cur_sh.write('%s\n\n' % cmd)
                        main_written += 1
                        written += 1
            if written:
                run_xpbs(out_sh, out_pbs, '%s.tx.sklrn.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                         run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                         run_params["mem_num"], run_params["mem_dim"],
                         chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Classify features using classify-sklearn', 'sh', run_pbs)


def run_barplot(i_datasets_folder: str, datasets: dict, taxonomies: dict,
                force: bool, prjct_nm: str, qiime_env: str,
                chmod: str, noloc: bool, run_params: dict, filt_raref: str) -> None:
    """
    barplot: Visualize taxonomy with an interactive bar plot

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param taxonomies: dataset -> [classification_method, tax_qza]
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'barplot')
    job_folder2 = get_job_folder(i_datasets_folder, 'barplot/chunks')

    written = 0
    run_pbs = '%s/1_run_barplot%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds_ in datasets.items():
            out_sh = '%s/run_barplot_%s%s.sh' % (job_folder2, dat, filt_raref)
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
                    odir = get_analysis_folder(i_datasets_folder, 'barplot/%s' % dat)
                    out_qzv = '%s/bar_%s_%s.qzv' % (odir, dat, method)
                    if force or not isfile(out_qzv):
                        write_barplots(out_qzv, qza, meta, tax_qza, cur_sh)
                        written += 1
            run_xpbs(out_sh, out_pbs, '%s.brplt.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc)
    if written:
        print_message('# Make sample compositions barplots', 'sh', run_pbs)


def get_precomputed_taxonomies(i_datasets_folder: str, datasets: dict,
                               datasets_filt_map: dict, taxonomies: dict,
                               method: str) -> None:
    """
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param taxonomies: dataset -> [classification_method, tax_qza]
    """
    for dat in datasets:

        if dat in datasets_filt_map:
            dat_tax = datasets_filt_map[dat]
        else:
            dat_tax = dat

        analysis_folder = get_analysis_folder(i_datasets_folder, 'taxonomy/%s' % dat_tax)

        tax_qza = '%s/tax_%s_%s.qza' % (analysis_folder, dat_tax, method)
        tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
        if isfile(tax_tsv):
            taxonomies[dat] = ['', tax_qza, tax_tsv]

        tax_qza = '%s/tax_%s.qza' % (analysis_folder, dat_tax)
        tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
        if isfile(tax_tsv):
            taxonomies[dat] = ['', tax_qza, tax_tsv]


def create_feature_metadata(i_datasets_folder: str, taxonomies: dict, q2_pd: pd.DataFrame):
    q2_pd = q2_pd.loc[(q2_pd.pair == 'no_pair') & (q2_pd.Pseudo_Q_squared > 0)]
    for dat, tax_fps in taxonomies.items():
        tax_tsv = tax_fps[-1]
        tax_pd = pd.read_table(tax_tsv, index_col=0)

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
            tax_sbs_pd = pd.concat([tax_pd, dat_sbs_pd], axis=1, sort=False)
            print(tax_sbs_pd.columns)
            odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % dat)
            fpo = '%s/sb_%s.tsv' % (odir, dat)
            tax_sbs_pd.reset_index().rename(
                columns={tax_sbs_pd.reset_index().columns.tolist()[0]: 'Feature ID'}
            ).to_csv(fpo, index=True, sep='\t')
