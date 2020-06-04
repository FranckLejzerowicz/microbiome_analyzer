# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, sys
import pandas as pd
from os.path import basename, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics,
    get_job_folder,
    get_analysis_folder,
    get_subsets,
    get_procrustes_dicts,
    get_raref_tab_meta_pds
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_diversity_beta,
    write_diversity_pcoa,
    write_diversity_biplot,
    write_emperor,
    write_emperor_biplot,
    get_subset,
)
from routine_qiime2_analyses._routine_q2_cmds import run_export, run_import


def run_beta(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
             datasets_read: dict, p_beta_subsets: str, trees: dict, force: bool,
             prjct_nm: str, qiime_env: str, chmod: str,
             noloc: bool, Bs: tuple, dropout: bool, filt_raref: str) -> dict:
    """
    Run beta: Beta diversity.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param datasets: list of datasets.
    :param datasets_phylo: phylogenetic decision for each dataset.
    :param trees: phylogenetic trees.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: deta divesity matrices.
    """

    beta_metrics = get_metrics('beta_metrics', Bs)
    beta_subsets = get_subsets(p_beta_subsets)
    job_folder = get_job_folder(i_datasets_folder, 'beta')
    job_folder2 = get_job_folder(i_datasets_folder, 'beta/chunks')

    betas = {}
    main_written = 0
    run_pbs = '%s/2_run_beta%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            written = 0
            tsv, meta = tsv_meta_pds
            if datasets_read[dat] == 'raref':
                if not isfile(tsv):
                    print('Must have run rarefaction to use it further...\nExiting')
                    sys.exit(0)
                tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                datasets_read[dat] = [tsv_pd, meta_pd]
            else:
                tsv_pd, meta_pd = datasets_read[dat]

            if dat not in betas:
                betas[dat] = {}
            out_sh = '%s/run_beta_%s5s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                divs = {}
                # divs = []
                for metric in beta_metrics:
                    if 'unifrac' in metric:
                        if not datasets_phylo[dat][0] or dat not in trees:
                            continue
                    qza = tsv.replace('.tsv', '.qza')
                    if metric not in divs:
                        divs[metric] = {}

                    odir = get_analysis_folder(i_datasets_folder, 'beta/%s' % dat)
                    out_fp = '%s/%s_%s_DM.qza' % (odir, basename(splitext(qza)[0]), metric)
                    if force or not os.path.isfile(out_fp):
                        write_diversity_beta(out_fp, datasets_phylo, trees,
                                             dat, qza, metric, cur_sh, False)
                        # if ret_continue:
                        #     continue
                        written += 1
                        main_written += 1
                    divs[metric][''] = (meta, qza, out_fp)
                    # divs.append(out_fp)

                if beta_subsets and dat in beta_subsets:
                    for subset, subset_regex in beta_subsets[dat].items():
                        subset_done = set()
                        odir = get_analysis_folder(i_datasets_folder, 'beta/%s/%s' % (dat, subset))
                        for metric in beta_metrics:
                            qza_to_subset = tsv.replace('.tsv', '.qza')
                            tsv_to_subset_pd = tsv_pd
                            if 'unifrac' in metric:
                                if not datasets_phylo[dat][0] or dat not in trees:
                                    continue
                                if datasets_phylo[dat][1]:
                                    qza_to_subset = trees[dat][0]
                                    tsv_to_subset = '%s.tsv' % splitext(qza_to_subset)[0]
                                    tsv_to_subset_pd = pd.read_csv(tsv_to_subset, header=0, index_col=0,
                                                             sep='\t', low_memory=False)
                            if dropout:
                                qza_subset = '%s/%s_%s.qza' % (odir, basename(splitext(qza)[0]), subset)
                            else:
                                qza_subset = '%s/%s_%s_noDropout.qza' % (odir, basename(splitext(qza)[0]), subset)
                            tsv_subset = '%s.tsv' % splitext(qza_subset)[0]
                            # meta_subset = '%s.meta' % splitext(qza_subset)[0]
                            subset_feats = get_subset(tsv_to_subset_pd, subset_regex)
                            if not len(subset_feats):
                                continue

                            if tsv_subset not in subset_done:
                                tsv_subset_pd = tsv_to_subset_pd.loc[[x for x in subset_feats if x in tsv_to_subset_pd.index],:].copy()
                                if dropout:
                                    tsv_subset_pd = tsv_subset_pd.loc[:,tsv_subset_pd.sum(0)>0]
                                tsv_subset_pd.to_csv(tsv_subset, index=True, sep='\t')
                                # write_filter_features(qza_to_subset, qza_subset, meta_subset, cur_sh)
                                cmd = run_import(tsv_subset, qza_subset, 'FeatureTable')
                                cur_sh.write('%s\n\n' % cmd)
                                subset_done.add(tsv_subset)
                            out_fp = '%s/%s__%s_DM.qza' % (odir, basename(splitext(qza_subset)[0]), metric)
                            if force or not isfile(out_fp):
                                write_diversity_beta(out_fp, {dat: [1, 0]}, trees,
                                                     dat, qza_subset, metric,
                                                     cur_sh, True)
                                # if ret_continue:
                                #     continue
                                written += 1
                                main_written += 1
                            divs[metric][subset] = (meta, qza_subset, out_fp)
                betas[dat] = divs
            run_xpbs(out_sh, out_pbs, '%s.bt.%s%s' % (prjct_nm, dat, filt_raref),
                     qiime_env, '24', '1', '1', '10', 'gb',
                     chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Calculate beta diversity indices', 'sh', run_pbs)
    return betas


def export_beta(i_datasets_folder: str, betas: dict,
                force: bool, prjct_nm: str, qiime_env: str,
                chmod: str, noloc: bool, filt_raref: str) -> None:
    """
    Export beta diverity matrices.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """

    job_folder = get_job_folder(i_datasets_folder, 'beta')
    out_sh = '%s/2x_run_beta_export%s.sh' % (job_folder, filt_raref)
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    written = 0
    with open(out_sh, 'w') as sh:
        for dat, metric_group_meta_dms in betas.items():
            for metric, group_meta_dms in metric_group_meta_dms.items():
                for group, (meta, qza, dm) in group_meta_dms.items():
                    mat_export = '%s.tsv' % splitext(dm)[0]
                    if force or not isfile(mat_export):
                        cmd = run_export(dm, mat_export, '')
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1
    run_xpbs(out_sh, out_pbs, '%s.xprt.bt%s' % (prjct_nm, filt_raref),
             qiime_env, '2', '1', '1', '1', 'gb',
             chmod, written, '# Export beta diversity matrices', None, noloc)


def run_pcoas(i_datasets_folder: str, betas: dict,
              force: bool, prjct_nm: str, qiime_env: str,
              chmod: str, noloc: bool, filt_raref: str) -> dict:
    """
    Run pcoa: Principal Coordinate Analysis.
    Run pcoa: Principal Coordinate Analysis Biplot¶
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: principal coordinates ordinations.
    """

    job_folder = get_job_folder(i_datasets_folder, 'pcoa')
    job_folder2 = get_job_folder(i_datasets_folder, 'pcoa/chunks')

    pcoas_d = {}
    main_written = 0
    run_pbs = '%s/3_run_pcoa%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, metric_groups_metas_dms in betas.items():
            written = 0
            odir = get_analysis_folder(i_datasets_folder, 'pcoa/%s' % dat)
            if dat not in pcoas_d:
                pcoas_d[dat] = []
            out_sh = '%s/run_PCoA_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for metric, group_meta_dms in metric_groups_metas_dms.items():
                    for group, (meta, qza, dm) in group_meta_dms.items():
                        out = '%s_PCoA.qza' % splitext(dm)[0].replace('/beta/', '/pcoa/')
                        out_tsv = '%s.tsv' % splitext(out)[0]
                        out_dir = os.path.dirname(out)
                        if not os.path.isdir(out_dir):
                            os.makedirs(out_dir)
                        pcoas_d[dat].append((meta, out))
                        if force or not isfile(out) or not isfile(out_tsv):
                            write_diversity_pcoa(dm, out, out_tsv, cur_sh)
                            written += 1
                            main_written += 1
            run_xpbs(out_sh, out_pbs, '%s.pc.%s%s' % (prjct_nm, dat, filt_raref),
                     qiime_env, '10', '1', '2', '2', 'gb',
                     chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Calculate principal coordinates', 'sh', run_pbs)
    return pcoas_d


def run_emperor(i_datasets_folder: str, pcoas_d: dict, prjct_nm: str,
                qiime_env: str, chmod: str, noloc: bool, filt_raref: str) -> None:
    """
    Run emperor.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param pcoas_d: principal coordinates ordinations.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'emperor')
    job_folder2 = get_job_folder(i_datasets_folder, 'emperor/chunks')

    main_written = 0
    first_print = 0
    run_pbs = '%s/4_run_emperor%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, metas_pcoas in pcoas_d.items():
            written = 0
            odir = get_analysis_folder(i_datasets_folder, 'emperor/%s' % dat)
            out_sh = '%s/run_emperor_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for meta_, pcoa in metas_pcoas:
                    meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                    if isfile(meta_alphas):
                        meta = meta_alphas
                    else:
                        meta = meta_
                        if not first_print:
                            print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                                  '\t(if you want alpha diversity as a variable in the PCoA)!')
                            first_print += 1
                    out_plot = '%s_emperor.qzv' % splitext(pcoa)[0].replace('/pcoa/', '/emperor/')
                    out_dir = os.path.dirname(out_plot)
                    if not os.path.isdir(out_dir):
                        os.makedirs(out_dir)
                    write_emperor(meta, pcoa, out_plot, cur_sh)
                    written += 1
                    main_written += 1
            run_xpbs(out_sh, out_pbs, '%s.mprr.%s%s' % (prjct_nm, dat, filt_raref),
                     qiime_env, '10', '1', '1', '1', 'gb',
                     chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Make EMPeror plots', 'sh', run_pbs)


def run_biplots(i_datasets_folder: str, datasets: dict, betas: dict,
                taxonomies: dict, force: bool, prjct_nm: str, qiime_env: str,
                chmod: str, noloc: bool, filt_raref: str) -> dict:
    """
    Run pcoa-biplot: Principal Coordinate Analysis Biplot¶
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa-biplot/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: principal coordinates ordinations.
    """

    job_folder = get_job_folder(i_datasets_folder, 'biplot')
    job_folder2 = get_job_folder(i_datasets_folder, 'biplot/chunks')

    biplots_d = {}
    main_written = 0
    run_pbs = '%s/3_run_biplot%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, metric_groups_metas_dms in betas.items():
            written = 0
            if dat in taxonomies:
                method, tax_qza = taxonomies[dat]
            else:
                tax_qza = 'missing'
            # tsv, meta = datasets[dat]
            odir = get_analysis_folder(i_datasets_folder, 'biplot/%s' % dat)
            if dat not in biplots_d:
                biplots_d[dat] = {}
            out_sh = '%s/run_biplot_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for metric, group_meta_dms in metric_groups_metas_dms.items():
                    for group, (meta, qza, dm) in group_meta_dms.items():
                        tsv = '%s.tsv' % splitext(qza)[0]
                        out_pcoa = '%s_PCoA.qza' % splitext(dm)[0].replace('/beta/', '/pcoa/')
                        out_biplot = '%s_biplot.qza' % splitext(dm)[0].replace('/beta/', '/biplot/')
                        out_dir = os.path.dirname(out_biplot)
                        if not os.path.isdir(out_dir):
                            os.makedirs(out_dir)
                        tsv_tax = '%s_tax.tsv' % splitext(out_biplot)[0]
                        if force or not isfile(out_biplot):
                            write_diversity_biplot(tsv, qza, out_pcoa, out_biplot,
                                                   tax_qza, tsv_tax, cur_sh)
                            written += 1
                            main_written += 1
                        biplots_d[dat].setdefault(meta, []).append((out_biplot, tsv_tax))
            run_xpbs(out_sh, out_pbs, '%s.bplt.%s%s' % (prjct_nm, dat, filt_raref),
                     qiime_env, '10', '1', '2', '2', 'gb',
                     chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Calculate principal coordinates (biplot)', 'sh', run_pbs)
    return biplots_d


def run_emperor_biplot(i_datasets_folder: str, biplots_d: dict, taxonomies: dict,
                prjct_nm: str, qiime_env: str, chmod: str, noloc: bool, filt_raref: str) -> None:
    """
    Run emperor.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param biplot_d: principal coordinates ordinations.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'emperor_biplot')
    job_folder2 = get_job_folder(i_datasets_folder, 'emperor_biplot/chunks')

    main_written = 0
    first_print = 0
    run_pbs = '%s/4_run_emperor_biplot%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, meta_biplots_taxs in biplots_d.items():
            written = 0
            if dat in taxonomies:
                method, tax_qza = taxonomies[dat]
                tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            else:
                tax_tsv = 'missing'
            odir = get_analysis_folder(i_datasets_folder, 'emperor_biplot/%s' % dat)
            out_sh = '%s/run_emperor_biplot_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for meta_, biplots_taxs in meta_biplots_taxs.items():
                    meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                    if isfile(meta_alphas):
                        meta = meta_alphas
                    else:
                        meta = meta_
                        if not first_print:
                            print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                                  '\t(if you want alpha diversity as a variable in the PCoA biplot)!')
                            first_print += 1
                    for biplot, tsv_tax in biplots_taxs:
                        out_plot = '%s_emperor_biplot.qzv' % splitext(biplot)[0].replace('/biplot/', '/emperor_biplot/')
                        out_dir = os.path.dirname(out_plot)
                        if not os.path.isdir(out_dir):
                            os.makedirs(out_dir)
                        if isfile(tsv_tax):
                            write_emperor_biplot(meta, biplot, out_plot, cur_sh, tsv_tax)
                        else:
                            write_emperor_biplot(meta, biplot, out_plot, cur_sh, tax_tsv)
                        written += 1
                        main_written += 1
            run_xpbs(out_sh, out_pbs, '%s.mprr.bplt.%s%s' % (prjct_nm, dat, filt_raref),
                     qiime_env, '10', '1', '1', '1', 'gb',
                     chmod, written, 'single', o, noloc)
    if main_written:
        print_message('# Make EMPeror biplots', 'sh', run_pbs)
