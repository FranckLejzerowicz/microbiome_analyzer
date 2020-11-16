# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, sys
import pandas as pd
from os.path import basename, dirname, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics,
    get_job_folder,
    get_analysis_folder,
    read_yaml_file,
    get_raref_tab_meta_pds
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_diversity_beta,
    write_diversity_pcoa,
    write_diversity_biplot,
    write_emperor,
    write_empress,
    write_emperor_biplot,
    write_empress_biplot,
    get_subset,
    run_export,
    run_import
)


def run_beta(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
             datasets_read: dict, datasets_rarefs: dict, p_beta_subsets: str,
             trees: dict, force: bool, prjct_nm: str, qiime_env: str, chmod: str,
             noloc: bool, Bs: tuple, dropout: bool, run_params: dict,
             filt_raref: str, eval_depths: dict, jobs: bool, chunkit: int) -> dict:
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
    evaluation = ''
    if len(eval_depths):
        evaluation = '_eval'
    beta_metrics = get_metrics('beta_metrics', Bs)
    beta_subsets = read_yaml_file(p_beta_subsets)
    job_folder = get_job_folder(i_datasets_folder, 'beta%s' % evaluation)
    job_folder2 = get_job_folder(i_datasets_folder, 'beta%s/chunks' % evaluation)

    betas = {}
    main_written = 0
    run_pbs = '%s/2_run_beta%s%s.sh' % (job_folder, evaluation, filt_raref)
    with open(run_pbs, 'w') as o:

        for dat, tsv_meta_pds_ in datasets.items():
            written = 0
            betas[dat] = []
            out_sh = '%s/run_beta%s_%s%s.sh' % (job_folder2, evaluation, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
                    tsv, meta = tsv_meta_pds
                    if not isinstance(datasets_read[dat][idx][0], pd.DataFrame) and datasets_read[dat][idx][0] == 'raref':
                        if not isfile(tsv):
                            print('Must have run rarefaction to use it further...\nExiting')
                            sys.exit(0)
                        tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                        datasets_read[dat][idx] = [tsv_pd, meta_pd]
                    else:
                        tsv_pd, meta_pd = datasets_read[dat][idx]
                    cur_raref = datasets_rarefs[dat][idx]
                    divs = {}
                    for metric in beta_metrics:
                        if 'unifrac' in metric:
                            if not datasets_phylo[dat][0] or dat not in trees:
                                continue
                        qza = tsv.replace('.tsv', '.qza')
                        if metric not in divs:
                            divs[metric] = {}

                        odir = get_analysis_folder(i_datasets_folder, 'beta%s/%s%s' % (evaluation, dat, cur_raref))
                        out_fp = '%s/%s_%s_DM.qza' % (odir, basename(splitext(qza)[0]), metric)
                        if force or not os.path.isfile(out_fp):
                            tree = write_diversity_beta(out_fp, datasets_phylo, trees,
                                                        dat, qza, metric, cur_sh, False)
                            written += 1
                            main_written += 1
                        else:
                            tree = ''
                            if 'unifrac' in metric:
                                tree = trees[dat][1]
                        divs[metric][''] = (meta, qza, out_fp, tree)

                    if beta_subsets and dat in beta_subsets:
                        for subset, subset_regex in beta_subsets[dat].items():
                            subset_done = set()
                            odir = get_analysis_folder(i_datasets_folder, 'beta%s/%s%s/%s' % (
                                evaluation, dat, cur_raref, subset))
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

                                subset_feats = get_subset(tsv_to_subset_pd, subset_regex)
                                if not len(subset_feats):
                                    continue

                                if tsv_subset not in subset_done:
                                    tsv_subset_pd = tsv_to_subset_pd.loc[[x for x in subset_feats if x in tsv_to_subset_pd.index],:].copy()
                                    if dropout:
                                        tsv_subset_pd = tsv_subset_pd.loc[:,tsv_subset_pd.sum(0)>0]
                                    tsv_subset_pd.to_csv(tsv_subset, index=True, sep='\t')

                                    cmd = run_import(tsv_subset, qza_subset, 'FeatureTable')
                                    cur_sh.write('%s\n\n' % cmd)
                                    subset_done.add(tsv_subset)
                                out_fp = '%s/%s__%s_DM.qza' % (odir, basename(splitext(qza_subset)[0]), metric)
                                if force or not isfile(out_fp):
                                    tree = write_diversity_beta(out_fp, {dat: [1, 0]}, trees,
                                                                dat, qza_subset, metric,
                                                                cur_sh, True)
                                    written += 1
                                    main_written += 1
                                else:
                                    tree = ''
                                    if 'unifrac' in metric:
                                        tree = trees[dat][1]

                                divs[metric][subset] = (meta, qza_subset, out_fp, tree)
                    betas[dat].append(divs)
            run_xpbs(out_sh, out_pbs, '%s.bt%s.%s%s' % (prjct_nm, evaluation, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Calculate beta diversity indices', 'sh', run_pbs, jobs)
    return betas


def export_beta(i_datasets_folder: str, betas: dict, datasets_rarefs: dict,
                force: bool, prjct_nm: str, qiime_env: str, chmod: str,
                noloc: bool, run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
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
        for dat, metric_group_meta_dms_ in betas.items():
            for idx, metric_group_meta_dms in enumerate(metric_group_meta_dms_):
                for metric, group_meta_dms in metric_group_meta_dms.items():
                    for group, (meta, qza, dm, tree) in group_meta_dms.items():
                        mat_export = '%s.tsv' % splitext(dm)[0]
                        if force or not isfile(mat_export):
                            cmd = run_export(dm, mat_export, '')
                            sh.write('echo "%s"\n' % cmd)
                            sh.write('%s\n\n' % cmd)
                            written += 1
    run_xpbs(out_sh, out_pbs, '%s.xprt.bt%s' % (prjct_nm, filt_raref), qiime_env,
             run_params["time"], run_params["n_nodes"], run_params["n_procs"],
             run_params["mem_num"], run_params["mem_dim"],
             chmod, written, '# Export beta diversity matrices', None, noloc, jobs)


def run_pcoas(i_datasets_folder: str, betas: dict, datasets_rarefs: dict,
              force: bool, prjct_nm: str, qiime_env: str, chmod: str,
              noloc: bool, run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> dict:
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
        for dat, metric_groups_metas_dms_ in betas.items():
            written = 0
            pcoas_d[dat] = []
            out_sh = '%s/run_PCoA_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, metric_groups_metas_dms in enumerate(metric_groups_metas_dms_):
                    dat_pcoas = []
                    cur_depth = datasets_rarefs[dat][idx]
                    get_analysis_folder(i_datasets_folder, 'pcoa/%s%s' % (dat, cur_depth))
                    for metric, group_meta_dms in metric_groups_metas_dms.items():
                        for group, (meta, qza, dm, tree) in group_meta_dms.items():
                            out = '%s_PCoA.qza' % splitext(dm)[0].replace('/beta/', '/pcoa/')
                            out_tsv = '%s.tsv' % splitext(out)[0]
                            out_dir = os.path.dirname(out)
                            if not os.path.isdir(out_dir):
                                os.makedirs(out_dir)
                            dat_pcoas.append((meta, out, qza, tree))
                            if force or not isfile(out) or not isfile(out_tsv):
                                write_diversity_pcoa(dm, out, out_tsv, cur_sh)
                                written += 1
                                main_written += 1
                    pcoas_d[dat].append(dat_pcoas)
            run_xpbs(out_sh, out_pbs, '%s.pc.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Calculate principal coordinates', 'sh', run_pbs, jobs)
    return pcoas_d


def run_emperor(i_datasets_folder: str, pcoas_d: dict, datasets_rarefs: dict,
                prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
                run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
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
        for dat, metas_pcoas_ in pcoas_d.items():
            written = 0
            out_sh = '%s/run_emperor_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, metas_pcoas in enumerate(metas_pcoas_):
                    cur_depth = datasets_rarefs[dat][idx]
                    get_analysis_folder(i_datasets_folder, 'emperor/%s%s' % (dat, cur_depth))
                    for meta_, pcoa, _, __ in metas_pcoas:
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
            run_xpbs(out_sh, out_pbs, '%s.mprr.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Make EMPeror plots', 'sh', run_pbs, jobs)


def run_biplots(i_datasets_folder: str, betas: dict, datasets_rarefs: dict,
                taxonomies: dict, force: bool, prjct_nm: str, qiime_env: str,
                chmod: str, noloc: bool, run_params: dict,
                filt_raref: str, jobs: bool, chunkit: int) -> (dict, dict):
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
    biplots_d2 = {}
    main_written = 0
    run_pbs = '%s/3_run_biplot%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, metric_groups_metas_dms_ in betas.items():
            written = 0
            if dat in taxonomies:
                method, tax_qza, tax_tsv = taxonomies[dat]
            else:
                tax_qza = 'missing'
            biplots_d[dat] = []
            biplots_d2[dat] = []
            out_sh = '%s/run_biplot_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, metric_groups_metas_dms in enumerate(metric_groups_metas_dms_):
                    dat_biplots = {}
                    dat_biplots2 = {}
                    cur_depth = datasets_rarefs[dat][idx]
                    get_analysis_folder(i_datasets_folder, 'biplot/%s%s' % (dat, cur_depth))
                    for metric, group_meta_dms in metric_groups_metas_dms.items():
                        for group, (meta, qza, dm, tree) in group_meta_dms.items():
                            tsv = '%s.tsv' % splitext(qza)[0]
                            out_pcoa = '%s_PCoA.qza' % splitext(dm)[0].replace('/beta/', '/pcoa/')
                            out_biplot = '%s_biplot.qza' % splitext(dm)[0].replace('/beta/', '/biplot/')
                            out_biplot2 = '%s_biplot_raw.qza' % splitext(dm)[0].replace('/beta/', '/biplot/')
                            out_dir = os.path.dirname(out_biplot)
                            if not os.path.isdir(out_dir):
                                os.makedirs(out_dir)
                            tsv_tax = '%s_tax.tsv' % splitext(out_biplot)[0]
                            if force or not isfile(out_biplot) or not isfile(out_biplot2):
                                write_diversity_biplot(tsv, qza, out_pcoa, out_biplot,
                                                       out_biplot2, tax_qza, tsv_tax, cur_sh)
                                written += 1
                                main_written += 1
                            dat_biplots.setdefault(meta, []).append((out_biplot, tsv_tax, qza, tree))
                            dat_biplots2.setdefault(meta, []).append((out_biplot2, tsv_tax, qza, tree))
                    biplots_d[dat].append(dat_biplots)
                    biplots_d2[dat].append(dat_biplots2)
            run_xpbs(out_sh, out_pbs, '%s.bplt.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Calculate principal coordinates (biplot)', 'sh', run_pbs, jobs)
    return biplots_d, biplots_d2


def run_emperor_biplot(i_datasets_folder: str, biplots_d: dict, biplots_d2: dict,
                       taxonomies: dict, split_taxa_pds: dict,  datasets_rarefs: dict,
                       prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
                       run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
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
        for dat, raref_meta_biplots_taxs_qzas_trees in biplots_d.items():
            raref_meta_biplots_taxs_qzas_trees2 = biplots_d2[dat]
            written = 0
            if dat in taxonomies:
                method, tax_qza, tax_tsv = taxonomies[dat]
                split_taxa_pd = split_taxa_pds[dat]
            else:
                tax_tsv = 'missing'
            out_sh = '%s/run_emperor_biplot_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, meta_biplots_taxs_qzas_trees in enumerate(raref_meta_biplots_taxs_qzas_trees):
                    meta_biplots_taxs_qzas_trees2 = raref_meta_biplots_taxs_qzas_trees2[idx]
                    cur_raref = datasets_rarefs[dat][idx]
                    get_analysis_folder(i_datasets_folder, 'emperor_biplot/%s%s' % (dat, cur_raref))
                    for meta_, biplots_taxs_qzas_trees in meta_biplots_taxs_qzas_trees.items():
                        biplots_taxs_qzas_trees2 = meta_biplots_taxs_qzas_trees2[meta_]
                        meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                        if isfile(meta_alphas):
                            meta = meta_alphas
                        else:
                            meta = meta_
                            if not first_print:
                                print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                                      '\t(if you want alpha diversity as a variable in the PCoA biplot)!')
                                first_print += 1
                        for bdx, (biplot, tsv_tax, _, __) in enumerate(biplots_taxs_qzas_trees):
                            biplot2, tsv_tax2 = biplots_taxs_qzas_trees2[bdx][:2]
                            out_plot = '%s_emperor_biplot.qzv' % splitext(biplot)[0].replace('/biplot/', '/emperor_biplot/')
                            out_dir = dirname(out_plot)
                            if not os.path.isdir(out_dir):
                                os.makedirs(out_dir)
                            if isfile(tsv_tax):
                                write_emperor_biplot(meta, biplot, out_plot, cur_sh, tsv_tax, split_taxa_pd)
                            else:
                                write_emperor_biplot(meta, biplot, out_plot, cur_sh, tax_tsv, {})
                            if isfile(tsv_tax2):
                                write_emperor_biplot(meta, biplot2, out_plot, cur_sh, tsv_tax2, split_taxa_pd)
                            else:
                                write_emperor_biplot(meta, biplot2, out_plot, cur_sh, tax_tsv, {})
                            written += 1
                            main_written += 1
            run_xpbs(out_sh, out_pbs, '%s.mprr.bplt.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Make EMPeror biplots', 'sh', run_pbs, jobs)


def run_empress(i_datasets_folder: str, pcoas_d: dict,
                trees: dict, datasets_phylo: dict, datasets_rarefs: dict,
                taxonomies: dict, prjct_nm: str, qiime_env: str, chmod: str,
                noloc: bool, run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
    """
    Run empress.
    https://docs.qiime2.org/2019.10/plugins/available/empress/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param pcoas_d: principal coordinates ordinations.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (default: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'empress')
    job_folder2 = get_job_folder(i_datasets_folder, 'empress/chunks')

    main_written = 0
    first_print = 0
    run_pbs = '%s/4_run_empress%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, metas_pcoas_qzas_trees_ in pcoas_d.items():
            written = 0
            out_sh = '%s/run_empress_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            if not datasets_phylo[dat][0] or dat not in trees:
                continue
            tax_qza = ''
            if dat in taxonomies:
                method, tax_qza, tax_tsv = taxonomies[dat]

            with open(out_sh, 'w') as cur_sh:
                for idx, metas_pcoas_qzas_trees in enumerate(metas_pcoas_qzas_trees_):
                    cur_depth = datasets_rarefs[dat][idx]

                    odir = get_analysis_folder(i_datasets_folder, 'songbird/%s%s' % (dat, cur_depth))
                    sb_qza = '%s/sb_%s.qza' % (odir, dat)
                    if not isfile(sb_qza):
                        sb_qza = ''

                    get_analysis_folder(i_datasets_folder, 'empress/%s%s' % (dat, cur_depth))
                    for meta_, pcoa, qza, tree in metas_pcoas_qzas_trees:
                        if tree:
                            meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                            if isfile(meta_alphas):
                                sam_meta = meta_alphas
                            else:
                                sam_meta = meta_
                                if not first_print:
                                    print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                                          '\t(if you want alpha diversity as a variable in the PCoA)!')
                                    first_print += 1
                            out_plot = '%s_empress.qzv' % splitext(pcoa)[0].replace('/pcoa/', '/empress/')
                            out_dir = os.path.dirname(out_plot)
                            if not os.path.isdir(out_dir):
                                os.makedirs(out_dir)
                            write_empress(sam_meta, qza, tax_qza, sb_qza, pcoa, tree, out_plot, cur_sh)
                            written += 1
                            main_written += 1
            run_xpbs(out_sh, out_pbs, '%s.mprss.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Make empress plots', 'sh', run_pbs, jobs)


def run_empress_biplot(i_datasets_folder: str, biplots_d: dict, biplots_d2: dict,
                       trees: dict, datasets_phylo: dict, taxonomies: dict,
                       datasets_rarefs: dict, prjct_nm: str, qiime_env: str,
                       chmod: str, noloc: bool, run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
    """
    Run empress.
    https://docs.qiime2.org/2019.10/plugins/available/empress/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param biplot_d: principal coordinates ordinations.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'empress_biplot')
    job_folder2 = get_job_folder(i_datasets_folder, 'empress_biplot/chunks')

    main_written = 0
    first_print = 0
    run_pbs = '%s/4_run_empress_biplot%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, raref_meta_biplots_taxs_qzas_trees in biplots_d.items():
            raref_meta_biplots_taxs_qzas_trees2 = biplots_d2[dat]
            written = 0
            if not datasets_phylo[dat][0] or dat not in trees:
                continue
            tax_qza = ''
            if dat in taxonomies:
                method, tax_qza, tax_tsv = taxonomies[dat]

            out_sh = '%s/run_empress_biplot_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for idx, meta_biplots_taxs_qzas_trees in enumerate(raref_meta_biplots_taxs_qzas_trees):
                    meta_biplots_taxs_qzas_trees2 = raref_meta_biplots_taxs_qzas_trees2[idx]
                    cur_raref = datasets_rarefs[dat][idx]
                    get_analysis_folder(i_datasets_folder, 'empress_biplot/%s%s' % (dat, cur_raref))

                    odir = get_analysis_folder(i_datasets_folder, 'songbird/%s%s' % (dat, cur_raref))
                    sb_qza = '%s/sb_%s.qza' % (odir, dat)
                    if not isfile(sb_qza):
                        sb_qza = ''

                    for meta_, biplots_taxs_qzas_trees in meta_biplots_taxs_qzas_trees.items():
                        biplots_taxs_qzas_trees2 = meta_biplots_taxs_qzas_trees2[meta_]
                        meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                        if isfile(meta_alphas):
                            meta = meta_alphas
                        else:
                            meta = meta_
                            if not first_print:
                                print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                                      '\t(if you want alpha diversity as a variable in the PCoA biplot)!')
                                first_print += 1
                        for bdx, (biplot, tsv_tax, qza, tree) in enumerate(biplots_taxs_qzas_trees):
                            biplot2, tsv_tax2, qza2, tree2 = biplots_taxs_qzas_trees2[idx]

                            if tree:
                                out_plot = '%s_empress_biplot.qzv' % splitext(biplot)[0].replace(
                                    '/biplot/', '/empress_biplot/')
                                out_dir = os.path.dirname(out_plot)
                                if not os.path.isdir(out_dir):
                                    os.makedirs(out_dir)
                                if isfile(tsv_tax):
                                    write_empress_biplot(meta, qza, tax_qza, sb_qza, biplot,
                                                         tree, out_plot, cur_sh)
                                else:
                                    write_empress_biplot(meta, qza, tax_qza, sb_qza, biplot,
                                                         tree, out_plot, cur_sh)
                                written += 1
                                main_written += 1
                            if tree2:
                                out_plot2 = '%s_empress_biplot_raw.qzv' % splitext(biplot2)[0].replace(
                                    '/biplot/', '/empress_biplot/')
                                out_dir2 = os.path.dirname(out_plot2)
                                if not os.path.isdir(out_dir2):
                                    os.makedirs(out_dir2)
                                if isfile(tsv_tax2):
                                    write_empress_biplot(meta, qza2, tax_qza, sb_qza, biplot2,
                                                         tree2, out_plot2, cur_sh)
                                else:
                                    write_empress_biplot(meta, qza2, tax_qza, sb_qza, biplot2,
                                                         tree2, out_plot2, cur_sh)
                                written += 1
                                main_written += 1

            run_xpbs(out_sh, out_pbs, '%s.mprss.bplt.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if main_written:
        print_message('# Make empress biplots', 'sh', run_pbs, jobs)
