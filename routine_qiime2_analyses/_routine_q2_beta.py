# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import basename, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics,
    get_job_folder,
    get_analysis_folder
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_diversity_beta,
    write_diversity_pcoa,
    write_diversity_biplot,
    write_emperor
)
from routine_qiime2_analyses._routine_q2_cmds import run_export


def run_beta(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
             trees: dict, force: bool, prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> dict:
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

    beta_metrics = get_metrics('beta_metrics')
    job_folder = get_job_folder(i_datasets_folder, 'beta')
    job_folder2 = get_job_folder(i_datasets_folder, 'beta/chunks')

    betas = {}
    written = 0
    run_pbs = '%s/2_run_beta.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            if dat not in betas:
                betas[dat] = {}
            out_sh = '%s/run_beta_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                qza = tsv.replace('.tsv', '.qza')
                divs = []
                for metric in beta_metrics:
                    odir = get_analysis_folder(i_datasets_folder, 'beta/%s' % dat)
                    out_fp = '%s/%s_%s_DM.qza' % (odir, basename(splitext(qza)[0]), metric)
                    if force or not os.path.isfile(out_fp):
                        ret_continue = write_diversity_beta(out_fp, datasets_phylo, trees,
                                                            dat, qza, metric, cur_sh)
                        if ret_continue:
                            continue
                        written += 1
                    divs.append(out_fp)
                betas[dat][meta] = divs
            run_xpbs(out_sh, out_pbs, '%s.bt.%s' % (prjct_nm, dat),
                     qiime_env, '24', '1', '1', '10', 'gb',
                     chmod, written, 'single', o, noloc)
    if written:
        print_message('# Calculate beta diversity indices', 'sh', run_pbs)
    return betas


def export_beta(i_datasets_folder: str, betas: dict,
                force: bool, prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> None:
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
    out_sh = '%s/2x_run_beta_export.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    written = 0
    with open(out_sh, 'w') as sh:
        for dat, meta_mats in betas.items():
            for meta, mats in meta_mats.items():
                for mat in mats:
                    mat_export = '%s.tsv' % splitext(mat)[0]
                    if force or not isfile(mat_export):
                        cmd = run_export(mat, mat_export, '')
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1
    run_xpbs(out_sh, out_pbs, '%s.xprt.bt' % prjct_nm,
             qiime_env, '2', '1', '1', '1', 'gb',
             chmod, written, '# Export beta diversity matrices', noloc)


def run_pcoas_biplots(i_datasets_folder: str, datasets: dict, betas: dict,
              force: bool, prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> tuple:
    """
    Run pcoa: Principal Coordinate Analysis.
    Run pcoa-biplot: Principal Coordinate Analysis Biplot¶
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa-biplot/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: principal coordinates ordinations.
    """

    job_folder_pcoa = get_job_folder(i_datasets_folder, 'pcoa')
    job_folder2_pcoa = get_job_folder(i_datasets_folder, 'pcoa/chunks')
    job_folder_biplot = get_job_folder(i_datasets_folder, 'biplot')
    job_folder2_biplot = get_job_folder(i_datasets_folder, 'biplot/chunks')

    pcoas_d = {}
    biplots_d = {}
    written_pcoa = 0
    written_biplot = 0
    run_pcoa_pbs = '%s/3_run_pcoa.sh' % job_folder_pcoa
    run_biplot_pbs = '%s/3_run_biplot.sh' % job_folder_biplot
    with open(run_pcoa_pbs, 'w') as o_pcoa, open(run_biplot_pbs, 'w') as o_biplot:
        for dat, meta_DMs in betas.items():
            tsv, meta = datasets[dat]
            qza = '%s.qza' % splitext(tsv)[0]
            odir_pcoa = get_analysis_folder(i_datasets_folder, 'pcoa/%s' % dat)
            odir_biplot = get_analysis_folder(i_datasets_folder, 'biplot/%s' % dat)
            if dat not in pcoas_d:
                pcoas_d[dat] = {}
                biplots_d[dat] = {}
            out_pcoa_sh = '%s/run_PCoA_%s.sh' % (job_folder2_pcoa, dat)
            out_pcoa_pbs = '%s.pbs' % splitext(out_pcoa_sh)[0]
            out_biplot_sh = '%s/run_biplot_%s.sh' % (job_folder2_biplot, dat)
            out_biplot_pbs = '%s.pbs' % splitext(out_biplot_sh)[0]
            with open(out_pcoa_sh, 'w') as cur_pcoa_sh, open(out_biplot_sh, 'w') as cur_biplot_sh:
                for meta, DMs in meta_DMs.items():
                    for DM in DMs:
                        out_pcoa = '%s_PCoA.qza' % splitext(DM)[0].replace('/beta/', '/pcoa/')
                        out_biplot = '%s_biplot.qza' % splitext(DM)[0].replace('/beta/', '/biplot/')
                        pcoas_d[dat].setdefault(meta, []).append(out_pcoa)
                        biplots_d[dat].setdefault(meta, []).append(out_biplot)
                        if force or not isfile(out_pcoa):
                            write_diversity_pcoa(DM, out_pcoa, cur_pcoa_sh)
                            written_pcoa += 1
                        if force or not isfile(out_biplot):
                            write_diversity_biplot(qza, out_pcoa, out_biplot, cur_biplot_sh)
                            written_biplot += 1
            run_xpbs(out_pcoa_sh, out_pcoa_pbs, '%s.pc.%s' % (prjct_nm, dat),
                     qiime_env, '10', '1', '2', '2', 'gb',
                     chmod, written_pcoa, 'single', o_pcoa, noloc)
            run_xpbs(out_biplot_sh, out_biplot_pbs, '%s.bplt.%s' % (prjct_nm, dat),
                     qiime_env, '10', '1', '2', '2', 'gb',
                     chmod, written_biplot, 'single', o_biplot, noloc)
    if written_pcoa:
        print_message('# Calculate principal coordinates', 'sh', run_pcoa_pbs)
    if written_biplot:
        print_message('# Calculate principal coordinates (biplot)', 'sh', run_biplot_pbs)
    return pcoas_d, biplots_d


def run_emperor(i_datasets_folder: str, pcoas_biplots_d: dict, taxonomies: dict,
                prjct_nm: str, qiime_env: str, chmod: str, biplot: bool, noloc: bool) -> None:
    """
    Run emperor.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param pcoas_biplot_d: principal coordinates ordinations.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    suffix = ''
    if biplot:
        suffix = '_biplot'
    job_folder = get_job_folder(i_datasets_folder, 'emperor%s' % suffix)
    job_folder2 = get_job_folder(i_datasets_folder, 'emperor%s/chunks' % suffix)

    written = 0
    first_print = 0
    run_pbs = '%s/4_run_emperor.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, meta_pcoas_biplots in pcoas_biplots_d.items():
            if dat in taxonomies:
                method, tax_qza = taxonomies[dat]
            odir = get_analysis_folder(i_datasets_folder, 'emperor%s/%s' % (suffix, dat))
            out_sh = '%s/run_emperor_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for meta_, pcoas_biplots in meta_pcoas_biplots.items():
                    meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                    if isfile(meta_alphas):
                        meta = meta_alphas
                    else:
                        meta = meta_
                        if not first_print:
                            print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                                  '\t(if you want alpha diversity as a variable in the PCoA)!')
                            first_print += 1
                    for pcoa_biplot in pcoas_biplots:
                        if biplot:
                            out_plot = '%s_emperor.qzv' % splitext(pcoa_biplot)[0].replace('/biplot/', '/emperor%s/' % suffix)
                        else:
                            out_plot = '%s_emperor.qzv' % splitext(pcoa_biplot)[0].replace('/pcoa/', '/emperor%s/' % suffix)
                        write_emperor(meta, pcoa_biplot, out_plot, cur_sh, tax_qza)
                        written += 1
            run_xpbs(out_sh, out_pbs, '%s.mprr.%s.%s' % (prjct_nm, suffix, dat),
                     qiime_env, '10', '1', '1', '1', 'gb',
                     chmod, written, 'single', o, noloc)
    if written:
        print_message('# Make EMPeror plots%s' % suffix.replace('_', ' '), 'sh', run_pbs)
