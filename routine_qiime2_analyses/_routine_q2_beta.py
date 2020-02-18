# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import basename, dirname, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics,
    get_job_folder,
    get_analysis_folder
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_diversity_beta,
    write_diversity_pcoa,
    write_emperor
)
from routine_qiime2_analyses._routine_q2_cmds import run_export


def run_beta(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
             trees: dict, force: bool, prjct_nm: str, qiime_env: str) -> dict:
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
                     qiime_env, '24', '1', '1', '10', 'gb', written, 'single', o)
    if written:
        print('# Calculate beta diversity indices')
        print('[TO RUN] sh', run_pbs)
    return betas


def export_beta(i_datasets_folder: str, betas: dict,
                force: bool, prjct_nm: str, qiime_env: str) -> None:
    """
    Export beta diverity matrices.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
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
             qiime_env, '2', '1', '1', '1', 'gb', written,
             '# Export beta diversity matrices')


def run_pcoas(i_datasets_folder: str, betas: dict,
              force: bool, prjct_nm: str, qiime_env: str) -> dict:
    """
    Run pcoa: Principal Coordinate Analysis.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param betas: beta diversity matrices.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :return: principal coordinates ordinations.
    """

    job_folder = get_job_folder(i_datasets_folder, 'pcoa')
    job_folder2 = get_job_folder(i_datasets_folder, 'pcoa/chunks')
    odir = get_analysis_folder(i_datasets_folder, 'pcoa')

    pcoas_d = {}
    written = 0
    run_pbs = '%s/3_run_pcoa.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, meta_DMs in betas.items():
            if dat not in pcoas_d:
                pcoas_d[dat] = {}
            out_sh = '%s/run_PCoA_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for meta, DMs in meta_DMs.items():
                    for DM in DMs:
                        out_pcoa = '%s_PCoA.qza' % splitext(DM)[0].replace('/beta/', '/pcoa/')
                        if force or not isfile(out_pcoa):
                            pcoas_d[dat].setdefault(meta, []).append(out_pcoa)
                            if not isdir(dirname(out_pcoa)):
                                os.makedirs(dirname(out_pcoa))
                            write_diversity_pcoa(DM, out_pcoa, cur_sh)
                            written += 1
            run_xpbs(out_sh, out_pbs, '%s.pc.%s' % (prjct_nm, dat),
                     qiime_env, '10', '1', '2', '2', 'gb', written, 'single', o)
    if written:
        print('# Calculate principal coordinates')
        print('[TO RUN] sh', run_pbs)
    return pcoas_d


def run_emperor(i_datasets_folder: str, pcoas_d: dict,
                prjct_nm: str, qiime_env: str) -> None:
    """
    Run emperor.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param pcoas_d: principal coordinates ordinations.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    """

    job_folder = get_job_folder(i_datasets_folder, 'emperor')
    job_folder2 = get_job_folder(i_datasets_folder, 'emperor/chunks')
    odir = get_analysis_folder(i_datasets_folder, 'emperor')

    written = 0
    first_print = 0
    run_pbs = '%s/4_run_emperor.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, meta_pcoas in pcoas_d.items():
            for meta_, pcoas in meta_pcoas.items():
                meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                if not isfile(meta_alphas):
                    meta = meta_
                    if not first_print:
                        print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                              '\t(if you want alpha diversity as a variable in the PCoA)!')
                        first_print += 1
                else:
                    meta = meta_alphas
                for pcoa in pcoas:
                    out_sh = '%s/run_emperor_%s_%s.sh' % (job_folder2, dat, basename(pcoa.split('_DM_')[0]))
                    out_pbs = '%s.pbs' % splitext(out_sh)[0]
                    with open(out_sh, 'w') as cur_sh:
                        out_plot = '%s_emperor.qzv' % splitext(pcoa)[0].replace('/pcoa/', '/emperor/')
                        write_emperor(meta, pcoa, out_plot, cur_sh)
                        written += 1
                    run_xpbs(out_sh, out_pbs, '%s.mprr.%s.%s' % (prjct_nm, dat, basename(pcoa.split('_DM_')[0])),
                             qiime_env, '10', '1', '1', '2', 'gb', written, 'single', o)
    if written:
        print('# Make EMPeror plots')
        print('[TO RUN] sh', run_pbs)
