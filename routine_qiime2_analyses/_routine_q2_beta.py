# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import basename, dirname, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder


def run_beta(i_datasets_folder: str, datasets: dict, beta_metrics: list,
             wol_trees: dict, force: bool, prjct_nm: str, qiime_env: str) -> dict:

    print('# Calculate beta diversity indices')
    job_folder = get_job_folder(i_datasets_folder, 'beta')
    job_folder2 = get_job_folder(i_datasets_folder, 'beta/chunks')

    betas = {}
    written = 0
    run_pbs = '%s/2_run_beta.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dataset, files_lists in datasets.items():
            if dataset not in betas:
                betas[dataset] = {}
            out_sh = '%s/run_beta_%s.sh' % (job_folder2, dataset)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                for files_list in files_lists:
                    qza = files_list[0].replace('.tsv', '.qza')
                    meta = files_list[1]
                    divs = []
                    for metric in beta_metrics:
                        odir = get_analysis_folder(i_datasets_folder, 'beta/%s' % dataset)
                        out_fp = '%s/%s_%s_DM.qza' % (
                            odir, basename(splitext(qza)[0]), metric)
                        divs.append(out_fp)
                        if force or not os.path.isfile(out_fp):
                            if 'unifrac' in metric:
                                cmd = [
                                    'qiime', 'diversity', 'beta-phylogenetic',
                                    '--i-table', qza,
                                    '--i-phylogeny', wol_trees[dataset],
                                    '--p-metric', metric,
                                    '--p-n-jobs', '1',
                                    '--o-distance-matrix', out_fp]
                            else:
                                cmd = [
                                    'qiime', 'diversity', 'beta',
                                    '--i-table', qza,
                                    '--p-metric', metric,
                                    '--p-n-jobs', '1',
                                    '--o-distance-matrix', out_fp]
                            sh.write('echo "%s"\n' % ' '.join(cmd))
                            sh.write('%s\n' % ' '.join(cmd))
                            written += 1
                    betas[dataset][meta] = divs
            if written:
                xpbs_call(out_sh, out_pbs, '%s.mg.bt.%s' % (prjct_nm, dataset), qiime_env,
                          '24', '1', '1', '10', 'gb')
                o.write('qsub %s\n' % out_pbs)
            else:
                print('\nNothing written in', out_sh, '--> removed')
                os.remove(out_sh)
    print('[TO RUN] sh', run_pbs)
    return betas


def export_beta(i_datasets_folder: str, betas: dict,
                force: bool, prjct_nm: str, qiime_env: str) -> None:

    print('# Export beta diversity matrices')
    job_folder = get_job_folder(i_datasets_folder, 'beta')
    out_sh = '%s/2x_run_beta_export.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    written = 0
    with open(out_sh, 'w') as sh:
        for dat, meta_mats in betas.items():
            for meta, mats in meta_mats.items():
                for mat in mats:
                    out_dir = splitext(mat)[0]
                    mat_export = '%s.tsv' % out_dir
                    if force or not isfile(mat_export):
                        sh.write('\nqiime tools export --input-path %s --output-path %s\n' % (mat, out_dir))
                        sh.write('mv %s/distance-matrix.tsv %s\n' % (out_dir, mat_export))
                        sh.write('rm -rf %s\n' % out_dir)
                        written += 1
    if written:
        xpbs_call(out_sh, out_pbs, '%s.xprt.bt' % prjct_nm, qiime_env,
                  '2', '1', '1', '1', 'gb')
        print('[TO RUN] qsub', out_pbs)
    else:
        print('\nNothing written in', out_sh, '--> removed')
        os.remove(out_sh)


def run_pcoas(i_datasets_folder: str, betas: dict,
              force: bool, prjct_nm: str, qiime_env: str) -> dict:

    print('# Calculate principal coordinates')
    job_folder = get_job_folder(i_datasets_folder, 'pcoa')
    job_folder2 = get_job_folder(i_datasets_folder, 'pcoa/chunks')
    odir = get_analysis_folder(i_datasets_folder, 'pcoa')

    pcoas_d = {}
    written = 0
    run_pbs = '%s/3_run_pcoa.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dataset, meta_DMs in betas.items():
            if dataset not in pcoas_d:
                pcoas_d[dataset] = {}
            out_sh = '%s/run_PCoA_%s.sh' % (job_folder2, dataset)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                for meta, DMs in meta_DMs.items():
                    for DM in DMs:
                        out_pcoa = '%s_PCoA.qza' % splitext(DM)[0].replace('/beta/', '/pcoa/')
                        if force or not isfile(out_pcoa):
                            pcoas_d[dataset].setdefault(meta, []).append(out_pcoa)
                            if not isdir(dirname(out_pcoa)):
                                os.makedirs(dirname(out_pcoa))
                            cmd = [
                                'qiime', 'diversity', 'pcoa',
                                '--i-distance-matrix', DM,
                                '--o-pcoa', out_pcoa]
                            sh.write('echo "%s"\n' % ' '.join(cmd))
                            sh.write('%s\n' % ' '.join(cmd))
                            written += 1
            if written:
                xpbs_call(out_sh, out_pbs, '%s.pc.%s' % (prjct_nm, dataset), qiime_env,
                          '10', '1', '2', '2', 'gb')
                o.write('qsub %s\n' % out_pbs)
            else:
                print('\nNothing written in', out_sh, '--> removed')
                os.remove(out_sh)
    if written:
        print('[TO RUN] sh', run_pbs)
    return pcoas_d


def run_emperor(i_datasets_folder: str, pcoas_d: dict,
                prjct_nm: str, qiime_env: str) -> None:

    print('# Make EMPeror plots')
    job_folder = get_job_folder(i_datasets_folder, 'emperor')
    job_folder2 = get_job_folder(i_datasets_folder, 'emperor/chunks')
    odir = get_analysis_folder(i_datasets_folder, 'emperor')

    first_print = 0
    run_pbs = '%s/4_run_emperor.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dataset, meta_pcoas in pcoas_d.items():
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
                    out_sh = '%s/run_emperor_%s_%s.sh' % (job_folder2, dataset, basename(pcoa.split('_DM_')[0]))
                    out_pbs = '%s.pbs' % splitext(out_sh)[0]
                    with open(out_sh, 'w') as sh:
                        out_plot = '%s_emperor.qzv' % splitext(pcoa)[0].replace('/pcoa/', '/emperor/')
                        if not isdir(dirname(out_plot)):
                            os.makedirs(dirname(out_plot))
                        emp_cmd = [
                            'qiime', 'emperor', 'plot',
                            '--m-metadata-file', meta,
                            '--i-pcoa', pcoa,
                            '--o-visualization', out_plot]
                        sh.write('echo "%s"\n' % ' '.join(emp_cmd))
                        sh.write('%s\n' % ' '.join(emp_cmd))

                    xpbs_call(out_sh, out_pbs,
                              '%s.mprr%s-%s' % (prjct_nm, dataset, basename(pcoa.split('_DM_')[0])), qiime_env,
                              '10', '1', '1', '2', 'gb')
                    o.write('qsub %s\n' % out_pbs)
    print('[TO RUN] sh', run_pbs)
