# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import multiprocessing
from os.path import basename, dirname, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics, get_job_folder, get_analysis_folder,
    run_export, write_main_sh, get_main_cases_dict
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import (
    get_case, write_alpha_group_significance_cmd,
    get_new_meta_pd, get_new_alpha_div, write_metadata_tabulate,
    write_diversity_alpha, write_diversity_alpha_correlation,
    write_longitudinal_volatility
)


def run_alpha(i_datasets_folder: str, datasets: dict, datasets_phylo: dict, trees: dict,
              force: bool, prjct_nm: str, qiime_env: str) -> dict:
    """
    Computes the alpha diversity vectors for each dataset.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param trees: to be update with tree to use for a dataset phylogenetic analyses.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :return: {'dataset1': [ 'meta', {'div_index1': '.qza', 'div_index2': '.qza', ... }],
              'dataset2': [ 'meta', {'div_index1': '.qza', 'div_index2': '.qza', ... }], '...'}
    """
    alpha_metrics = get_metrics('alpha_metrics')
    job_folder = get_job_folder(i_datasets_folder, 'alpha')
    job_folder2 = get_job_folder(i_datasets_folder, 'alpha/chunks')
    written = 0
    diversities = {}
    run_pbs = '%s/1_run_alpha.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            if dat not in diversities:
                diversities[dat] = [meta, {}]
            out_sh = '%s/run_alpha_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                qza = '%s.qza' % splitext(tsv)[0]
                divs = []
                for metric in alpha_metrics:
                    odir = get_analysis_folder(i_datasets_folder, 'alpha/%s' % dat)
                    out_fp = '%s/%s_%s.qza' % (odir, basename(splitext(qza)[0]), metric)
                    out_tsv = '%s.tsv' % splitext(out_fp)[0]
                    if force or not isfile(out_fp):
                        ret_continue = write_diversity_alpha(out_fp, datasets_phylo, trees,
                                                             dat, qza, metric, cur_sh)
                        if ret_continue:
                            continue
                        cmd = run_export(out_fp, out_tsv, '')
                        cur_sh.write('echo "%s"\n' % cmd)
                        cur_sh.write('%s\n\n' % cmd)
                        written += 1
                    divs.append(out_fp)
                diversities[dat][1][qza] = divs
            run_xpbs(out_sh, out_pbs, '%s.mg.lph.%s' % (prjct_nm, dat),
                     qiime_env, '4', '1', '1', '1', 'gb', written, 'single', o)
    if written:
        print('# Calculate alpha diversity indices')
        print('[TO RUN] sh', run_pbs)
    return diversities


def merge_meta_alpha(i_datasets_folder: str, diversities: dict,
                     force: bool, prjct_nm: str, qiime_env: str):
    """
    Computes the alpha diversity vectors for each dataset.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param diversities: paths to [metadata, [alpha_divs]]
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :return:
    """
    job_folder = get_job_folder(i_datasets_folder, 'alpha')
    job_folder2 = get_job_folder(i_datasets_folder, 'alpha/chunks')
    written = 0
    to_export = []
    run_pbs = '%s/2_run_merge_alphas.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, (meta, qza_divs) in diversities.items():
            out_sh = '%s/run_merge_alpha_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for qza, divs in qza_divs.items():
                    rad = splitext(qza)[0]
                    out_fp = '%s_alphas.qzv' % rad.replace('/data/', '/metadata/').replace('/tab_', '/meta_')
                    if force or not isfile(out_fp):
                        if not isdir(dirname(out_fp)):
                            os.makedirs(dirname(out_fp))
                        to_export.append(out_fp)
                        write_metadata_tabulate(out_fp, divs, meta, cur_sh)
                        written += 1
            run_xpbs(out_sh, out_pbs, '%s.mrg.lph.%s' % (prjct_nm, dat),
                     qiime_env, '2', '1', '1', '150', 'mb', written, 'single', o)
    if written:
        print('# Merge alpha diversity indices to metadata')
        print('[TO RUN] sh', run_pbs)
    return to_export


def export_meta_alpha(i_datasets_folder: str, to_export: list,
                      force: bool, prjct_nm: str, qiime_env: str):
    """
    Export the alpha diversity vectors.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param to_export: files to export.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :return:
    """

    job_folder = get_job_folder(i_datasets_folder, 'alpha')
    out_sh = '%s/3_run_merge_alpha_export.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    written = 0
    with open(out_sh, 'w') as sh:
        for export in to_export:
            out_fp = '%s.tsv' % splitext(export)[0]
            if force or not isfile(out_fp):
                cmd = run_export(export, out_fp, '')
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n\n' % cmd)
                written += 1
    run_xpbs(out_sh, out_pbs, '%s.xprt.lph' % prjct_nm,
             qiime_env, '2', '1', '1', '150', 'mb', written,
             '# Export alpha diversity indices to metadata')


def run_correlations(i_datasets_folder: str, datasets: dict, diversities: dict,
                     force: bool, prjct_nm: str, qiime_env: str):
    """
    Run alpha-correlation: Alpha diversity correlation
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-correlation/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param diversities: alpha diversity qiime2 Arfetact per dataset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :return:
    """
    job_folder = get_job_folder(i_datasets_folder, 'alpha_correlations')
    job_folder2 = get_job_folder(i_datasets_folder, 'alpha_correlations/chunks')
    written = 0
    run_pbs = '%s/4_run_alpha_correlation.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            out_sh = '%s/run_alpha_correlation_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            odir = get_analysis_folder(i_datasets_folder, 'alpha_correlations/%s' % dat)
            with open(out_sh, 'w') as cur_sh:
                for method in ['spearman', 'pearson']:
                    for qza in diversities[dat][1]['%s.qza' % splitext(tsv)[0]]:
                        out_fp = qza.replace('.qza', '_%s.qzv' % method).replace('/alpha/', '/alpha_correlations/')
                        if force or not isfile(out_fp):
                            write_diversity_alpha_correlation(out_fp, qza, method, meta, cur_sh)
                            written += 1
            run_xpbs(out_sh, out_pbs, '%s.lphcrr.%s' % (prjct_nm, dat),
                     qiime_env, '10', '1', '1', '1', 'gb', written, 'single', o)
    if written:
        print('# Correlate numeric metadata variables with alpha diversity indices')
        print('[TO RUN] sh', run_pbs)


def run_volatility(i_datasets_folder: str, datasets: dict, p_longi_column: str,
                   force: bool, prjct_nm: str, qiime_env: str) -> None:
    """
    Run volatility: Generate interactive volatility plot.
    https://docs.qiime2.org/2019.10/plugins/available/longitudinal/volatility/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param p_longi_column: metadata column that is the time stratification.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :return:
    """
    job_folder = get_job_folder(i_datasets_folder, 'longitudinal')
    job_folder2 = get_job_folder(i_datasets_folder, 'longitudinal/chunks')
    written = 0
    first_print = 0
    first_print2 = 0
    run_pbs = '%s/5_run_volatility.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            meta_alphas = '%s_alphas.tsv' % splitext(meta)[0]
            if not isfile(meta_alphas):
                if not first_print:
                    print('\nWarning: First make sure you run alpha -> alpha merge -> alpha export'
                          ' before running volatility\n\t(if you need the alpha as a response variable)!')
                    first_print += 1
                continue
            with open(meta) as f:
                for line in f:
                    break
            time_point = [x for x in line.strip().split('\t') if p_longi_column in x][0]
            if not time_point:
                if not first_print2:
                    print('Variable %s not in metadata %s\n' % (p_longi_column, meta_alphas))
                    first_print2 += 1
                continue
            out_sh = '%s/run_volatility_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                odir = get_analysis_folder(i_datasets_folder, 'longitudinal/%s' % dat)
                out_fp = '%s/%s_volatility.qzv' % (odir, dat)
                if force or not isfile(out_fp):
                    write_longitudinal_volatility(out_fp, meta_alphas, time_point, cur_sh)
                    written += 1
            run_xpbs(out_sh, out_pbs, '%s.vltlt.%s' % (prjct_nm, dat),
                     qiime_env, '2', '1', '1', '100', 'mb', written, 'single', o)
    if written:
        print('# Longitudinal change in alpha diversity indices')
        print('[TO RUN] sh', run_pbs)


def run_multi_kw(odir: str, meta_pd: pd.DataFrame, diversities: dict,
                 alpha_metrics: list, cases_dict: dict, out_sh: str,
                 out_pbs: str, dat: str, force: bool, prjct_nm: str,
                 qiime_env: str) -> None:
    """
    Run alpha-group-significance: Alpha diversity comparisons.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-group-significance/
    (in-loop function).

    :param odir: output analysis directory.
    :param meta_pd: metadata table.
    :param diversities: alpha diversity qiime2 Arfetact per dataset.
    :param alpha_metrics: list of alpha diversity metrics.
    :param cases_dict: groups to test.
    :param out_sh: input bash script file.
    :param out_pbs: output torque script file.
    :param dat: current dataset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :return:
    """
    written = 0
    with open(out_sh, 'w') as cur_sh:
        for qza, divs in diversities[dat][1].items():
            for div_qza in divs:
                for metric in alpha_metrics:
                    if metric in div_qza:
                        break
                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        case = get_case(case_vals, metric, case_var)
                        cur_rad = odir + '/' + basename(div_qza).replace('.qza', '_%s' % case)
                        new_qzv = '%s_kruskal-wallis.qzv' % cur_rad
                        if force or not isfile(new_qzv):
                            new_meta = '%s.meta' % cur_rad
                            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                            new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                            new_div = get_new_alpha_div(case, div_qza, cur_rad, new_meta_pd, cur_sh)
                            write_alpha_group_significance_cmd(new_div, new_meta, new_qzv, cur_sh)
                            written += 1
    run_xpbs(out_sh, out_pbs, '%s.kv.%s' % (prjct_nm, dat),
             qiime_env, '2', '1', '1', '1', 'gb', written, '', None)


def run_alpha_group_significance(i_datasets_folder: str, diversities: dict, p_perm_groups: str,
                                 force: bool, prjct_nm: str, qiime_env: str) -> None:
    """
    Run alpha-group-significance: Alpha diversity comparisons.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-group-significance/
    Main per-dataset looper for the Kruskal-Wallis tests on alpha diversity vectors.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param diversities: alpha diversity qiime2 Arfetact per dataset.
    :param p_perm_groups: path to the subsets file.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    """
    alpha_metrics = get_metrics('alpha_metrics')
    main_cases_dict = get_main_cases_dict(p_perm_groups)

    jobs = []
    all_sh_pbs = []
    first_print = 0
    job_folder2 = get_job_folder(i_datasets_folder, 'alpha_group_significance/chunks')
    for dat in diversities:
        odir = get_analysis_folder(i_datasets_folder, 'alpha_group_significance/%s' % dat)
        out_sh = '%s/run_alpha_group_significance_%s.sh' % (job_folder2, dat)
        out_pbs = '%s.pbs' % splitext(out_sh)[0]
        all_sh_pbs.append((out_sh, out_pbs))

        presence_mat = [1 for (qza, divs) in diversities[dat][1] for div in divs if isfile(div)]
        if not presence_mat:
            if not first_print:
                print('Alpha diversity must be measured already to automatise Kruskal-Wallis tests\n'
                      '\t(re-run this after step "1_run_alpha.sh" is done)')
                first_print += 1
            continue

        meta = diversities[dat][0]
        meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')
        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict))

        p = multiprocessing.Process(
            target=run_multi_kw,
            args=(odir, meta_pd, diversities, alpha_metrics, cases_dict,
                  out_sh, out_pbs, dat, force, prjct_nm, qiime_env))
        p.start()
        jobs.append(p)

    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_datasets_folder, 'alpha_group_significance')
    main_sh = write_main_sh(job_folder, '6_run_alpha_group_significance', all_sh_pbs)
    if main_sh:
        if p_perm_groups:
            print("# Kruskal-Wallis on alpha diversity (groups config in %s)" % p_perm_groups)
        else:
            print("# Kruskal-Wallis on alpha diversity")
        print('[TO RUN] sh', main_sh)
