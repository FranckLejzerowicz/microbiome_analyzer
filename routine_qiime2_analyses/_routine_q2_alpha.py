# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, sys
import pandas as pd
import multiprocessing
from os.path import basename, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics, get_job_folder, get_analysis_folder,
    write_main_sh, get_main_cases_dict, read_meta_pd,
    get_alpha_subsets, get_raref_tab_meta_pds
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import (
    get_case, write_alpha_group_significance_cmd,
    get_new_meta_pd, get_new_alpha_div, write_metadata_tabulate,
    write_diversity_alpha, write_diversity_alpha_correlation,
    write_longitudinal_volatility, get_metric, get_subset,
    write_filter_features, run_export
)


def run_alpha(i_datasets_folder: str, datasets: dict, datasets_read: dict,
              datasets_phylo: dict, p_alpha_subsets: str, trees: dict,
              force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> dict:
    """
    Computes the alpha diversity vectors for each dataset.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param p_alpha_subsets: Subsets for alpha diversity.
    :param trees: to be update with tree to use for a dataset phylogenetic analyses.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: {'dataset1': [ 'meta', {'div_index1': '.qza', 'div_index2': '.qza', ... }],
              'dataset2': [ 'meta', {'div_index1': '.qza', 'div_index2': '.qza', ... }], '...'}
    """
    # alpha_subsets_deletions = []
    alpha_metrics = get_metrics('alpha_metrics')
    alpha_subsets = get_alpha_subsets(p_alpha_subsets)
    job_folder = get_job_folder(i_datasets_folder, 'alpha')
    job_folder2 = get_job_folder(i_datasets_folder, 'alpha/chunks')
    written = 0
    diversities = {}
    run_pbs = '%s/1_run_alpha.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds

            if datasets_read[dat] == 'raref':
                if not isfile(tsv):
                    print('Must have run rarefcation to use it further...\nExiting')
                    sys.exit(0)
                tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                datasets_read[dat] = [tsv_pd, meta_pd]
            else:
                tsv_pd, meta_pd = datasets_read[dat]

            out_sh = '%s/run_alpha_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                qza = '%s.qza' % splitext(tsv)[0]
                divs = {}
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
                    divs.setdefault('', []).append(out_fp)

                if alpha_subsets and dat in alpha_subsets:
                    for subset, subset_regex in alpha_subsets[dat].items():
                        odir = get_analysis_folder(i_datasets_folder, 'alpha/%s/%s' % (dat, subset))
                        qza_subset = '%s/%s_%s.qza' % (odir, basename(splitext(qza)[0]),  subset)
                        meta_subset = '%s.meta' % splitext(qza_subset)[0]
                        nfeats = get_subset(tsv_pd, subset, meta_subset, subset_regex)
                        if not nfeats:
                            continue
                        write_filter_features(qza, qza_subset, meta_subset, cur_sh)

                        for metric in alpha_metrics:
                            out_fp = '%s/%s_%s__%s.qza' % (odir, basename(splitext(qza)[0]), metric, subset)
                            out_tsv = '%s.tsv' % splitext(out_fp)[0]
                            # alpha_subsets_deletions.extend([out_fp, out_tsv, meta_subset])
                            if force or not isfile(out_fp):
                                ret_continue = write_diversity_alpha(out_fp, datasets_phylo, trees,
                                                                     dat, qza_subset, metric, cur_sh)
                                if ret_continue:
                                    continue
                                cmd = run_export(out_fp, out_tsv, '')
                                cur_sh.write('echo "%s"\n' % cmd)
                                cur_sh.write('%s\n\n' % cmd)
                                written += 1
                            divs.setdefault(subset, []).append(out_fp)

                diversities[dat] = divs
            run_xpbs(out_sh, out_pbs, '%s.mg.lph.%s' % (prjct_nm, dat),
                     qiime_env, '4', '1', '1', '1', 'gb',
                     chmod, written, 'single', o)
    if written:
        print_message('# Calculate alpha diversity indices', 'sh', run_pbs)
    return diversities


def merge_meta_alpha(i_datasets_folder: str, datasets: dict, diversities: dict,
                     force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> dict:
    """
    Computes the alpha diversity vectors for each dataset.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param diversities: paths to [alpha_divs]
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    :return:
    """
    job_folder = get_job_folder(i_datasets_folder, 'tabulate')
    job_folder2 = get_job_folder(i_datasets_folder, 'tabulate/chunks')
    written = 0
    to_export = {}
    run_pbs = '%s/2_run_merge_alphas.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, group_divs in diversities.items():
            tsv, meta = datasets[dat]
            base = basename(splitext(tsv)[0]).lstrip('tab_')
            out_sh = '%s/run_merge_alpha_%s.sh' % (job_folder2, base)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for group, divs in group_divs.items():
                    if group:
                        output_folder = get_analysis_folder(i_datasets_folder, 'tabulate/%s/%s' % (dat, group))
                    else:
                        output_folder = get_analysis_folder(i_datasets_folder, 'tabulate/%s' % dat)
                    out_fp = '%s/%s_alphas__%s.qzv' % (output_folder, base, group)
                    out_fp_tsv = '%s.tsv' % splitext(out_fp)[0]
                    to_export.setdefault(dat, []).append(out_fp_tsv)
                    if force or not isfile(out_fp):
                        write_metadata_tabulate(out_fp, divs, meta, cur_sh)
                        cmd = run_export(out_fp, out_fp_tsv, '')
                        cur_sh.write('echo "%s"\n' % cmd)
                        cur_sh.write('%s\n\n' % cmd)
                        written += 1
            run_xpbs(out_sh, out_pbs, '%s.mrg.lph.%s' % (prjct_nm, base),
                     qiime_env, '2', '1', '1', '150', 'mb',
                     chmod, written, 'single', o)
    if written:
        print_message('# Merge and export alpha diversity indices', 'sh', run_pbs)
    return to_export


def export_meta_alpha(datasets: dict, to_export: dict) -> None:
    """
    Export the alpha diversity vectors.

    :param datasets: list of datasets.
    :param to_export: files to export per dataset.
    """
    first_print = True
    for dat, meta_alphas_fps in to_export.items():
        tsv, meta = datasets[dat]
        meta_alphas_fps_exist = [x for x in meta_alphas_fps if isfile(x)]
        if len(meta_alphas_fps_exist) != len(meta_alphas_fps):
            if not first_print:
                print('\nWarning: First make sure you run alpha -> alpha merge/export (2_run_merge_alphas.sh) '
                      ' before running volatility\n\t(if you need the alpha as a response variable)!')
                first_print = False
            continue

        meta_alphas_pds = []
        for meta_alpha_fp in meta_alphas_fps_exist:
            meta_alpha_pd = pd.read_csv(meta_alpha_fp, header=0, dtype=str, sep='\t')
            meta_alpha_pd.rename(columns={meta_alpha_pd.columns[0]: 'sample_name'}, inplace=True)
            meta_alpha_pd.set_index('sample_name', inplace=True)

            group = meta_alpha_fp.split('_alphas__')[-1].split('.tsv')[0]
            if group != '':
                replace_cols = dict((x, '%s_%s' % (group, x)) for x in meta_alpha_pd.columns)
                meta_alpha_pd.rename(columns=replace_cols, inplace=True)
            meta_alphas_pds.append(meta_alpha_pd)
        meta_alphas_pd = pd.concat(meta_alphas_pds, axis=1)
        if meta_alphas_pd.index.tolist()[0] == '#q2:types':
            meta_alphas_pd = meta_alphas_pd.iloc[1:,:]
        meta_alphas_pd = meta_alphas_pd.reset_index()
        meta_alphas_pd.rename(columns={meta_alphas_pd.columns[0]: 'sample_name'}, inplace=True)
        meta_pd = read_meta_pd(meta)

        meta_alphas_pd = meta_pd.merge(meta_alphas_pd, on='sample_name', how='left')
        meta_alpha_fpo = '%s_alphas.tsv' % splitext(meta)[0]
        meta_alphas_pd.to_csv(meta_alpha_fpo, index=False, sep='\t')
        if os.getcwd().startswith('/panfs'):
            meta_alpha_fpo = meta_alpha_fpo.replace(os.getcwd(), '')
        print(' -> Written:', meta_alpha_fpo)


def run_correlations(i_datasets_folder: str, datasets: dict, diversities: dict,
                     force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    Run alpha-correlation: Alpha diversity correlation
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-correlation/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param diversities: alpha diversity qiime2 Arfetact per dataset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'alpha_correlations')
    job_folder2 = get_job_folder(i_datasets_folder, 'alpha_correlations/chunks')
    written = 0
    run_pbs = '%s/4_run_alpha_correlation.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            if dat not in diversities:
                continue
            tsv, meta = tsv_meta_pds
            out_sh = '%s/run_alpha_correlation_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            odir = get_analysis_folder(i_datasets_folder, 'alpha_correlations/%s' % dat)
            with open(out_sh, 'w') as cur_sh:
                # for method in ['spearman', 'pearson']:
                for method in ['spearman']:
                    for qza in diversities[dat]['']:
                        out_fp = qza.replace('.qza', '_%s.qzv' % method).replace('/alpha/', '/alpha_correlations/')
                        if force or not isfile(out_fp):
                            write_diversity_alpha_correlation(out_fp, qza, method, meta, cur_sh)
                            written += 1
            run_xpbs(out_sh, out_pbs, '%s.lphcrr.%s' % (prjct_nm, dat),
                     qiime_env, '10', '1', '1', '1', 'gb',
                     chmod, written, 'single', o)
    if written:
        print_message('# Correlate numeric metadata variables with alpha diversity indices', 'sh', run_pbs)


def run_volatility(i_datasets_folder: str, datasets: dict, p_longi_column: str,
                   force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    Run volatility: Generate interactive volatility plot.
    https://docs.qiime2.org/2019.10/plugins/available/longitudinal/volatility/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param p_longi_column: metadata column that is the time stratification.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
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
                    print('\nWarning: First make sure you run alpha -> alpha merge/export (2_run_merge_alphas.sh) '
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
                     qiime_env, '2', '1', '1', '100', 'mb',
                     chmod, written, 'single', o)
    if written:
        print_message('# Longitudinal change in alpha diversity indices', 'sh', run_pbs)


def run_multi_kw(odir: str, meta_pd: pd.DataFrame, div_qza: str, case_vals_list: list,
                 metric: str, case_var: str, cur_sh: str, force: bool) -> None:
    """
    Run alpha-group-significance: Alpha diversity comparisons.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-group-significance/
    (in-loop function).

    :param odir: output analysis directory.
    :param meta_pd: metadata table.
    :param div_qza:
    :param case_vals_list:
    :param metric:
    :param case_var:
    :param cur_sh: input bash script file.
    :param force: Force the re-writing of scripts for all commands.
    """
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        for case_vals in case_vals_list:
            case = get_case(case_vals, metric, case_var)
            cur_rad = odir + '/' + basename(div_qza).replace('.qza', '_%s' % case)
            new_qzv = '%s_kruskal-wallis.qzv' % cur_rad
            if force or not isfile(new_qzv):
                new_meta = '%s.meta' % cur_rad
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                new_div = get_new_alpha_div(case, div_qza, cur_rad, new_meta_pd, cur_sh_o)
                write_alpha_group_significance_cmd(new_div, new_meta, new_qzv, cur_sh_o)
                remove = False
    if remove and isfile(cur_sh):
        os.remove(cur_sh)


def run_alpha_group_significance(i_datasets_folder: str, datasets: dict, diversities: dict,
                                 p_perm_groups: str, force: bool, prjct_nm: str,
                                 qiime_env: str, chmod: str) -> None:
    """
    Run alpha-group-significance: Alpha diversity comparisons.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-group-significance/
    Main per-dataset looper for the Kruskal-Wallis tests on alpha diversity vectors.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param diversities: alpha diversity qiime2 Arfetact per dataset.
    :param p_perm_groups: path to the subsets file.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """

    job_folder2 = get_job_folder(i_datasets_folder, 'alpha_group_significance/chunks')
    alpha_metrics = get_metrics('alpha_metrics')
    main_cases_dict = get_main_cases_dict(p_perm_groups)

    jobs = []
    all_sh_pbs = {}
    first_print = 0
    for dat in diversities:
        presence_mat = [1 for qza in diversities[dat][''] if isfile(qza)]
        if not presence_mat:
            if not first_print:
                print('Alpha diversity must be measured already to automatise Kruskal-Wallis tests\n'
                      '\t(re-run this after step "1_run_alpha.sh" is done)')
                first_print += 1
            continue

        meta = datasets[dat][1]
        meta_pd = read_meta_pd(meta)
        meta_pd = meta_pd.set_index('sample_name')
        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'alpha Kruskal-Wallis')

        odir = get_analysis_folder(i_datasets_folder, 'alpha_group_significance/%s' % dat)
        for qza in diversities[dat]['']:
            metric = get_metric(alpha_metrics, qza)
            div_tsv = '%s.tsv' % splitext(qza)[0]
            if not isfile(div_tsv) or not isfile(div_tsv):
                print('  [KRUSKAL-WALLIS] metric %s not calculated\nSkipping it...' % metric)
                continue
            out_sh = '%s/run_alpha_group_significance_%s_%s.sh' % (job_folder2, dat, metric)
            for case_var, case_vals_list in cases_dict.items():
                cur_sh = '%s/run_adonis_%s_%s_%s.sh' % (
                    job_folder2, dat, metric, case_var)
                cur_sh = cur_sh.replace(' ', '-')
                all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                p = multiprocessing.Process(
                    target=run_multi_kw,
                    args=(odir, meta_pd, qza, case_vals_list,
                          metric, case_var, cur_sh, force))
                p.start()
                jobs.append(p)
    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_datasets_folder, 'alpha_group_significance')
    main_sh = write_main_sh(job_folder, '6_run_alpha_group_significance', all_sh_pbs,
                            '%s.kv' % prjct_nm, '2', '1', '1', '1', 'gb',
                            qiime_env, chmod)
    if main_sh:
        if p_perm_groups:
            print("# Kruskal-Wallis on alpha diversity (groups config in %s)" % p_perm_groups)
        else:
            print("# Kruskal-Wallis on alpha diversity")
        print_message('', 'sh', main_sh)
