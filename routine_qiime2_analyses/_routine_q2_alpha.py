# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import yaml
import pandas as pd
import multiprocessing
from os.path import basename, dirname, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics, get_job_folder, get_analysis_folder,
    run_import, run_export
)


def run_alpha(i_folder: str, datasets: dict, trees: dict,
              force: bool, prjct_nm: str, qiime_env: str) -> dict:

    alpha_metrics = get_metrics('alpha_metrics')

    job_folder = get_job_folder(i_folder, 'alpha')
    job_folder2 = get_job_folder(i_folder, 'alpha/chunks')

    written = 0
    diversities = {}
    run_pbs = '%s/1_run_alpha.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dataset, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            if dataset not in diversities:
                diversities[dataset] = {}
            out_sh = '%s/run_alpha_%s.sh' % (job_folder2, dataset)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                qza = '%s.qza' % splitext(tsv)[0]
                divs = [meta]
                for metric in alpha_metrics:
                    odir = get_analysis_folder(i_folder, 'alpha/%s' % dataset)
                    out_fp = '%s/%s_%s.qza' % (odir, basename(splitext(qza)[0]), metric)
                    out_tsv = '%s.tsv' % splitext(out_fp)[0]
                    if force or not isfile(out_fp):
                        if metric in ['faith_pd']:
                            if dataset not in trees:
                                continue
                            cmd = 'qiime diversity alpha-phylogenetic \\ \n'
                            cmd += '--i-table %s \\ \n' % trees[dataset][0]
                            cmd += '--i-phylogeny %s \\ \n' % trees[dataset][1]
                            cmd += '--p-metric %s \\ \n' % metric
                            cmd += '--o-alpha-diversity %s\n' % out_fp
                        else:
                            cmd = 'qiime diversity alpha \\ \n'
                            cmd += '--i-table %s \\ \n' % qza
                            cmd += '--p-metric %s \\ \n' % metric
                            cmd += '--o-alpha-diversity %s\n' % out_fp
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1
                    if force or not isfile(out_tsv):
                        cmd = run_export(out_fp, out_tsv, '')
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1
                    divs.append(out_fp)
                diversities[dataset][qza] = divs
            if written:
                xpbs_call(out_sh, out_pbs, '%s.mg.lph.%s' % (prjct_nm, dataset), qiime_env,
                          '4', '1', '1', '1', 'gb')
                o.write('qsub %s\n' % out_pbs)
            else:
                os.remove(out_sh)
    if written:
        print('# Calculate alpha diversity indices')
        print('[TO RUN] sh', run_pbs)
    return diversities


def merge_meta_alpha(i_folder: str, diversities: dict,
                     force: bool, prjct_nm: str, qiime_env: str):
    job_folder = get_job_folder(i_folder, 'alpha')
    job_folder2 = get_job_folder(i_folder, 'alpha/chunks')

    written = 0
    to_export = []
    run_pbs = '%s/2_run_merge_alphas.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dataset in diversities:
            out_sh = '%s/run_merge_alpha_%s.sh' % (job_folder2, dataset)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                for qza, divs in diversities[dataset].items():
                    meta = divs[0]
                    rad = splitext(qza)[0]
                    divs = divs[1:]
                    out_fp = '%s_alphas.qzv' % rad.replace('/data/', '/metadata/').replace('/tab_', '/meta_')
                    if force or not isfile(out_fp):
                        if not isdir(dirname(out_fp)):
                            os.makedirs(dirname(out_fp))
                        to_export.append(out_fp)
                        cmd = 'qiime metadata tabulate \\ \n'
                        cmd += '--o-visualization %s \\ \n' % out_fp
                        for div in divs:
                            cmd += '--m-input-file %s \\ \n' % div
                        cmd += '--m-input-file %s\n' % meta
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1
            if written:
                xpbs_call(out_sh, out_pbs, '%s.mg.mrg.lph.%s' % (prjct_nm, dataset), qiime_env,
                          '2', '1', '1', '150', 'mb')
                o.write('qsub %s\n' % out_pbs)
            else:
                os.remove(out_sh)
    if written:
        print('# Merge alpha diversity indices to metadata')
        print('[TO RUN] sh', run_pbs)
    return to_export


def export_meta_alpha(i_folder: str, to_export: list,
                      force: bool, prjct_nm: str, qiime_env: str):

    job_folder = get_job_folder(i_folder, 'alpha')
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
    if written:
        xpbs_call(out_sh, out_pbs, '%s.xprt.lph' % prjct_nm, qiime_env,
                  '2', '1', '1', '150', 'mb')
        print('# Export alpha diversity indices to metadata')
        print('[TO RUN] qsub', out_pbs)
    else:
        os.remove(out_sh)


def run_correlations(i_folder: str, datasets: dict, diversities: dict,
                     force: bool, prjct_nm: str, qiime_env: str):

    job_folder = get_job_folder(i_folder, 'alpha_correlations')
    job_folder2 = get_job_folder(i_folder, 'alpha_correlations/chunks')
    odir = get_analysis_folder(i_folder, 'alpha_correlations')
    written = 0
    run_pbs = '%s/4_run_alpha_correlation.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            out_sh = '%s/run_alpha_correlation_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                for method in ['spearman', 'pearson']:
                    for qza in diversities[dat]['%s.qza' % splitext(tsv)[0]][1:]:
                        out_fp = qza.replace('.qza', '_%s.qzv' % method).replace('/alpha/', '/alpha_correlations/')
                        if force or not isfile(out_fp):
                            if not isdir(dirname(out_fp)):
                                os.makedirs(dirname(out_fp))
                            cmd = 'qiime diversity alpha-correlation \\ \n'
                            cmd += '--i-alpha-diversity %s \\ \n' % qza
                            cmd += '--p-method %s \\ \n' % method
                            cmd += '--m-metadata-file %s \\ \n' % meta
                            cmd += '--o-visualization %s\n' % out_fp
                            sh.write('echo "%s"\n' % cmd)
                            sh.write('%s\n\n' % cmd)
                            written += 1
            if written:
                xpbs_call(out_sh, out_pbs, '%s.lphcrr.%s' % (prjct_nm, dat), qiime_env,
                          '10', '1', '1', '1', 'gb')
                o.write('qsub %s\n' % out_pbs)
            else:
                os.remove(out_sh)
    if written:
        print('# Correlate numeric metadata variables with alpha diversity indices')
        print('[TO RUN] sh', run_pbs)


def run_volatility(i_folder: str, datasets: dict, p_longi_column: str,
                   force: bool, prjct_nm: str, qiime_env: str) -> None:

    job_folder = get_job_folder(i_folder, 'longitudinal')
    job_folder2 = get_job_folder(i_folder, 'longitudinal/chunks')
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
            with open(out_sh, 'w') as sh:
                odir = get_analysis_folder(i_folder, 'longitudinal/%s' % dat)
                out_fp = '%s/%s_volatility.qzv' % (odir, dat)
                if force or not isfile(out_fp):
                    cmd = 'qiime longitudinal volatility \\ \n'
                    cmd += '--m-metadata-file %s \\ \n' % meta_alphas
                    cmd += '--p-state-column "%s" \\ \n' % time_point
                    cmd += '--p-individual-id-column "host_subject_id"'
                    cmd += '--o-visualization %s\n' % out_fp
                    sh.write('echo "%s"\n' % cmd)
                    sh.write('%s\n\n' % cmd)
                    written += 1
            if written:
                xpbs_call(out_sh, out_pbs, '%s.vltlt.%s' % (prjct_nm, dat), qiime_env,
                          '2', '1', '1', '100', 'mb')
                o.write('qsub %s\n' % out_pbs)
            else:
                os.remove(out_sh)
    if written:
        print('# Longitudinal change in alpha diversity indices')
        print('[TO RUN] sh', run_pbs)


def run_alpha_group_significance(i_folder: str, diversities: dict, p_perm_groups: str,
                                 force: bool, prjct_nm: str, qiime_env: str):

    alpha_metrics = get_metrics('alpha_metrics')

    def run_multi_kw(odir, meta_pd, div_qza, case, case_var, case_vals, cur_sh):

        cur_rad = odir + '/' + basename(div_qza).replace('.qza', '_%s' % case)
        new_qzv = '%s_kruskal-wallis.qzv' % cur_rad
        new_meta = '%s.meta' % cur_rad
        new_div = '%s.qza' % cur_rad

        with open(cur_sh, 'w') as cur_o:
            if force or not isfile(new_qzv):
                if 'ALL' in case:
                    new_meta_pd = meta_pd.copy()
                    cur_o.write('echo "cp %s %s"\n' % (div_qza, new_div))
                    cur_o.write('cp %s %s\n' % (div_qza, new_div))
                else:
                    new_tsv = '%s.tsv' % cur_rad
                    if len([x for x in case_vals if '>' in x or '<' in x]):
                        new_meta_pd = meta_pd.copy()
                        for case_val in case_vals:
                            if case_val[0] == '>':
                                new_meta_pd = new_meta_pd[new_meta_pd[case_var] >= float(case_val[1:])].copy()
                            elif case_val[0] == '<':
                                new_meta_pd = new_meta_pd[new_meta_pd[case_var] <= float(case_val[1:])].copy()
                    else:
                        new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()
                    new_tsv_pd = pd.read_csv(div_qza.replace('.qza', '.tsv'), header=0, sep='\t')
                    new_tsv_pd.rename(columns={new_tsv_pd.columns.tolist()[0]: 'Feature ID'}, inplace=True)
                    new_tsv_pd.set_index('Feature ID', inplace=True)
                    new_tsv_pd = new_tsv_pd.loc[new_meta_pd.index.tolist(),:]
                    new_tsv_pd.reset_index().to_csv(new_tsv, index=False, sep='\t')
                    cmd = run_import(new_tsv, new_div, 'SampleData[AlphaDiversity]')
                    cur_o.write('echo "%s"\n' % cmd)
                    cur_o.write('%s\n' % cmd)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')

                cmd = 'qiime diversity alpha-group-significance \\ \n'
                cmd += '--i-alpha-diversity %s \\ \n' % new_div
                cmd += '--m-metadata-file %s \\ \n' % new_meta
                cmd += '--o-visualization %s\n' % new_qzv
                cur_o.write('echo "%s"\n' % cmd)
                cur_o.write('%s\n\n' % cmd)

    job_folder = get_job_folder(i_folder, 'alpha_group_significance')
    job_folder2 = get_job_folder(i_folder, 'alpha_group_significance/chunks')

    if p_perm_groups:
        with open(p_perm_groups) as handle:
            # cases_dict = yaml.load(handle)
            cases_dict = yaml.load(handle, Loader=yaml.FullLoader)
    cases_dict.update({'ALL': [[]]})

    jobs = []
    first_print = 0
    main_sh = '%s/6_run_alpha_group_significance.sh' % job_folder
    with open(main_sh, 'w') as main_o:
        for dataset in diversities:
            odir = get_analysis_folder(i_folder, 'alpha_group_significance/%s' % dataset)
            out_sh = '%s/run_alpha_group_significance_%s.sh' % (job_folder2, dataset)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                for qza, divs in diversities[dataset].items():
                    meta = divs[0]
                    divs = divs[1:]
                    meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')
                    for div_qza in divs:
                        for metric in alpha_metrics:
                            if metric in div_qza:
                                break
                        if not isfile(div_qza):
                            if first_print:
                                print('Alpha diversity must be measured already to automatise Kruskal-Wallis tests\n'
                                      '\t(re-run this after step "1_run_alpha.sh" is done)')
                                first_print += 1
                            continue

                        for case_var, case_vals_list in cases_dict.items():
                            for case_vals in case_vals_list:
                                if len(case_vals):
                                    case = '%s_%s_%s' % (metric, case_var, '-'.join(
                                        [x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
                                else:
                                    case = '%s_%s' % (metric, case_var)
                                cur_sh = '%s/run_alpha_group_significance_%s_%s.sh' % (job_folder2, dataset, case)
                                sh.write('sh %s\n' % cur_sh)
                                p = multiprocessing.Process(
                                    target=run_multi_kw,
                                    args=(odir, meta_pd, div_qza, case,
                                          case_var, case_vals, cur_sh)
                                )
                                p.start()
                                jobs.append(p)
            xpbs_call(out_sh, out_pbs, '%s.perm.%s' % (prjct_nm, dataset), qiime_env,
                      '2', '1', '1', '2', 'gb')
            main_o.write('qsub %s\n' % out_pbs)
    for j in jobs:
        j.join()

    if p_perm_groups:
        print("# Kruskal-Wallis on alpha diversity (groups config in %s)" % p_perm_groups)
    else:
        print("# Kruskal-Wallis on alpha diversity")
    print('[TO RUN] sh', main_sh)
