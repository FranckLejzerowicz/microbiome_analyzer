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
from os.path import basename, isfile, isdir
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder, run_import


def run_deicode(i_folder: str, datasets: dict, p_perm_groups: str,
                force: bool, prjct_nm: str, qiime_env: str):

    print('# DEICODE (groups config in %s)' % p_perm_groups)
    def run_multi_deicode(
            rt, job_folder2, dat, tsv, tsv_pd, meta_pd,
            case, case_var, case_vals, force, prjct_nm, qiime_env):

        written = 0
        out_sh = '%s/run_beta_deicode_%s.sh' % (job_folder2, case)
        out_pbs = out_sh.replace('.sh', '.pbs')
        with open(out_sh, 'w') as sh:

            cur_rad = '/'.join([rt, basename(tsv).replace('.tsv', '_%s' % case)])
            new_tsv = '%s.tsv' % cur_rad
            new_meta = new_tsv.replace('/tab_', '/meta_')
            new_mat_qza = '%s_DM.qza' % cur_rad
            new_qza = '%s.qza' % cur_rad
            ordi_qza = '%s_deicode_ordination.qza' % cur_rad
            ordi_qzv = '%s_deicode_ordination_biplot.qzv' % cur_rad

            if force or not isfile(ordi_qzv):

                if 'ALL' in case:
                    new_meta_pd = meta_pd.copy()
                elif len([x for x in case_vals if '>' in x or '<' in x]):
                    new_meta_pd = meta_pd.copy()
                    for case_val in case_vals:
                        if case_val[0] == '>':
                            new_meta_pd = new_meta_pd[new_meta_pd[case_var] >= float(case_val[1:])].copy()
                        elif case_val[0] == '<':
                            new_meta_pd = new_meta_pd[new_meta_pd[case_var] <= float(case_val[1:])].copy()
                else:
                    new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')

                new_tsv_pd = tsv_pd[new_meta_pd.index.tolist()].copy()
                new_tsv_pd = new_tsv_pd.loc[new_tsv_pd.sum(1) > 0, :]
                new_tsv_pd = new_tsv_pd.loc[:, new_tsv_pd.sum(0) > 0]
                new_tsv_pd.reset_index().to_csv(new_tsv, index=False, sep='\t')

                cmd = run_import(new_tsv, new_qza, "FeatureTable[Frequency]")
                sh.write('%s\n' % ' '.join(cmd))

                cmd = [
                    'qiime', 'deicode', 'rpca',
                    '--i-table', new_qza,
                    '--p-min-feature-count', '10',
                    '--p-min-sample-count', '500',
                    '--o-biplot', ordi_qza,
                    '--o-distance-matrix', new_mat_qza
                ]
                sh.write('echo "%s"\n' % ' '.join(cmd))
                sh.write('%s\n' % ' '.join(cmd))
                cmd = [
                    'qiime', 'emperor', 'biplot',
                    '--i-biplot', ordi_qza,
                    '--m-sample-metadata-file', new_meta,
                    '--o-visualization', ordi_qzv,
                    '--p-number-of-features', '20'
                ]
                sh.write('echo "%s"\n' % ' '.join(cmd))
                sh.write('%s\n' % ' '.join(cmd))
                written += 1

        if written:
            xpbs_call(out_sh, out_pbs, '%s.dcd.%s' % (prjct_nm, case), qiime_env,
                      '2', '1', '1', '2', 'gb')
        else:
            os.remove(out_sh)

    job_folder = get_job_folder(i_folder, 'deicode')
    job_folder2 = get_job_folder(i_folder, 'deicode/chunks')

    with open(p_perm_groups) as handle:
        cases_dict = yaml.load(handle)
        # cases_dict = yaml.load(handle, Loader=yaml.FullLoader)
    cases_dict.update({'ALL': [[]]})

    jobs = []
    main_sh = '%s/3_run_beta_deicode.sh' % job_folder
    with open(main_sh, 'w') as o:

        for dat, tsvs_metas_list in datasets.items():
            rt = get_analysis_folder(i_folder, 'deicode/%s' % dat)

            for tsv_meta in tsvs_metas_list:
                tsv, meta = tsv_meta
                tsv_pd = pd.read_csv(tsv, header=0, index_col=0, sep='\t')
                meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')

                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        if len(case_vals):
                            case = '%s_%s' % (
                                case_var, '-'.join([x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
                        else:
                            case = case_var
                        out_pbs = '%s/run_beta_deicode_%s.pbs' % (job_folder2, case)
                        o.write('qsub %s\n' % out_pbs)
                        p = multiprocessing.Process(
                            target=run_multi_deicode,
                            args=(
                                rt, job_folder2, dat, tsv, tsv_pd, meta_pd,
                                case, case_var, case_vals, force, prjct_nm, qiime_env,
                            )
                        )
                        p.start()
                        jobs.append(p)
    for j in jobs:
        j.join()
    print('sh', main_sh)
