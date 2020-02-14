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
from os.path import basename, isfile, splitext
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder, run_import


def run_deicode(i_folder: str, datasets: dict, p_perm_groups: str,
                force: bool, prjct_nm: str, qiime_env: str):

    def run_multi_deicode(rt, out_sh, out_pbs, tsv, meta_pd,
                          case, case_var, case_vals,
                          force, prjct_nm, qiime_env):
        written = 0
        with open(out_sh, 'w') as sh:
            qza = '%s.qza' % splitext(tsv)[0]
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

                cmd = 'qiime feature-table filter-samples \\ \n'
                cmd += '--i-table %s \\ \n' % qza
                cmd += '--m-metadata-file %s \\ \n' % new_meta
                cmd += '--o-filtered-table %s \\ \n' % new_qza
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n' % cmd)

                cmd = 'qiime deicode rpca \\ \n'
                cmd += '--i-table %s \\ \n' % new_qza
                cmd += '--p-min-feature-count 10 \\ \n'
                cmd += '--p-min-sample-count 500 \\ \n'
                cmd += '--o-biplot %s \\ \n' % ordi_qza
                cmd += '--o-distance-matrix %s\n' % new_mat_qza
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n' % cmd)

                cmd = 'qiime emperor biplot'
                cmd += '--i-biplot %s \\ \n' % ordi_qza
                cmd += '--m-sample-metadata-file %s \\ \n' % new_meta
                cmd += '--o-visualization %s \\ \n' % ordi_qzv
                cmd += '--p-number-of-features 20\n'
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n' % cmd)
                written += 1

        if written:
            xpbs_call(out_sh, out_pbs, '%s.dcd.%s' % (prjct_nm, case), qiime_env,
                      '2', '1', '1', '2', 'gb')
        else:
            os.remove(out_sh)

    job_folder = get_job_folder(i_folder, 'deicode')
    job_folder2 = get_job_folder(i_folder, 'deicode/chunks')

    if p_perm_groups:
        with open(p_perm_groups) as handle:
            # cases_dict = yaml.load(handle)
            cases_dict = yaml.load(handle, Loader=yaml.FullLoader)
    cases_dict.update({'ALL': [[]]})

    jobs = []
    all_shs = []
    main_sh = '%s/3_run_beta_deicode.sh' % job_folder
    with open(main_sh, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            rt = get_analysis_folder(i_folder, 'deicode/%s' % dat)
            meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')
            for case_var, case_vals_list in cases_dict.items():
                for case_vals in case_vals_list:
                    if len(case_vals):
                        case = '%s_%s' % (
                            case_var, '-'.join([x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
                    else:
                        case = case_var
                    out_sh = '%s/run_beta_deicode_%s.sh' % (job_folder2, case)
                    all_shs.append(out_sh)
                    out_pbs = out_sh.replace('.sh', '.pbs')
                    o.write('qsub %s\n' % out_pbs)
                    p = multiprocessing.Process(
                        target=run_multi_deicode,
                        args=(
                            rt, out_sh, out_pbs, tsv, meta_pd,
                            case, case_var, case_vals,
                            force, prjct_nm, qiime_env,
                        )
                    )
                    p.start()
                    jobs.append(p)
    for j in jobs:
        j.join()

    if len([1 for sh in all_shs for line in open(sh).readlines() if len(line.strip())]):
        if p_perm_groups:
            print('# DEICODE (groups config in %s)' % p_perm_groups)
        else:
            print('# DEICODE')
        print('sh', main_sh)
