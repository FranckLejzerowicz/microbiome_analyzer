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

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_main_cases_dict,
    write_main_sh
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import (
    write_deicode_biplot,
    get_case, get_new_meta_pd
)


def run_multi_deicode(odir, tsv, meta_pd, cases_dict, dat, out_sh, out_pbs, force, prjct_nm, qiime_env):
    written = 0
    qza = '%s.qza' % splitext(tsv)[0]
    with open(out_sh, 'w') as cur_sh:
        for case_var, case_vals_list in cases_dict.items():
            for case_vals in case_vals_list:
                case = get_case(case_vals, '', case_var)
                cur_rad = '/'.join([odir, basename(tsv).replace('.tsv', '_%s' % case)])
                new_tsv = '%s.tsv' % cur_rad
                new_meta = new_tsv.replace('/tab_', '/meta_')
                new_mat_qza = '%s_DM.qza' % cur_rad
                new_qza = '%s.qza' % cur_rad
                ordi_qza = '%s_deicode_ordination.qza' % cur_rad
                ordi_qzv = '%s_deicode_ordination_biplot.qzv' % cur_rad
                if force or not isfile(ordi_qzv):
                    new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                    new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                    write_deicode_biplot(qza, new_meta, new_qza, ordi_qza, new_mat_qza, ordi_qzv, cur_sh)
                    written += 1
    run_xpbs(out_sh, out_pbs, '%s.dcd.%s' % (prjct_nm, dat),
             qiime_env, '2', '1', '1', '200', 'mb', written, False, None)


def run_deicode(i_datasets_folder: str, datasets: dict, p_perm_groups: str,
                force: bool, prjct_nm: str, qiime_env: str):

    job_folder2 = get_job_folder(i_datasets_folder, 'deicode/chunks')
    main_cases_dict = get_main_cases_dict(p_perm_groups)

    jobs = []
    all_sh_pbs = []
    for dat, tsv_meta_pds in datasets.items():

        odir = get_analysis_folder(i_datasets_folder, 'deicode/%s' % dat)
        out_sh = '%s/run_beta_deicode_%s.sh' % (job_folder2, dat)
        out_pbs = out_sh.replace('.sh', '.pbs')
        all_sh_pbs.append((out_sh, out_pbs))

        tsv, meta = tsv_meta_pds
        meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')

        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict))

        p = multiprocessing.Process(
            target=run_multi_deicode,
            args=(odir, tsv, meta_pd, cases_dict, dat, out_sh, out_pbs, force, prjct_nm, qiime_env))
        p.start()
        jobs.append(p)

    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_datasets_folder, 'deicode')
    main_sh = write_main_sh(job_folder, '3_run_beta_deicode', all_sh_pbs)
    if main_sh:
        if p_perm_groups:
            print('# DEICODE (groups config in %s)' % p_perm_groups)
        else:
            print('# DEICODE')
        print('sh', main_sh)
