# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import yaml
import pandas as pd
from os.path import basename, isfile, splitext
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics, get_job_folder,
    get_analysis_folder,
    get_main_cases_dict,
    write_main_sh
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict,
    check_metadata_formulas,
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_diversity_adonis,
    get_metric, check_absence_mat
)


def run_multi_adonis(odir: str, beta_metrics: list, mat_qzas: list,
                     tsv: str, meta_pd: pd.DataFrame, cases_dict: dict,
                     formulas: dict, dat: str, out_sh: str, out_pbs: str,
                     force: bool, prjct_nm: str, qiime_env: str) -> None:
    written = 0
    qza = '%s.qza' % splitext(tsv)[0]
    with open(out_sh, 'w') as cur_sh:
        for mat_qza in mat_qzas:
            metric = get_metric(beta_metrics, mat_qza)
            for form, formula in formulas.items():
                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        case = get_case(case_vals, metric, case_var, form)
                        cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s' % case)
                        new_tsv = '%s.tsv' % cur_rad
                        new_meta = new_tsv.replace('/tab_', '/meta_')
                        new_qza = '%s.qza' % cur_rad
                        new_qzv = '%s_adonis.qzv' % cur_rad
                        new_mat_qza = '%s/%s.tsv' % (odir, basename(mat_qza).replace('.qza', '_%s.qza' % case))
                        if force or not isfile(new_qzv):
                            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                            new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                            write_diversity_adonis(new_meta, mat_qza, new_mat_qza, qza,
                                                   new_qza, formula, new_qzv, cur_sh)
                            written += 1
    run_xpbs(out_sh, out_pbs, '%s.dns.%s' % (prjct_nm, dat),
             qiime_env, '2', '1', '1', '1', 'gb', written, False, None)


def run_adonis(p_formulas: str, i_datasets_folder: str, datasets: dict, betas: dict,
               p_perm_groups: str, force: bool, prjct_nm: str, qiime_env: str):

    job_folder2 = get_job_folder(i_datasets_folder, 'adonis/chunks')
    beta_metrics = get_metrics('beta_metrics')

    main_cases_dict = get_main_cases_dict(p_perm_groups)

    with open(p_formulas) as handle:
        formulas = yaml.load(handle, Loader=yaml.FullLoader)

    jobs = []
    all_sh_pbs = []
    first_print = 0
    for dat, tsv_meta_pds in datasets.items():

        odir = get_analysis_folder(i_datasets_folder, 'adonis/%s' % dat)
        out_sh = '%s/run_adonis_%s.sh' % (job_folder2, dat)
        out_pbs = '%s.pbs' % splitext(out_sh)[0]
        all_sh_pbs.append((out_sh, out_pbs))

        tsv, meta = tsv_meta_pds
        meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')
        mat_qzas = betas[dat][meta]

        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict))
        formulas = check_metadata_formulas(meta, meta_pd, dict(formulas))

        absence_mat = check_absence_mat(mat_qzas, first_print, 'ADONIS')
        if absence_mat:
            continue

        p = multiprocessing.Process(
            target=run_multi_adonis,
            args=(odir, beta_metrics, mat_qzas, tsv, meta_pd, cases_dict,
                  formulas, dat, out_sh, out_pbs, force, prjct_nm, qiime_env))
        p.start()
        jobs.append(p)

    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_datasets_folder, 'adonis')
    main_sh = write_main_sh(job_folder, '3_run_adonis', all_sh_pbs)
    if main_sh:
        if p_perm_groups:
            print("# Run Adonis (groups config in %s)" % p_perm_groups)
        else:
            print("# Run Adonis")
        print('[TO RUN] sh', main_sh)
