# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

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
    check_metadata_testing_groups
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_diversity_beta_group_significance,
    add_q2_types_to_meta, get_metric,
    check_absence_mat
)


def run_multi_perm(odir: str, beta_metrics: list, mat_qzas: list,
                   tsv: str, meta_pd: pd.DataFrame, cases_dict: dict,
                   testing_groups: tuple, dat: str, out_sh: str, out_pbs: str,
                   force: bool, prjct_nm: str, qiime_env: str) -> None:
    written = 0
    qza = '%s.qza' % splitext(tsv)[0]
    with open(out_sh, 'w') as cur_sh:
        for mat_qza in mat_qzas:
            metric = get_metric(beta_metrics, mat_qza)
            for case_var, case_vals_list in cases_dict.items():
                testing_groups_case_var = list(testing_groups) + [case_var]
                for case_vals in case_vals_list:
                    case_ = get_case(case_vals, metric, case_var)
                    for testing_group in testing_groups_case_var:
                        if testing_group == 'ALL':
                            continue
                        case = '%s__%s' % (case_, testing_group)
                        cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s' % case)
                        new_tsv = '%s.tsv' % cur_rad
                        new_meta = new_tsv.replace('/tab_', '/meta_')
                        new_qza = '%s.qza' % cur_rad
                        new_qzv = '%s_permanova.qzv' % cur_rad
                        new_mat_qza = '%s/%s.tsv' % (odir, basename(mat_qza).replace('.qza', '_%s.qza' % case))
                        if force or not isfile(new_qzv):
                            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                            if new_meta_pd[testing_group].unique().size > 1:
                                add_q2_types_to_meta(new_meta_pd, new_meta)
                                write_diversity_beta_group_significance(new_meta, mat_qza, new_mat_qza,
                                                                        qza, new_qza, testing_group,
                                                                        new_qzv, cur_sh)
                                written += 1
    run_xpbs(out_sh, out_pbs, '%s.perm.%s' % (prjct_nm, dat),
             qiime_env, '2', '1', '1', '1', 'gb', written, False, None)


def run_permanova(i_datasets_folder: str, datasets: dict, betas: dict,
                  main_testing_groups: tuple, p_perm_groups: str,
                  force: bool, prjct_nm: str, qiime_env: str):

    job_folder2 = get_job_folder(i_datasets_folder, 'permanova/chunks')
    beta_metrics = get_metrics('beta_metrics')

    main_cases_dict = get_main_cases_dict(p_perm_groups)

    jobs = []
    all_sh_pbs = []
    first_print = 0
    for dat, tsv_meta_pds in datasets.items():

        odir = get_analysis_folder(i_datasets_folder, 'permanova/%s' % dat)
        out_sh = '%s/run_beta_group_significance_%s.sh' % (job_folder2, dat)
        out_pbs = '%s.pbs' % splitext(out_sh)[0]
        all_sh_pbs.append((out_sh, out_pbs))

        tsv, meta = tsv_meta_pds
        meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')
        mat_qzas = betas[dat][meta]

        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict))
        testing_groups = check_metadata_testing_groups(meta, meta_pd, main_testing_groups)

        absence_mat = check_absence_mat(mat_qzas, first_print, 'PERMANOVA')
        if absence_mat:
            continue

        p = multiprocessing.Process(
            target=run_multi_perm,
            args=(odir, beta_metrics, mat_qzas, tsv, meta_pd, cases_dict,
                  testing_groups, dat, out_sh, out_pbs, force, prjct_nm, qiime_env))
        p.start()
        jobs.append(p)

    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_datasets_folder, 'permanova')
    main_sh = write_main_sh(job_folder, '3_run_beta_group_significance', all_sh_pbs)
    if main_sh:
        if p_perm_groups:
            print("# PERMANOVA (groups config in %s)" % p_perm_groups)
        else:
            print("# PERMANOVA")
        print('[TO RUN] sh', main_sh)