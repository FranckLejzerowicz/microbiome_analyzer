# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
from os.path import basename, isfile, splitext
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics, get_job_folder,
    get_analysis_folder,
    get_main_cases_dict,
    get_formulas_dict,
    write_main_sh,
    read_meta_pd
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


def run_multi_adonis(odir: str, case_vals_list: list, metric: str, case_var: str,
                     form: str, formula: str, tsv: str, mat_qza: str,
                     meta_pd: pd.DataFrame, cur_sh: str, force: bool) -> None:
    """
    Run adonis: adonis PERMANOVA test for beta group significance.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/adonis/
    (in-loop function).

    :param odir: output analysis directory.
    :param case_vals_list:
    :param metric:
    :param case_var:
    :param form:
    :param formula:
    :param tsv: features table input to the beta diversity matrix.
    :param mat_qza:
    :param meta_pd: metadata table.
    :param cur_sh: input bash script file.
    :param force: Force the re-writing of scripts for all commands.
    """
    remove = True
    qza = '%s.qza' % splitext(tsv)[0]
    with open(cur_sh, 'w') as cur_sh_o:
        for case_vals in case_vals_list:
            case = '%s__%s' % (metric, get_case(case_vals, case_var, form))
            cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s_%s' % (metric, case))
            new_meta = '%s.meta' % cur_rad
            new_qza = '%s.qza' % cur_rad
            new_qzv = '%s_adonis.qzv' % cur_rad
            new_mat_qza = '%s/%s' % (odir, basename(mat_qza).replace('.qza', '_%s.qza' % case))
            if force or not isfile(new_qzv):
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                write_diversity_adonis(new_meta, mat_qza, new_mat_qza, qza,
                                       new_qza, formula, new_qzv, cur_sh_o)
                remove = False
    if remove and isfile(cur_sh):
        os.remove(cur_sh)


def run_adonis(p_formulas: str, i_data_sets_folder: str, data_sets: dict, betas: dict,
               p_perm_groups: str, force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    Run beta-group-significance: Beta diversity group significance.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-group-significance/
    Main per-dataset looper for the ADONIS tests on beta diversity matrices.

    :param p_formulas: formulas to test.
    :param i_data_sets_folder: Path to the folder containing the data/metadata subfolders.
    :param data_sets: list of datasets.
    :param betas: beta diversity matrices.
    :param p_perm_groups: groups to subset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """

    job_folder2 = get_job_folder(i_data_sets_folder, 'adonis/chunks')
    beta_metrics = get_metrics('beta_metrics')

    main_cases_dict = get_main_cases_dict(p_perm_groups)
    formulas = get_formulas_dict(p_formulas)

    jobs = []
    all_sh_pbs = {}
    first_print = 0
    for dat, tsv_meta_pds in data_sets.items():

        tsv, meta = tsv_meta_pds
        mat_qzas = betas[dat][meta]

        absence_mat = check_absence_mat(mat_qzas, first_print, 'ADONIS')
        if absence_mat:
            continue

        meta_pd = read_meta_pd(meta)
        meta_pd = meta_pd.set_index('sample_name')
        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'ADONIS')
        formulas = check_metadata_formulas(meta, meta_pd, dict(formulas), 'ADONIS')

        odir = get_analysis_folder(i_data_sets_folder, 'adonis/%s' % dat)
        for mat_qza in mat_qzas:
            metric = get_metric(beta_metrics, mat_qza)
            out_sh = '%s/run_adonis_%s_%s.sh' % (job_folder2, dat, metric)
            for form, formula in formulas.items():
                for case_var, case_vals_list in cases_dict.items():
                    cur_sh = '%s/run_adonis_%s_%s_%s_%s.sh' % (
                        job_folder2, dat, metric, form, case_var)
                    cur_sh = cur_sh.replace(' ', '-')
                    all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                    p = multiprocessing.Process(
                        target=run_multi_adonis,
                        args=(odir, case_vals_list, metric, case_var, form,
                              formula, tsv, mat_qza, meta_pd, out_sh, force))
                    p.start()
                    jobs.append(p)
    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_data_sets_folder, 'adonis')
    main_sh = write_main_sh(job_folder, '3_run_adonis', all_sh_pbs,
                            '%s.dns' % prjct_nm, '2', '1', '1', '1', 'gb',
                            qiime_env, chmod)
    if main_sh:
        if p_perm_groups:
            print("# Run Adonis (groups config in %s)" % p_perm_groups)
        else:
            print("# Run Adonis")
        print_message('', 'sh', main_sh)
