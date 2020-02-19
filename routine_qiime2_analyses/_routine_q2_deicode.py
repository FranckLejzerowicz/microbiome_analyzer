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


def run_multi_deicode(odir: str, tsv: str, meta_pd: pd.DataFrame, case_var: str,
                      case_vals_list: list, cur_sh: str, force: bool) -> None:
    """
    Performs robust center log-ratio transform robust PCA and
    ranks the features by the loadings of the resulting SVD.
    https://library.qiime2.org/plugins/deicode/19/
    (in-loop function).

    :param odir: output analysis directory.
    :param tsv: features table input to the beta diversity matrix.
    :param meta_pd: metadata table.
    :param case_var: metadata variable to make groups from.
    :param case_vals_list: groups for the metadata variable.
    :param cur_sh: input bash script file.
    :param force: Force the re-writing of scripts for all commands.
    """

    remove = True
    qza = '%s.qza' % splitext(tsv)[0]
    with open(cur_sh, 'w') as cur_sh_o:
        for case_vals in case_vals_list:
            case = get_case(case_vals, '', case_var)
            cur_rad = '/'.join([odir, basename(tsv).replace('.tsv', '_%s' % case)])
            new_meta = '%s.meta' % cur_rad
            new_mat_qza = '%s_DM.qza' % cur_rad
            new_qza = '%s.qza' % cur_rad
            ordi_qza = '%s_deicode_ordination.qza' % cur_rad
            ordi_qzv = '%s_deicode_ordination_biplot.qzv' % cur_rad
            if force or not isfile(ordi_qzv):
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                write_deicode_biplot(qza, new_meta, new_qza, ordi_qza,
                                     new_mat_qza, ordi_qzv, cur_sh_o)
                remove = False
    if remove and isfile(cur_sh):
        os.remove(cur_sh)


def run_deicode(i_data_sets_folder: str, data_sets: dict, p_perm_groups: str,
                force: bool, prjct_nm: str, qiime_env: str) -> None:
    """
    Performs robust center log-ratio transform robust PCA and
    ranks the features by the loadings of the resulting SVD.
    https://library.qiime2.org/plugins/deicode/19/
    Main per-dataset looper for the ADONIS tests on beta diversity matrices.

    :param i_data_sets_folder: Path to the folder containing the data/metadata subfolders.
    :param data_sets: list of data_sets.
    :param p_perm_groups: groups to subset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    """
    job_folder2 = get_job_folder(i_data_sets_folder, 'deicode/chunks')
    main_cases_dict = get_main_cases_dict(p_perm_groups)

    jobs = []
    all_sh_pbs = {}
    for dat, tsv_meta_pds in data_sets.items():

        tsv, meta = tsv_meta_pds
        meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')
        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'DEICODE')
        out_sh = '%s/run_adonis_%s.sh' % (job_folder2, dat)
        odir = get_analysis_folder(i_data_sets_folder, 'deicode/%s' % dat)
        for case_var, case_vals_list in cases_dict.items():
            cur_sh = '%s/run_beta_deicode_%s_%s.sh' % (job_folder2, dat, case_var)
            cur_sh = cur_sh.replace(' ', '-')
            all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
            p = multiprocessing.Process(
                target=run_multi_deicode,
                args=(odir, tsv, meta_pd, case_var, case_vals_list, out_sh, force))
            p.start()
            jobs.append(p)
    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_data_sets_folder, 'deicode')
    main_sh = write_main_sh(job_folder, '3_run_beta_deicode', all_sh_pbs,
                            '%s.dcd' % prjct_nm, '2', '1', '1', '200', 'mb', qiime_env)
    if main_sh:
        if p_perm_groups:
            print('# DEICODE (groups config in %s)' % p_perm_groups)
        else:
            print('# DEICODE')
        print('sh', main_sh)
