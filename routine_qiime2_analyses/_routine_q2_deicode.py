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
    get_job_folder,
    get_analysis_folder,
    get_main_cases_dict,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import (
    write_deicode_biplot,
    get_case, get_new_meta_pd
)


def run_single_deicode(odir: str, tsv: str, meta_pd: pd.DataFrame, case_var: str,
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
            cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s' % case)
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
    if remove:
        os.remove(cur_sh)


def run_deicode(i_datasets_folder: str, datasets: dict, datasets_rarefs: dict,
                p_perm_groups: str, force: bool, prjct_nm: str, qiime_env: str,
                chmod: str, noloc: bool, run_params: dict,
                filt_raref: str) -> None:
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
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder2 = get_job_folder(i_datasets_folder, 'deicode/chunks')
    main_cases_dict = get_main_cases_dict(p_perm_groups)
    # jobs = []
    all_sh_pbs = {}
    for dat, tsv_meta_pds_ in datasets.items():
        out_sh = '%s/run_deicode_%s%s.sh' % (job_folder2, dat, filt_raref)
        for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
            cur_raref = datasets_rarefs[dat][idx]
            tsv, meta = tsv_meta_pds
            meta_pd = read_meta_pd(meta)
            meta_pd = meta_pd.set_index('sample_name')
            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'DEICODE')
            odir = get_analysis_folder(i_datasets_folder, 'deicode/%s%s' % (dat, cur_raref))
            for case_var, case_vals_list in cases_dict.items():
                cur_sh = '%s/run_beta_deicode_%s%s_%s%s.sh' % (job_folder2, dat, cur_raref, case_var, filt_raref)
                cur_sh = cur_sh.replace(' ', '-')
                all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                run_single_deicode(odir, tsv, meta_pd, case_var,
                                   case_vals_list, cur_sh, force)

    job_folder = get_job_folder(i_datasets_folder, 'deicode')
    main_sh = write_main_sh(job_folder, '3_run_beta_deicode%s' % filt_raref, all_sh_pbs,
                            '%s.dcd%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_perm_groups:
            if p_perm_groups.startswith('/panfs'):
                p_perm_groups = p_perm_groups.replace(os.getcwd(), '')
            print('# DEICODE (groups config in %s)' % p_perm_groups)
        else:
            print('# DEICODE')
        print_message('', 'sh', main_sh)
