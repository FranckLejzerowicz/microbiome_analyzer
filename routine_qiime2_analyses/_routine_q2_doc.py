# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
from os.path import isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_doc_config,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import get_case, get_new_meta_pd, write_doc


def run_single_doc(odir: str, tsv: str, meta_pd: pd.DataFrame, case_var: str,
                   doc_params: dict, case_vals_list: list, cur_sh: str,
                   force: bool, filt: str, cur_raref: str, fp: str, fa: str,
                   n_nodes: str, n_procs: str) -> None:
    remove = True
    qza = '%s.qza' % splitext(tsv)[0]
    with open(cur_sh, 'w') as cur_sh_o:
        for case_vals in case_vals_list:
            case = get_case(case_vals, '', case_var)
            cur_rad = '%s/%s_%s%s' % (odir, case.strip('_'), filt, cur_raref)
            if not isdir(cur_rad):
                os.makedirs(cur_rad)
            new_meta = '%s/meta.tsv' % cur_rad
            new_qza = '%s/tab.qza' % cur_rad
            new_tsv = '%s/tab.tsv' % cur_rad
            if force or not isfile('%s/DO.tsv' % cur_rad):
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                write_doc(qza, new_meta, new_qza, new_tsv, fp, fa,
                          cur_rad, n_nodes, n_procs, doc_params, cur_sh_o)
                remove = False
    if remove:
        os.remove(cur_sh)


def run_doc(i_datasets_folder: str, datasets: dict, p_doc_config: str,
            datasets_rarefs: dict, force: bool, prjct_nm: str,
            qiime_env: str, chmod: str, noloc: bool, run_params: dict,
            filt_raref: str, eval_depths: dict) -> None:

    evaluation = ''
    if len(eval_depths):
        evaluation = '_eval'

    job_folder2 = get_job_folder(i_datasets_folder, 'doc%s/chunks' % evaluation)
    doc_config, doc_params, main_cases_dict = get_doc_config(p_doc_config)

    all_sh_pbs = {}
    for dat, tsv_meta_pds_ in datasets.items():
        if doc_config and 'filtering' in doc_config and dat in doc_config['filtering']:
            filters = doc_config['filtering'][dat]
        else:
            filters = {'0_0': ['0', '0']}
        for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
            tsv, meta = tsv_meta_pds
            meta_pd = read_meta_pd(meta)
            meta_pd = meta_pd.set_index('sample_name')
            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'DOC')
            cur_raref = datasets_rarefs[dat][idx]
            out_sh = '%s/run_doc%s_%s%s%s.sh' % (job_folder2, evaluation, dat, filt_raref, cur_raref)
            odir = get_analysis_folder(i_datasets_folder, 'doc%s/%s' % (evaluation, dat))
            for case_var, case_vals_list in cases_dict.items():
                for filt, (fp, fa) in filters.items():
                    cur_sh = '%s/run_doc%s_%s_%s%s%s_%s.sh' % (
                        job_folder2, evaluation, dat, case_var, filt_raref, cur_raref, filt)
                    cur_sh = cur_sh.replace(' ', '-')
                    all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                    run_single_doc(odir, tsv, meta_pd, case_var, doc_params,
                                   case_vals_list, cur_sh, force, filt, cur_raref, fp, fa,
                                   run_params["n_nodes"], run_params["n_procs"])
    job_folder = get_job_folder(i_datasets_folder, 'doc%s' % evaluation)
    main_sh = write_main_sh(job_folder, '3_run_doc%s%s' % (evaluation, filt_raref), all_sh_pbs,
                            '%s.doc%s%s' % (prjct_nm, evaluation, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_doc_config:
            if p_doc_config.startswith('/panfs'):
                p_doc_config = p_doc_config.replace(os.getcwd(), '')
            print('# DOC (groups config in %s)' % p_doc_config)
        else:
            print('# DOC')
        print_message('', 'sh', main_sh)
