# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import pandas as pd
from os.path import isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_phate_dicts,
    write_main_sh,
    read_meta_pd,
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case, write_phate_cmd
)


def run_single_phate(dat: str, odir: str, tsv: str, meta_pd: pd.DataFrame, case_var: str,
                     phate_labels: list, phate_params: dict, run_params: dict,
                     case_vals_list: list, cur_sh: str, cur_import_sh: str, force: bool,
                     filt: str, cur_raref: str, fp: str, fa: str) -> dict:

    remove = True
    qza = '%s.qza' % splitext(tsv)[0]
    cases = {}
    with open(cur_sh, 'w') as cur_sh_o, open(cur_import_sh, 'w') as cur_import_sh_o:
        for case_vals in case_vals_list:
            case = get_case(case_vals, '', case_var)
            cur_rad = '%s/%s_%s%s' % (odir, case.strip('_'), filt, cur_raref)
            if not isdir(cur_rad):
                os.makedirs(cur_rad)
            new_meta = '%s/meta.tsv' % cur_rad
            new_qza = '%s/tab.qza' % cur_rad
            new_tsv = '%s/tab.tsv' % cur_rad
            phate_html = '%s/phate_%s_%s_%s.html' % (cur_rad, dat, filt, case)
            phate_tsv = '%s_xphate.tsv' % splitext(phate_html)[0]
            if len(glob.glob('%s/TOO_FEW.*' % cur_rad)):
                continue
            cases[case] = phate_tsv
            if force or not isfile(phate_html) or not isfile(phate_tsv):
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                write_phate_cmd(qza, new_qza, new_tsv, new_meta, fp, fa,
                                phate_html, phate_labels, phate_params,
                                run_params["n_nodes"], run_params["n_procs"],
                                cur_sh_o, cur_import_sh_o)
                remove = False
    if remove:
        os.remove(cur_sh)
    return cases


def run_phate(p_phate_config: str, i_datasets_folder: str, datasets: dict,
              datasets_rarefs: dict, force: bool, prjct_nm: str,
              qiime_env: str, chmod: str, noloc: bool, split: bool,
              run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> dict:

    job_folder2 = get_job_folder(i_datasets_folder, 'phate/chunks')
    phate_dicts = get_phate_dicts(p_phate_config)
    phate_filtering, phate_labels, phate_params, main_cases_dict = phate_dicts

    phates = {}
    all_sh_pbs = {}
    all_import_sh_pbs = {}
    for dat, tsv_meta_pds_ in datasets.items():
        phates[dat] = []
        if dat in phate_filtering:
            filters = phate_filtering[dat]
        else:
            filters = {'0_0': ['0', '0']}
        for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
            tsv, meta = tsv_meta_pds
            meta_pd = read_meta_pd(meta)
            meta_pd = meta_pd.set_index('sample_name')
            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'phate')
            cur_raref = datasets_rarefs[dat][idx]
            if not split:
                out_sh = '%s/run_phate_%s%s%s.sh' % (
                    job_folder2, dat, filt_raref, cur_raref)
                out_import_sh = '%s/run_import_phate_%s%s%s.sh' % (
                    job_folder2, dat, filt_raref, cur_raref)
            odir = get_analysis_folder(i_datasets_folder, 'phate/%s' % dat)
            raref_phates = {}
            for filt, (fp, fa) in filters.items():
                raref_phates[filt] = {}
                if split:
                    out_sh = '%s/run_phate_%s%s%s%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, filt)
                    out_import_sh = '%s/run_import_phate_%s%s%s%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, filt)
                for case_var, case_vals_list in cases_dict.items():
                    cur_sh = '%s/run_phate_%s%s%s_%s_%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, case_var, filt)
                    cur_sh = cur_sh.replace(' ', '-')
                    cur_import_sh = '%s/run_import_phate_%s%s%s_%s_%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, case_var, filt)
                    cur_import_sh = cur_import_sh.replace(' ', '-')
                    all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                    all_import_sh_pbs.setdefault((dat, out_import_sh), []).append(cur_import_sh)
                    phate = run_single_phate(
                        dat, odir, tsv, meta_pd, case_var, phate_labels,
                        phate_params, run_params, case_vals_list,
                        cur_sh, cur_import_sh, force, filt, cur_raref, fp, fa)
                    raref_phates[filt][case_var] = phate
            phates[dat].append(raref_phates)

    job_folder = get_job_folder(i_datasets_folder, 'phate')
    main_sh = write_main_sh(job_folder, '3_run_import_phate%s' % filt_raref, all_import_sh_pbs,
                            '%s.mrt.pht%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        if p_phate_config:
            if p_phate_config.startswith('/panfs'):
                p_phate_config = p_phate_config.replace(os.getcwd(), '')
            print('# Import for PHATE (groups config in %s)' % p_phate_config)
        else:
            print('# Import for PHATE')
        print_message('', 'sh', main_sh, jobs)

    main_sh = write_main_sh(job_folder, '3_run_phate%s' % filt_raref, all_sh_pbs,
                            '%s.pht%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            'xphate', chmod, noloc, jobs)
    if main_sh:
        if p_phate_config:
            if p_phate_config.startswith('/panfs'):
                p_phate_config = p_phate_config.replace(os.getcwd(), '')
            print('# PHATE (groups config in %s)' % p_phate_config)
        else:
            print('# PHATE')
        print_message('', 'sh', main_sh, jobs)
    return phates
