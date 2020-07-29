# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import time
import random
import numpy as np
import pandas as pd
from os.path import isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_doc_config,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import get_case, get_new_meta_pd, write_doc


def run_single_doc(i_dataset_folder: str, odir: str, tsv: str,
                   meta_pd: pd.DataFrame, case_var: str,
                   doc_params: dict, case_vals_list: list, cur_sh: str,
                   cur_import_sh: str, force: bool, filt: str, cur_raref: str,
                   fp: str, fa: str, n_nodes: str, n_procs: str) -> list:
    remove = True
    qza = '%s.qza' % splitext(tsv)[0]
    cases = []
    with open(cur_sh, 'w') as cur_sh_o, open(cur_import_sh, 'w') as cur_import_sh_o:
        for case_vals in case_vals_list:
            token = ''.join([str(random.choice(range(100))) for x in range(3)])
            case = get_case(case_vals, '', case_var)
            cur_rad = '%s/%s_%s%s' % (odir, case.strip('_'), filt, cur_raref)
            cases.append(cur_rad)
            cur_rad_r = '%s/R' % cur_rad
            cur_rad_token = '%s/tmp/%s' % (i_dataset_folder, token)
            if not isdir(cur_rad_r):
                os.makedirs(cur_rad_r)
            if not isdir(cur_rad_token):
                os.makedirs(cur_rad_token)
            new_meta = '%s/meta.tsv' % cur_rad
            new_qza = '%s/tab.qza' % cur_rad
            new_tsv = '%s/tab.tsv' % cur_rad
            new_tsv_token = '%s/tab.tsv' % cur_rad_token
            if force or not isfile('%s/DO.tsv' % cur_rad):
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                write_doc(qza, fp, fa, new_meta, new_qza, new_tsv,
                          cur_rad, new_tsv_token, cur_rad_token,
                          n_nodes, n_procs, doc_params,
                          cur_sh_o, cur_import_sh_o)
                remove = False
    if remove:
        os.remove(cur_sh)
    return cases


def run_doc(i_datasets_folder: str, datasets: dict, p_doc_config: str,
            datasets_rarefs: dict, force: bool, prjct_nm: str,
            qiime_env: str, chmod: str, noloc: bool, run_params: dict,
            filt_raref: str, eval_depths: dict, split: bool) -> None:

    evaluation = ''
    if len(eval_depths):
        evaluation = '_eval'

    job_folder2 = get_job_folder(i_datasets_folder, 'doc%s/chunks' % evaluation)
    doc_config, doc_params, main_cases_dict = get_doc_config(p_doc_config)

    all_sh_pbs = {}
    all_import_sh_pbs = {}
    dat_cases_tabs = {}
    for dat, tsv_meta_pds_ in datasets.items():
        dat_cases_tabs[dat] = {}
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
            dat_cases_tabs[dat][cur_raref] = {}
            if not split:
                out_sh = '%s/run_doc%s_%s%s%s.sh' % (job_folder2, evaluation, dat, filt_raref, cur_raref)
                out_import_sh = '%s/run_import_doc%s_%s%s%s.sh' % (job_folder2, evaluation, dat, filt_raref, cur_raref)
            odir = get_analysis_folder(i_datasets_folder, 'doc%s/%s' % (evaluation, dat))
            for case_var, case_vals_list in cases_dict.items():
                if split:
                    out_sh = '%s/run_doc%s_%s%s%s_%s.sh' % (
                        job_folder2, evaluation, dat, filt_raref, cur_raref, case_var)
                    out_import_sh = '%s/run_import_doc%s_%s%s%s_%s.sh' % (
                        job_folder2, evaluation, dat, filt_raref, cur_raref, case_var)
                for filt, (fp, fa) in filters.items():
                    cur_sh = '%s/run_doc%s_%s_%s%s%s_%s.sh' % (
                        job_folder2, evaluation, dat, case_var, filt_raref, cur_raref, filt)
                    cur_sh = cur_sh.replace(' ', '-')
                    cur_import_sh = '%s/run_import_doc%s_%s_%s%s%s_%s.sh' % (
                        job_folder2, evaluation, dat, case_var, filt_raref, cur_raref, filt)
                    cur_import_sh = cur_import_sh.replace(' ', '-')
                    all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                    all_import_sh_pbs.setdefault((dat, out_import_sh), []).append(cur_import_sh)
                    cases = run_single_doc(i_datasets_folder, odir, tsv, meta_pd, case_var, doc_params,
                                           case_vals_list, cur_sh, cur_import_sh, force, filt, cur_raref,
                                           fp, fa, run_params["n_nodes"], run_params["n_procs"])
                    dat_cases_tabs[dat][cur_raref].setdefault(case_var, []).extend(cases)
    job_folder = get_job_folder(i_datasets_folder, 'doc%s' % evaluation)

    main_sh = write_main_sh(job_folder, '3_run_import_doc%s%s' % (evaluation, filt_raref),
                            all_import_sh_pbs, '%s.doc.mpt%s%s' % (prjct_nm, evaluation, filt_raref),
                            "1", "1", "1", "250", "mb", qiime_env, chmod, noloc)
    if main_sh:
        if p_doc_config:
            if p_doc_config.startswith('/panfs'):
                p_doc_config = p_doc_config.replace(os.getcwd(), '')
            print('# Import for DOC (groups config in %s)' % p_doc_config)
        else:
            print('# Import DOC')
        print_message('', 'sh', main_sh)

    main_sh = write_main_sh(job_folder, '3_run_doc%s%s' % (evaluation, filt_raref), all_sh_pbs,
                            '%s.doc%s%s' % (prjct_nm, evaluation, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, '~/.')
    if main_sh:
        if p_doc_config:
            if p_doc_config.startswith('/panfs'):
                p_doc_config = p_doc_config.replace(os.getcwd(), '')
            print('# DOC (groups config in %s)' % p_doc_config)
        else:
            print('# DOC')
        print_message('', 'sh', main_sh)

    do_r = 1
    if do_r and not len(eval_depths):
        job_folder = get_job_folder(i_datasets_folder, 'doc/R')
        job_folder2 = get_job_folder(i_datasets_folder, 'doc/R/chunks')
        main_written = 0
        main_sh = '%s/run_R_doc%s.sh' % (job_folder, filt_raref)
        with open(main_sh, 'w') as main_o:
            for dat, raref_case_var_cases in dat_cases_tabs.items():
                shs = []
                written = 0
                for raref, case_var_cases in raref_case_var_cases.items():
                    for case_var, cases in case_var_cases.items():
                        for cdx, case in enumerate(cases):
                            plot = '%s_%s_%s_%s' % (dat, raref, case_var, cdx)
                            if not isfile('%s/R/DO.tsv' % case):
                                cur_r = '%s/run_R_doc_%s_%s_%s_vanilla.R' % (job_folder2, dat, case_var, cdx)
                                shs.append('R -f %s --vanilla\n' % cur_r)
                                with open(cur_r, 'w') as o:
                                    o.write("library(DOC)\n")
                                    o.write("library(ggplot2)\n")
                                    o.write("otu <- read.table('%s/tab.tsv', header=T, sep='\\t', comment.char='', check.names=F, nrows=2)\n" % case)
                                    o.write("index_name <- colnames(otu)[1]\n")
                                    o.write("otu <- read.table('%s/tab.tsv', header=T, sep='\\t', comment.char='', check.names=F, row.names=index_name)\n" % case)
                                    o.write("if (dim(otu)[1] > 100) {\n")
                                    o.write("    res <- DOC(otu)\n")
                                    o.write("    res.null <- DOC.null(otu)\n")
                                    o.write("    write.table(x=res$DO, file='%s/R/DO.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res$LME, file='%s/R/LME.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    colnames(res$NEG) <- c('Neg_Slope', 'Data')\n")
                                    o.write("    write.table(x=res$NEG, file='%s/R/NEG.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res$FNS, file='%s/R/FNS.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res$BOOT, file='%s/R/BOOT.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res$CI, file='%s/R/CI.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res.null$DO, file='%s/R/null_DO.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res.null$LME, file='%s/R/null_LME.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    colnames(res.null$NEG) <- c('Neg_Slope', 'Data')\n")
                                    o.write("    write.table(x=res.null$NEG, file='%s/R/null_NEG.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res.null$FNS, file='%s/R/null_FNS.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res.null$BOOT, file='%s/R/null_BOOT.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    write.table(x=res.null$CI, file='%s/R/null_CI.tsv', sep='\\t', quote=F, row.names=F)\n" % case)
                                    o.write("    pdf('%s/R/plot.pdf')\n" % case)
                                    o.write("    merged <- DOC.merge(list(s_%s = res, s_%s = res.null))\n" % (plot, plot))
                                    o.write("    plot(merged)\n")
                                    o.write("    dev.off()\n")
                                    o.write("}\n")
                                    main_written += 1
                                    written += 1
                if written:
                    if split and len(shs) >= 3:
                        chunks = [list(x) for x in np.array_split(np.array(shs), 3)]
                    else:
                        chunks = [shs]
                    for cdx, chunk in enumerate(chunks):
                        out_sh = '%s/run_R_doc_%s%s_%s.sh' % (job_folder2, dat, filt_raref, cdx)
                        out_pbs = '%s.pbs' % splitext(out_sh)[0]
                        with open(out_sh, 'w') as o:
                            for c in chunk:
                                o.write('echo "%s"\n\n' % c)
                                o.write('%s\n\n' % c)
                        run_xpbs(out_sh, out_pbs, '%s.doc.R.%s%s_%s' % (prjct_nm, dat, filt_raref, cdx),
                                 'xdoc', run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                                 run_params["mem_num"], run_params["mem_dim"],
                                 chmod, written, 'single', main_o, noloc)
        if main_written:
            print_message('# DOC (R)', 'sh', main_sh)

