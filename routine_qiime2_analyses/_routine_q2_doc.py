# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
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
                   fp: str, fa: str, n_nodes: str, n_procs: str,
                   dat_phates: dict, doc_phate: bool, need_to_run_phate: list,
                   need_to_run_less_phate: list) -> list:
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

            # run DOC on each cluster from PHATE
            if doc_phate and filt in dat_phates and case in dat_phates[filt]:
                # get the clusters
                xphate_tsv = dat_phates[filt][case]
                print('xphate_tsv:', xphate_tsv)
                if not isfile(xphate_tsv):
                    if not need_to_run_phate:
                        print('Unable to run DOC on a set of PHATE clusters:\n'
                              '(Be sure to run PHATE first...)')
                    need_to_run_phate.append(
                        xphate_tsv.replace('%s/qiime/phate' % i_dataset_folder, '...'))
                    continue
                xphate_pd = pd.read_csv(xphate_tsv, header=0, sep='\t', dtype={'sample_name': str})
                xphate_pd = xphate_pd.loc[xphate_pd['variable'].str.contains('cluster_k')]
                if len(xphate_pd[['knn', 'decay', 't']].drop_duplicates()) > 5:
                    if not need_to_run_less_phate:
                        print('Warning: PHATE has been for multiple parameters combinations:\n'
                              ' --> It may be unwise to let DOC run on every combination...\n'
                              ' --> Be sure to run PHATE using few, desired sets of parameters!')
                    need_to_run_less_phate.append(
                        xphate_tsv.replace('%s/qiime/phate' % i_dataset_folder, '...'))
                cols = ['sample_name', 'knn', 'decay', 't', 'variable', 'factor']
                xphate_clusters = dict(xphate_pd[cols].groupby(
                    ['knn', 'decay', 't', 'variable', 'factor']
                ).apply(func=lambda x: x.sample_name.tolist()))
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                # repeat DOC command for the clusters
                for (knn, decay, t, k, cluster), samples_phate in xphate_clusters.items():
                    token = ''.join([str(random.choice(range(100))) for x in range(3)])
                    cur_rad_phate = '%s/phate/%s_%s_%s_k%s_clust%s' % (cur_rad, knn, decay, t, k, cluster)
                    cases.append(cur_rad_phate)
                    cur_rad_phate_r = '%s/R' % cur_rad_phate
                    cur_rad_token = '%s/tmp/%s' % (i_dataset_folder, token)
                    if not isdir(cur_rad_phate_r):
                        os.makedirs(cur_rad_phate_r)
                    new_meta = '%s/meta.tsv' % cur_rad_phate
                    new_qza = '%s/tab.qza' % cur_rad_phate
                    new_tsv = '%s/tab.tsv' % cur_rad_phate
                    new_tsv_token = '%s/tab.tsv' % cur_rad_token
                    if force or not isfile('%s/DO.tsv' % cur_rad_phate):
                        new_meta_pd_phate = new_meta_pd.loc[samples_phate, :].copy()
                        new_meta_pd_phate.reset_index().to_csv(new_meta, index=False, sep='\t')
                        write_doc(qza, fp, fa, new_meta, new_qza, new_tsv,
                                  cur_rad_phate, new_tsv_token, cur_rad_token,
                                  n_nodes, n_procs, doc_params,
                                  cur_sh_o, cur_import_sh_o)
                        remove = False
    if remove:
        os.remove(cur_sh)
    return cases


def run_doc(i_datasets_folder: str, datasets: dict, p_doc_config: str,
            datasets_rarefs: dict, force: bool, prjct_nm: str,
            qiime_env: str, chmod: str, noloc: bool, run_params: dict,
            filt_raref: str, phates: dict, doc_phate: bool, split: bool) -> None:

    job_folder2 = get_job_folder(i_datasets_folder, 'doc/chunks')
    doc_filtering, doc_params, main_cases_dict = get_doc_config(p_doc_config)

    all_sh_pbs = {}
    all_import_sh_pbs = {}
    dat_cases_tabs = {}
    need_to_run_phate = []
    need_to_run_less_phate = []
    for dat, tsv_meta_pds_ in datasets.items():
        dat_cases_tabs[dat] = {}
        if dat in doc_filtering:
            filters = doc_filtering[dat]
        else:
            filters = {'0_0': ['0', '0']}
        for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
            dat_phates = []
            if dat in phates:
                dat_phates = phates[dat][idx]
            tsv, meta = tsv_meta_pds
            meta_pd = read_meta_pd(meta)
            meta_pd = meta_pd.set_index('sample_name')
            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'DOC')
            cur_raref = datasets_rarefs[dat][idx]
            dat_cases_tabs[dat][cur_raref] = {}
            if not split:
                out_sh = '%s/run_doc_%s%s%s.sh' % (job_folder2, dat, filt_raref, cur_raref)
                out_import_sh = '%s/run_import_doc_%s%s%s.sh' % (job_folder2, dat, filt_raref, cur_raref)
            odir = get_analysis_folder(i_datasets_folder, 'doc/%s' % dat)
            for filt, (fp, fa) in filters.items():
                if split:
                    out_sh = '%s/run_doc_%s%s%s_%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, filt)
                    out_import_sh = '%s/run_import_doc_%s%s%s_%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, filt)
                for case_var, case_vals_list in cases_dict.items():
                    cur_sh = '%s/run_doc_%s_%s%s%s_%s.sh' % (
                        job_folder2, dat, case_var, filt_raref, cur_raref, filt)
                    cur_sh = cur_sh.replace(' ', '-')
                    cur_import_sh = '%s/run_import_doc_%s_%s%s%s_%s.sh' % (
                        job_folder2, dat, case_var, filt_raref, cur_raref, filt)
                    cur_import_sh = cur_import_sh.replace(' ', '-')
                    all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                    all_import_sh_pbs.setdefault((dat, out_import_sh), []).append(cur_import_sh)
                    cases = run_single_doc(i_datasets_folder, odir, tsv, meta_pd, case_var, doc_params,
                                           case_vals_list, cur_sh, cur_import_sh, force, filt, cur_raref,
                                           fp, fa, run_params["n_nodes"], run_params["n_procs"],
                                           dat_phates, doc_phate, need_to_run_phate, need_to_run_less_phate)
                    dat_cases_tabs[dat][cur_raref].setdefault(case_var, []).extend(cases)

    for need_to_run in need_to_run_phate:
        print(' -', need_to_run)

    job_folder = get_job_folder(i_datasets_folder, 'doc')
    main_sh = write_main_sh(job_folder, '3_run_import_doc%s' % filt_raref,
                            all_import_sh_pbs, '%s.doc.mpt%s' % (prjct_nm, filt_raref),
                            "4", "1", "1", "500", "mb", qiime_env, chmod, noloc)
    if main_sh:
        if p_doc_config:
            if p_doc_config.startswith('/panfs'):
                p_doc_config = p_doc_config.replace(os.getcwd(), '')
            print('# Import for DOC (groups config in %s)' % p_doc_config)
        else:
            print('# Import DOC')
        print_message('', 'sh', main_sh)

    main_sh = write_main_sh(job_folder, '3_run_doc%s' % filt_raref, all_sh_pbs,
                            '%s.doc%s' % (prjct_nm, filt_raref),
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
    if do_r:
        job_folder = get_job_folder(i_datasets_folder, 'doc/R')
        job_folder2 = get_job_folder(i_datasets_folder, 'doc/R/chunks')
        main_written = 0
        main_sh = '%s/run_R_doc%s.sh' % (job_folder, filt_raref)
        with open(main_sh, 'w') as main_o:
            for dat, raref_case_var_cases in dat_cases_tabs.items():

                shs = []
                written = 0
                odir = get_analysis_folder(i_datasets_folder, 'doc/%s' % dat)
                log_error = '%s/log.error' % odir
                for raref, case_var_cases in raref_case_var_cases.items():
                    for case_var, cases in case_var_cases.items():
                        for cdx, case in enumerate(cases):
                            plot = '%s_%s_%s_%s' % (dat, raref, case_var, cdx)
                            case_r = '%s/R' % case
                            pdf = '%s/plot.pdf' % case_r
                            if not isfile(pdf):
                                cur_r = '%s/run_R_doc_%s_%s_%s_vanilla.R' % (job_folder2, dat, case_var, cdx)
                                cur_sh = 'echo "*** %s" >> %s\n' % (plot, log_error)
                                cur_sh += 'R -f %s --vanilla 2>> %s\n' % (cur_r, log_error)
                                cur_sh += 'echo "end" >> %s\n' % log_error
                                shs.append(cur_sh)
                                with open(cur_r, 'w') as o:
                                    o.write("library(DOC)\n")
                                    o.write("library(ggplot2)\n")
                                    o.write("otu <- read.table('%s/tab.tsv', header=T, sep='\\t', comment.char='', check.names=F, nrows=2)\n" % case)
                                    o.write("index_name <- colnames(otu)[1]\n")
                                    o.write("otu <- read.table('%s/tab.tsv', header=T, sep='\\t', comment.char='', check.names=F, row.names=index_name)\n" % case)
                                    o.write("if (dim(otu)[1] > 100) {\n")
                                    o.write("    res <- DOC(otu)\n")
                                    o.write("    res.null <- DOC.null(otu)\n")
                                    o.write("    write.table(x=res$DO, file='%s/DO.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res$LME, file='%s/LME.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    colnames(res$NEG) <- c('Neg_Slope', 'Data')\n")
                                    o.write("    write.table(x=res$NEG, file='%s/NEG.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res$FNS, file='%s/FNS.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res$BOOT, file='%s/BOOT.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res$CI, file='%s/CI.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res.null$DO, file='%s/null_DO.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res.null$LME, file='%s/null_LME.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    colnames(res.null$NEG) <- c('Neg_Slope', 'Data')\n")
                                    o.write("    write.table(x=res.null$NEG, file='%s/null_NEG.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res.null$FNS, file='%s/null_FNS.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res.null$BOOT, file='%s/null_BOOT.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    write.table(x=res.null$CI, file='%s/null_CI.tsv', sep='\\t', quote=F, row.names=F)\n" % case_r)
                                    o.write("    colnames(res$NEG) <- c('Neg.Slope', 'Data')\n")
                                    o.write("    colnames(res.null$NEG) <- c('Neg.Slope', 'Data')\n")
                                    o.write("    pdf('%s')\n" % pdf)
                                    o.write("    merged <- DOC.merge(list(s_%s = res, s_%s=res.null))\n" % (plot, plot))
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

