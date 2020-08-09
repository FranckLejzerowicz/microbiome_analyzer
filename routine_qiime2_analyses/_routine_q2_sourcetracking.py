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
from os.path import isdir, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_sourcetracking_config,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import (
    get_case, get_new_meta_pd, write_sourcetracking)


def run_single_sourcetracking(
        odir: str, tsv: str, meta_pd: pd.DataFrame, case_var: str,
        sourcetracking_params: dict, method: str, sourcetracking_sourcesink: dict,
        case_vals_list: list, cur_sh: str, cur_import_sh: str,
        force: bool, filt: str, cur_raref: str, fp: str,  fa: str,
        n_nodes: str, n_procs: str) -> list:

    cases = []
    remove = True
    qza = '%s.qza' % splitext(tsv)[0]
    with open(cur_sh, 'w') as cur_sh_o, open(cur_import_sh, 'w') as cur_import_sh_o:
        for case_vals in case_vals_list:
            case = get_case(case_vals, '', case_var)
            for sourcesink_name, sourcesink_d in sourcetracking_sourcesink.items():
                column = sourcesink_d['column']
                sink = sourcesink_d['sink']
                sources = ['']
                if 'source' in sourcesink_d:
                    sources = sourcesink_d['source']
                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                if column not in new_meta_pd.columns:
                    raise IOError('"%s" not in metadata...' % column)
                if sink not in set(new_meta_pd[column].unique()):
                    raise IOError('All sinks "%s" not in metadata column "%s"' % (sink, column))
                if sources != [''] and not set(sources).issubset(set(new_meta_pd[column].unique())):
                    raise IOError('All sources "%s" not in metadata column "%s"' % (sources, column))

                cur_rad = '%s/%s_%s%s/%s' % (odir, case.strip('_'), filt, cur_raref, sourcesink_name)
                if not isdir(cur_rad):
                    os.makedirs(cur_rad)

                folder = '%s/%s-%s_%s' % (
                    cur_rad,
                    column,
                    sink.replace(
                        '/', '').replace(
                        '(', '').replace(
                        ')', '').replace(
                        ' ', ''),
                    '_'.join([source.replace(
                        '/', '').replace(
                        '(', '').replace(
                        ')', '').replace(
                        ' ', '') for source in sources])
                )

                new_meta = '%s/meta.tsv' % cur_rad
                new_qza = '%s/tab.qza' % cur_rad
                new_tsv = '%s/tab.tsv' % cur_rad
                new_meta_pd[[column]].reset_index().to_csv(new_meta, index=False, sep='\t')

                folder_method = folder + '/' + method
                if method == 'q2':
                    outs = folder_method + '/t0/r0/*/predictions.tsv'
                elif method == 'feast':
                    outs = folder_method + '/t0/out.r0*'
                elif method == 'sourcetracker':
                    outs = folder_method + '/t0/r0/mixing_proportions.txt'

                if force or not len(glob.glob(outs)):
                    write_sourcetracking(
                        qza, new_qza, new_tsv, new_meta, method, fp, fa,
                        cur_rad, column, sink, sources, sourcetracking_params,
                        n_nodes, n_procs, cur_sh_o, cur_import_sh_o)
                    cur_sh_o.write('echo "sh %s/cmd_%s.sh"\n' % (folder_method, method))
                    cur_sh_o.write('sh %s/cmd_%s.sh\n\n\n' % (folder_method, method))
                    remove = False
    if remove:
        os.remove(cur_sh)
    return cases


def run_sourcetracking(i_datasets_folder: str, datasets: dict, p_sourcetracking_config: str,
                       datasets_rarefs: dict, force: bool, prjct_nm: str, qiime_env: str,
                       chmod: str, noloc: bool, run_params: dict, filt_raref: str, split: bool) -> None:

    job_folder2 = get_job_folder(i_datasets_folder, 'sourcetracking/chunks')
    sourcetracking_dicts = get_sourcetracking_config(p_sourcetracking_config)
    sourcetracking_sourcesink = sourcetracking_dicts[0]
    sourcetracking_filtering = sourcetracking_dicts[1]
    sourcetracking_params = sourcetracking_dicts[2]
    main_cases_dict = sourcetracking_dicts[3]

    all_sh_pbs = {}
    all_import_sh_pbs = {}
    for dat, tsv_meta_pds_ in datasets.items():
        if dat in sourcetracking_filtering:
            filters = sourcetracking_filtering[dat]
        else:
            filters = {'0_0': ['0', '0']}
        for idx, tsv_meta_pds in enumerate(tsv_meta_pds_):
            tsv, meta = tsv_meta_pds
            meta_pd = read_meta_pd(meta)
            meta_pd = meta_pd.set_index('sample_name')
            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'sourcetracking')
            cur_raref = datasets_rarefs[dat][idx]
            if not split:
                out_sh = '%s/run_sourcetracking_%s%s%s.sh' % (job_folder2, dat, filt_raref, cur_raref)
                out_import_sh = '%s/run_import_sourcetracking_%s%s%s.sh' % (job_folder2, dat, filt_raref, cur_raref)
            odir = get_analysis_folder(i_datasets_folder, 'sourcetracking/%s' % dat)
            for method in sourcetracking_params['method']:
                if split:
                    out_sh = '%s/run_sourcetracking_%s%s%s_%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, method)
                    out_import_sh = '%s/run_import_sourcetracking_%s%s%s_%s.sh' % (
                        job_folder2, dat, filt_raref, cur_raref, method)
                for case_var, case_vals_list in cases_dict.items():
                    for filt, (fp, fa) in filters.items():
                        cur_sh = '%s/run_sourcetracking_%s_%s%s%s_%s.sh' % (
                            job_folder2, dat, case_var, filt_raref, cur_raref, filt)
                        cur_sh = cur_sh.replace(' ', '-')
                        cur_import_sh = '%s/run_import_sourcetracking_%s_%s%s%s_%s.sh' % (
                            job_folder2, dat, case_var, filt_raref, cur_raref, filt)
                        cur_import_sh = cur_import_sh.replace(' ', '-')
                        all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                        all_import_sh_pbs.setdefault((dat, out_import_sh), []).append(cur_import_sh)
                        run_single_sourcetracking(
                            odir, tsv, meta_pd, case_var, sourcetracking_params, method,
                            sourcetracking_sourcesink, case_vals_list, cur_sh, cur_import_sh, force,
                            filt, cur_raref, fp, fa, run_params["n_nodes"], run_params["n_procs"])

    job_folder = get_job_folder(i_datasets_folder, 'sourcetracking')
    main_sh = write_main_sh(job_folder, '3_run_import_sourcetracking%s' % filt_raref,
                            all_import_sh_pbs, '%s.mpt.srctrk%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, '~/.')
    if main_sh:
        if p_sourcetracking_config:
            if p_sourcetracking_config.startswith('/panfs'):
                p_sourcetracking_config = p_sourcetracking_config.replace(os.getcwd(), '')
            print('# import sourcetracking (groups config in %s)' % p_sourcetracking_config)
        else:
            print('# import sourcetracking')
        print_message('', 'sh', main_sh)

    main_sh = write_main_sh(job_folder, '3_run_sourcetracking%s' % filt_raref, all_sh_pbs,
                            '%s.srctrk%s' % (prjct_nm, filt_raref), run_params["time"],
                            run_params["n_nodes"], run_params["n_procs"], run_params["mem_num"],
                            run_params["mem_dim"], 'sourcetracker', chmod, noloc, '~/.')
    if main_sh:
        if p_sourcetracking_config:
            if p_sourcetracking_config.startswith('/panfs'):
                p_sourcetracking_config = p_sourcetracking_config.replace(os.getcwd(), '')
            print('# sourcetracking (groups config in %s)' % p_sourcetracking_config)
        else:
            print('# sourcetracking')
        print_message('', 'sh', main_sh)
