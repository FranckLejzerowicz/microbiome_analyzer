# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, glob
import pandas as pd
from os.path import basename, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    read_yaml_file,
    read_meta_pd,
    write_main_sh
)
from routine_qiime2_analyses._routine_q2_cmds import (
    filter_feature_table,
    run_add_metadata,
    write_nestedness,
    get_new_meta_pd,
    get_case, run_export
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict


def run_single_nestedness(odir: str, group: str, meta_pd: pd.DataFrame, nodfs: list,
                          nulls: list, modes: list, cur_sh: str, qza: str, case: str,
                          case_var: str, case_vals: list, binary: str, force: bool) -> dict:
    res = {}
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        if group:
            cur_rad = '%s/%s_%s_%s' % (odir, splitext(basename(qza))[0], group, case)
        else:
            cur_rad = '%s/%s_%s' % (odir, splitext(basename(qza))[0], case)
        new_meta = '%s.meta' % cur_rad
        new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
        cols = set()
        lat_lon_date = ['latitude', 'longitude', 'datetime']
        nodfs_valid = []
        for col in (nodfs + lat_lon_date):
            if col not in set(new_meta_pd.columns):
                continue
            if new_meta_pd[col].unique().size == 1:
                continue
            if col not in lat_lon_date and min(new_meta_pd[col].value_counts()) == 1:
                continue
            cols.add(col)
            if col in nodfs:
                nodfs_valid.append(col)
        new_meta_pd = new_meta_pd[sorted(cols)].reset_index()
        new_meta_pd.columns = (['#SampleID'] + sorted(cols))
        new_meta_pd.to_csv(new_meta, index=False, sep='\t')
        new_qza = '%s.qza' % cur_rad
        new_biom = '%s.biom' % cur_rad
        new_tsv = '%s.tsv' % cur_rad
        new_biom_meta = '%s_w-md.biom' % cur_rad

        if not isfile(new_biom):
            cmd = filter_feature_table(qza, new_qza, new_meta)
            cmd += run_export(new_qza, new_tsv, 'FeatureTable')
            cur_sh_o.write('echo "%s"\n' % cmd)
            cur_sh_o.write(cmd)

        cmd = run_add_metadata(new_biom, new_biom_meta, new_meta)
        cur_sh_o.write('echo "%s"\n' % cmd)
        cur_sh_o.write(cmd)

        for null in nulls:
            for mode in modes:
                odir = '%s/null-%s/mode-%s' % (cur_rad, null, mode)
                graphs = '%s/graphs.csv' % odir
                fields = '%s/fields.txt' % odir
                res[(null, mode)] = odir
                if not isdir(odir):
                    os.makedirs(odir)
                if not isfile(graphs) and not len(glob.glob('%s/*comparisons.csv' % odir)):
                    write_nestedness(new_biom_meta, odir, graphs, fields, binary,
                                     nodfs_valid, null, mode, cur_sh_o)
                    remove = False
    if remove:
        os.remove(cur_sh)
    return res


def get_nestedness_config(nestedness_config: dict) -> (dict, dict, dict, dict, dict):
    subsets = {'ALL': [[]]}
    if 'subsets' in nestedness_config:
        subsets.update(nestedness_config['subsets'])
    nodfs = []
    if 'nodfs' in nestedness_config:
        nodfs.extend(nestedness_config['nodfs'])
    colors = []
    if 'colors' in nestedness_config:
        colors.extend(nestedness_config['colors'])
    nulls = ['equiprobablefixed']
    if 'nulls' in nestedness_config:
        nulls = nestedness_config['nulls']
    modes = ['betweeneachpairoftypes']
    if 'modes' in nestedness_config:
        modes = nestedness_config['modes']
    return subsets, nodfs, colors, nulls, modes


def run_nestedness(i_datasets_folder: str, betas: dict, p_nestedness_groups: str,
                   datasets_rarefs: dict, force: bool, prjct_nm: str, qiime_env: str,
                   chmod: str, noloc: bool, split: bool, run_params: dict,
                   filt_raref: str, jobs: bool) -> dict:

    job_folder2 = get_job_folder(i_datasets_folder, 'nestedness/chunks')
    nestedness_config = read_yaml_file(p_nestedness_groups)
    if 'soft' not in nestedness_config:
        print('Must provide the path to the Nestedness soft (containing bin/Autocorrelation.jar)')
        return {}
    if nestedness_config['soft'].endswith('Autocorrelation.jar') and isfile(nestedness_config['soft']):
        binary = nestedness_config['soft']
    else:
        binary = '%s/bin/Autocorrelation.jar' % nestedness_config['soft']
        if not isfile(binary):
            print('Must provide the path to the Nestedness soft (containing bin/Autocorrelation.jar)')
            return {}
    subsets, nodfs, colors, nulls, modes = get_nestedness_config(nestedness_config)

    all_sh_pbs = {}
    nestedness_res = {}
    for dat, rarefs_metrics_groups_metas_qzas_dms_trees in betas.items():
        if not split:
            out_sh = '%s/run_nestedness_%s%s.sh' % (job_folder2, dat, filt_raref)
        nestedness_res[dat] = []
        for idx, metrics_groups_metas_qzas_dms_trees in enumerate(rarefs_metrics_groups_metas_qzas_dms_trees):
            nestedness_raref = {}
            cur_raref = datasets_rarefs[dat][idx]
            odir = get_analysis_folder(i_datasets_folder, 'nestedness/%s%s' % (dat, cur_raref))
            if split:
                out_sh = '%s/run_nestedness_%s%s%s.sh' % (job_folder2, dat, cur_raref, filt_raref)
            for _, groups_metas_qzas_dms_trees in metrics_groups_metas_qzas_dms_trees.items():
                for group, (meta, qza, __, ___) in groups_metas_qzas_dms_trees.items():
                    meta_pd = read_meta_pd(meta).set_index('sample_name')
                    cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(subsets), 'nestedness')
                    for case_var, case_vals_list in cases_dict.items():
                        for case_vals in case_vals_list:
                            case = get_case(case_vals, case_var).replace(' ', '_')
                            cur_sh = '%s/run_nestedness_%s%s_%s_%s%s.sh' % (
                                job_folder2, dat, cur_raref, group, case, filt_raref)
                            cur_sh = cur_sh.replace(' ', '-')
                            all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                            res = run_single_nestedness(odir, group, meta_pd, nodfs, nulls,
                                                        modes, cur_sh, qza, case, case_var,
                                                        case_vals, binary, force)
                            nestedness_raref[(group, case)] = res
                break
            nestedness_res[dat].append(nestedness_raref)

    job_folder = get_job_folder(i_datasets_folder, 'nestedness')
    main_sh = write_main_sh(job_folder, '3_run_nestedness%s' % filt_raref, all_sh_pbs,
                            '%s.prm%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs)
    if main_sh:
        if p_nestedness_groups:
            print("# nestedness (config in %s)" % p_nestedness_groups)
        else:
            print("# nestedness")
        print_message('', 'sh', main_sh, jobs)

    return nestedness_res


def nestedness_figure(nestedness_res, datasets_rarefs):
    for dat, nestedness_rarefs in nestedness_res.items():
        for idx, nestedness_raref in enumerate(nestedness_rarefs):
            cur_raref = datasets_rarefs[dat][idx]
            for (group, case), res in nestedness_raref.items():
                for (null, mode), odir in res.items():
                    print(dat)
                    print(cur_raref)
                    print(group, case)
                    print(null, mode)
                    print(odir)
                    print(glob.glob('%s/*o' % dir))
                    break
                break
            break
        break
