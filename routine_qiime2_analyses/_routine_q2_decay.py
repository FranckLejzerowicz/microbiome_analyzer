# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
from os.path import basename, isfile, splitext, isdir
import matplotlib.pyplot as plt
import seaborn as sns

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    read_yaml_file,
    read_meta_pd,
    write_main_sh
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_distance_decay,
    get_case,
    get_new_meta_pd
)


def get_decay_config(decay_config: dict) -> (dict, dict):

    subsets = {'ALL': [[]]}
    if 'subsets' in decay_config:
        subsets.update(decay_config['subsets'])

    modes = {'individual': ['']}
    if 'modes' in decay_config:
        modes.update(decay_config['modes'])

    params = {'step': 10, 'iteration': 10}
    if 'params' in decay_config:
        params.update(decay_config['params'])

    return subsets, modes, params


def run_distance_decay(i_datasets_folder: str, betas: dict, p_distance_decay: str, 
                       datasets_rarefs: dict, force: bool, prjct_nm: str, qiime_env: str,
                       chmod: str, noloc: bool, slurm: bool, split: bool,
                       run_params: dict,
                       filt_raref: str, jobs: bool, chunkit: int) -> (dict, list):

    job_folder2 = get_job_folder(i_datasets_folder, 'decay/chunks')
    decay_config = read_yaml_file(p_distance_decay)
    subsets, modes, params = get_decay_config(decay_config)

    all_sh_pbs = {}
    decay_res = {}
    for dat, rarefs_metrics_groups_metas_qzas_dms_trees in betas.items():
        if not split:
            out_sh = '%s/run_decay_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
        decay_res[dat] = []
        for idx, metrics_groups_metas_qzas_dms_trees in enumerate(rarefs_metrics_groups_metas_qzas_dms_trees):
            decay_raref = {}
            cur_raref = datasets_rarefs[dat][idx]
            odir = get_analysis_folder(i_datasets_folder, 'decay/%s%s' % (dat, cur_raref))
            if split:
                out_sh = '%s/run_decay_%s_%s%s%s.sh' % (job_folder2, prjct_nm, dat, cur_raref, filt_raref)
            for metric, groups_metas_qzas_dms_trees in metrics_groups_metas_qzas_dms_trees.items():
                for group, metas_qzas_mat_qzas_trees in groups_metas_qzas_dms_trees.items():
                    for (meta, qza, mat_qza, tree) in metas_qzas_mat_qzas_trees:
                        meta_pd = read_meta_pd(meta).set_index('sample_name')
                        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(subsets), 'decay')
                        for case_var, case_vals_list in cases_dict.items():
                            for case_vals in case_vals_list:
                                case = get_case(case_vals, case_var).replace(' ', '_')
                                cur_sh = '%s/run_decay_%s%s_%s_%s_%s%s.sh' % (
                                    job_folder2, dat, cur_raref, metric, group, case, filt_raref)
                                cur_sh = cur_sh.replace(' ', '-')
                                all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                                new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
                                res = run_single_decay(odir, group, new_meta_pd, cur_sh, mat_qza,
                                                       case, modes, force, run_params["n_nodes"],
                                                       run_params["n_procs"], int(params['iteration']),
                                                       int(params['step']))
                                decay_raref[(metric, group, case)] = res
            decay_res[dat].append(decay_raref)

    job_folder = get_job_folder(i_datasets_folder, 'decay')
    main_sh = write_main_sh(job_folder, '3_run_decay_%s%s' % (prjct_nm, filt_raref), all_sh_pbs,
                            '%s.prm%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, slurm, jobs, chunkit)
    if main_sh:
        if p_distance_decay:
            print("# decay (config in %s)" % p_distance_decay)
        else:
            print("# decay")
        print_message('', 'sh', main_sh, jobs)

    return decay_res


def run_single_decay(odir: str, group: str, new_meta_pd: pd.DataFrame,
                     cur_sh: str, mat_qza: str, case: str, modes: dict,
                     force: bool, n_nodes: str, n_procs: str,
                     iteration: int, step: int) -> dict:
    res = {}
    remove = True
    mode_group_source = ''
    mode_group_target = ''
    mode_value_source = ''
    mode_value_target = ''
    with open(cur_sh, 'w') as cur_sh_o:
        for mode, mode_groups in modes.items():
            for mode_group in mode_groups:
                if mode == 'individual':
                    cur_rad = '%s/%s/%s_%s_%s' % (
                        odir, mode, splitext(basename(mat_qza))[0], group, case)
                    mode_meta_pd = new_meta_pd.iloc[:, :2].reset_index()
                elif 'targeted' in mode:
                    continue
                else:
                    if group:
                        cur_rad = '%s/%s_%s/%s_%s_%s' % (
                            odir, mode, mode_group, splitext(basename(mat_qza))[0], group, case)
                    else:
                        cur_rad = '%s/%s_%s/%s_%s' % (
                            odir, mode, mode_group, splitext(basename(mat_qza))[0], case)
                    if mode_group not in set(new_meta_pd.columns):
                        continue
                    if new_meta_pd[mode_group].unique().size == 1:
                        continue
                    if min(new_meta_pd[mode_group].value_counts()) < 25:
                        continue
                    if str(new_meta_pd[mode_group].dtype) != 'object':
                        continue
                    mode_meta_pd = new_meta_pd[[mode_group]].reset_index()

                if not isdir(cur_rad):
                    os.makedirs(cur_rad)

                new_meta = '%s.meta' % cur_rad
                mode_meta_pd.columns = ['#SampleID'] + mode_meta_pd.columns.tolist()[1:]
                mode_meta_pd.to_csv(new_meta, index=False, sep='\t')
                mat_qza_filt = '%s.qza' % cur_rad
                new_qza = '%s_decay.qza' % cur_rad
                new_tsv = '%s_decay.tsv' % cur_rad
                if not isfile(new_tsv) or force:
                    write_distance_decay(mat_qza, mat_qza_filt, new_qza, new_tsv,
                                         new_meta, mode, mode_group, mode_group_source,
                                         mode_group_target, mode_value_source,
                                         mode_value_target, iteration, step,
                                         n_nodes, n_procs, cur_sh_o)
                    remove = False
                res[mode, mode_group] = new_tsv
    if remove:
        os.remove(cur_sh)
    return res


def distance_decay_figure(i_datasets_folder: str,
                          distance_decay_res: dict,
                          datasets_rarefs: dict,
                          filt_raref: str) -> None:

    odir = get_analysis_folder(i_datasets_folder, 'decay')
    for dat, metrics_group_subset_mode_modegrps in distance_decay_res.items():
        not_exists_yet = 0
        decays = []
        rarefs = set()
        groups = set()
        subsets = set()
        modes = set()
        metrics = set()
        for idx, metrics_group_subset_mode_modegrp in enumerate(metrics_group_subset_mode_modegrps):
            raref = datasets_rarefs[dat][idx]
            if raref:
                cur_raref = raref.split('_raref')[-1]
            else:
                cur_raref = 'raw'
            for (metric, group, subset), mode_modegrp in metrics_group_subset_mode_modegrp.items():
                metrics.add(metric)
                groups.add(group)
                subsets.add(subset)
                for (mode, mode_group), fp in mode_modegrp.items():

                    if not isfile(fp):
                        not_exists_yet += 1
                        continue

                    if mode_group:
                        cur_mode = '%s (%s)' % (mode, mode_group)
                    else:
                        cur_mode = mode
                    modes.add(cur_mode)

                    decay = pd.read_csv(fp)
                    decay['rarefaction'] = cur_raref
                    decay['metric'] = metric
                    decay['metric / rarefaction'] = decay['metric'] + ' / ' + decay['rarefaction']
                    decay['features group'] = group
                    decay['samples subset'] = subset
                    decay['analysis mode'] = cur_mode

                    decays.append(decay)

        decays_pd = pd.concat(decays)
        if 'aitchison' in set(decays_pd.metric):
            decays_pds = {
                '_aitchison': decays_pd.loc[decays_pd.metric == 'aitchison'],
                '': decays_pd.loc[decays_pd.metric != 'aitchison']
            }
        else:
            decays_pds = {'': decays_pd.copy()}

        for aitchison, decays_pd in decays_pds.items():
            row = 'features group'
            col = 'samples subset'
            if groups == {''}:
                row = ''
            if subsets == {'ALL'}:
                col = ''

            title = dat

            style = 'analysis mode'
            hue = 'metric / rarefaction'
            if rarefs == {'raw'}:
                title += ' - no rarefaction'
                hue = 'metric'
            if modes == {'analysis mode'}:
                if col:
                    style = col
                    col = ''
                elif row:
                    style = row
                    row = ''

            if row:
                height = decays_pd[row].unique().size * 3
            else:
                height = 4

            if row and col:
                g = sns.relplot(data=decays_pd, x='step', y='min', hue=hue,
                                style=style, row=row, col=col, kind='line',
                                facet_kws={'sharex': False}, height=height)
                g.set_titles('{row_name}\n{col_name}')
            elif row:
                g = sns.relplot(data=decays_pd, x='step', y='min', hue=hue,
                                style=style, col=row, kind='line', col_wrap=3,
                                facet_kws={'sharex': False}, height=height)
                g.set_titles('{col_name}')
            elif col:
                g = sns.relplot(data=decays_pd, x='step', y='min', hue=hue,
                                style=style, col=col, kind='line', col_wrap=3,
                                facet_kws={'sharex': False}, height=height)
                g.set_titles('{col_name}')
            else:
                g = sns.relplot(data=decays_pd, x='step', y='min', hue=hue,
                                style=style, kind='line', height=height)
            plt.suptitle(title, fontsize=12)
            plt.subplots_adjust(top=0.93)
            fig_o = '%s/%s_decays%s.pdf' % (odir, dat, aitchison)
            plt.savefig(fig_o, bbox_inches='tight')
            print('    (decay) Written figure: %s' % fig_o)
            if not_exists_yet:
                print('            -> it has %s decay analyses missing (need [TO RUN])' % not_exists_yet)

