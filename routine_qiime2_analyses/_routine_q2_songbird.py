# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import random
import itertools
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_songbird_dicts,
    write_main_sh,
    read_meta_pd,
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict,
    check_metadata_models,
    rename_duplicate_columns
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_songbird_cmd
)
from routine_qiime2_analyses._routine_q2_mmbird import get_mmvec_outputs
from routine_qiime2_analyses._routine_q2_mmvec import (
    make_filtered_and_common_dataset,
    check_filtered_and_common_dataset
)


def get_train_column(new_meta_pd, meta_vars, train, new_meta_ct):
    if train.isdigit() or train.replace('.', '').isdigit():
        train_column = 'TrainTest'
        if train.isdigit():
            train_int = int(train)
            if train_int < (0.1 * new_meta_pd.shape[0]):
                train_perc = 0.1
            else:
                train_perc = train_int / new_meta_pd.shape[0]
        else:
            train_float = float(train)
            if 0 < train_float < 1:
                train_perc = train_float
            else:
                raise IOError('\t\t\t[SONGBIRD] Float passed as percent of samples for'
                              ' training not valid (must be in range 0-1)')
        new_meta_pd['concat_cols'] = new_meta_pd[meta_vars].apply(
            func=lambda x: '__'.join([str(y) for y in x]), axis=1)

        vc = new_meta_pd['concat_cols'].value_counts()
        # if len(meta_vars) == 1 and str(new_meta_pd['concat_cols'].dtype) == 'object' and min(vc) > 1:
        X = np.array(new_meta_pd.values)
        y = new_meta_pd.index.tolist()
        _, __, test_samples, train_samples = train_test_split(
            X, y, test_size=train_perc,
            stratify=new_meta_pd['concat_cols'].tolist()
        )
        # else:
        #     train_samples = random.sample(
        #         new_meta_pd.index.tolist(),
        #         k=int(train_perc * new_meta_pd.shape[0])
        #     )
        new_meta_pd[train_column] = ['Train' if x in train_samples else
                                     'Test' for x in new_meta_pd.index]
        ct = pd.crosstab(new_meta_pd[train_column], new_meta_pd['concat_cols']).T
        ct.to_csv(new_meta_ct, sep='\t')
    else:
        if train in new_meta_pd.columns:
            if {'Train', 'Test'}.issubset(new_meta_pd[train]):
                train_column = train
                new_meta_pd = new_meta_pd.loc[
                    new_meta_pd[train].isin(['Train', 'Test'])
                ]
            else:
                raise IOError('\t\t\t[SONGBIRD] Columns passed for training do '
                              'not have "Train" and "Test" factors')
        else:
            raise IOError('\t\t\t[SONGBIRD] Columns passed for training not exists')
    return new_meta_pd, train_column


def run_single_songbird(odir: str, odir_base: str, qza: str, new_qza: str,
                        new_meta: str, cur_sh: str, force: bool, batch: str,
                        learn: str, epoch: str, diff_prior: str, thresh_feat: str,
                        thresh_sample: str, formula: str, train_column: str, metadatas: dict,
                        baselines: dict, model_baseline: str, baseline_formula: str) -> (str, str):
    """
    Run songbird: Vanilla regression methods for microbiome differential abundance analysis.
    https://github.com/biocore/songbird
    (in-loop function).
    :param odir: output analysis directory.
    :param qza: features table input.
    :param meta_pd: metadata table.
    :param cur_sh: input bash script file.
    :param case:
    :param formula:
    :param case_var:
    :param case_vals:
    :param force:
    :param batch:
    :param learn:
    :param epoch:
    :param diff_prior:
    :param thresh_feat:
    :param thresh_sample:
    :param train:
    :param force: Force the re-writing of scripts for all commands.
    """
    remove = True
    diffs = '%s/differentials.tsv' % odir
    diffs_qza = '%s/differentials.qza' % odir
    stats = '%s/differentials-stats.qza' % odir
    plot = '%s/differentials-biplot.qza' % odir
    if model_baseline in baselines:
        base_diff_qza = ''
        base_stats = baselines[model_baseline]
        base_plot = ''
    else:
        base_diff_qza = '%s/differentials-baseline.qza' % odir_base
        base_stats = '%s/differentials-stats-baseline.qza' % odir_base
        base_plot = '%s/differentials-biplot-baseline.qza' % odir_base
        baselines[model_baseline] = base_stats
    tensor = '%s/tensorboard.qzv' % odir_base
    tensor_html = '%s/tensorboard.html' % odir_base
    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(tensor_html):
            write_songbird_cmd(
                qza, new_qza, new_meta, formula, epoch, batch, diff_prior,
                learn, thresh_sample, thresh_feat, train_column, metadatas,
                diffs, diffs_qza, stats, plot, base_diff_qza, base_stats,
                base_plot, baseline_formula, tensor, tensor_html, cur_sh_o)
            remove = False
    if remove:
        os.remove(cur_sh)
    return diffs, tensor_html


# def get_songbird_metadata_train_test(meta_pd, meta_vars_, meta_var, new_meta,
#                                      train, case, case_var, case_vals, drop):
#     new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
#     if train in new_meta_pd.columns:
#         meta_vars_.append(train)
#     new_meta_pd.columns = [x.lower() for x in new_meta_pd.columns]
#     if meta_var:
#         meta_vars = list(set(list(meta_vars_) + [meta_var]))
#     else:
#         meta_vars = list(meta_vars_)
#     new_meta_pd = new_meta_pd[meta_vars]
#     new_meta_pd = new_meta_pd.loc[~new_meta_pd.isna().any(1)]
#     new_meta_pd = rename_duplicate_columns(new_meta_pd)
#     if len(drop):
#         new_meta_pd = new_meta_pd.loc[(~new_meta_pd[meta_var.lower()].isin(drop)), :]
#     new_meta_pd_ = new_meta_pd.copy()
#     new_meta_pd_['tmptmptmp'] = [''.join(map(str, x)) for x in new_meta_pd_.values if str(x) != 'nan']
#     if 1 in new_meta_pd_.tmptmptmp.value_counts():
#         return None
#
#     new_meta_pd, train_column = get_train_column(new_meta_pd, meta_vars, train)
#     new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
#     return train_column

def get_metadata_train_test(meta_pd, meta_vars, new_meta, train, drop, new_meta_ct):

    if train in meta_pd.columns:
        meta_vars.append(train)

    new_meta_pd = meta_pd[meta_vars]
    new_meta_pd = new_meta_pd.loc[~new_meta_pd.isna().any(1)]
    new_meta_pd = rename_duplicate_columns(new_meta_pd)

    if drop:
        to_remove = pd.concat([
            new_meta_pd[meta_var.lower()].isin(var_drop) for meta_var, var_drop in drop.items()],
            axis=1
        ).any(axis=1)
        new_meta_pd = new_meta_pd.loc[~to_remove]

    new_meta_pd_ = new_meta_pd.copy()
    new_meta_pd_['tmptmptmp'] = [''.join(map(str, x)) for x in new_meta_pd_.values if str(x) != 'nan']
    if 1 in new_meta_pd_.tmptmptmp.value_counts():
        return None
    new_meta_pd, train_column = get_train_column(new_meta_pd, meta_vars, train, new_meta_ct)
    new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
    return train_column


def get_unique_filterings(songbird_filtering):
    unique_filterings = {}
    for filt_name, dat_d in songbird_filtering[''].items():
        for dat, (preval, abund) in dat_d.items():
            unique_filterings.setdefault(dat, set()).add((filt_name, preval, abund))
    return unique_filterings


# def run_songbird(p_diff_models: str, i_datasets_folder: str, datasets: dict,
#                  datasets_read: dict, datasets_filt: dict, input_to_filtered: dict,
#                  mmvec_outputs: list, force: bool, prjct_nm: str, qiime_env: str,
#                  chmod: str, noloc: bool, split: bool,
#                  run_params: dict, filt_raref: str) -> list:
#     """
#     Run songbird: Vanilla regression methods for microbiome differential abundance analysis.
#     https://github.com/biocore/songbird
#     Main per-dataset looper for the songbird datasets.
#
#     :param p_diff_models: Formulas for multinomial regression-based differential abundance ranking.
#     :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
#     :param datasets: list of datasets.
#     :param force: Force the re-writing of scripts for all commands.
#     :param prjct_nm: Nick name for your project.
#     :param qiime_env: qiime2-xxxx.xx conda environment.
#     :param chmod: whether to change permission of output files (default: 775).
#     """
#     job_folder = get_job_folder(i_datasets_folder, 'songbird')
#     job_folder2 = get_job_folder(i_datasets_folder, 'songbird/chunks')
#     songbird_dicts = get_songbird_dicts(p_diff_models)
#     songbird_models = songbird_dicts[0]
#     songbird_filtering = songbird_dicts[1]
#     unique_filtering = get_unique_filterings(songbird_filtering)
#
#     params = songbird_dicts[2]
#     models_baselines = songbird_dicts[3]
#     songbird_datasets = songbird_dicts[4]
#     main_cases_dict = songbird_dicts[5]
#
#     trains = params['train']
#     batches = params['batches']
#     learns = params['learns']
#     epochs = params['epochs']
#     thresh_feats = params['thresh_feats']
#     thresh_samples = params['thresh_samples']
#     diff_priors = params['diff_priors']
#
#     filt_datasets_done, common_datasets_done = check_filtered_and_common_dataset(
#         i_datasets_folder, datasets, datasets_filt, songbird_datasets,
#         {}, songbird_filtering, unique_filtering,
#         'songbird', input_to_filtered, {})
#
#     already_computed = {}
#     filt_datasets, common_datasets = make_filtered_and_common_dataset(
#         i_datasets_folder, datasets, datasets_filt,
#         datasets_read, songbird_datasets, {}, songbird_filtering,
#         unique_filtering, job_folder, force, prjct_nm, qiime_env,
#         chmod, noloc, 'songbird', filt_raref, filt_datasets_done,
#         common_datasets_done, input_to_filtered, already_computed)
#
#     songbird_models.update(dict((input_to_filtered[x], y)
#         for x, y in songbird_models.items() if x in input_to_filtered))
#
#     songbirds = {}
#     for dat, filts_files in filt_datasets.items():
#         for filts, files in filts_files.items():
#             songbirds.setdefault(dat[0], []).append([filts, files[0], files[2], ''])
#
#     if mmvec_outputs:
#         mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
#         mmvec_outputs_pd.to_csv('~/testtest.txt', index=False, sep='\t')
#         for r, row in mmvec_outputs_pd.iterrows():
#             pair = row['pair']
#             omic1 = row['omic1']
#             omic2 = row['omic2']
#             filt1 = row['filt1']
#             filt2 = row['filt2']
#             omic1_common_fp = row['omic1_common_fp']
#             omic2_common_fp = row['omic2_common_fp']
#             meta_common_fp = row['meta_common_fp']
#             songbirds.setdefault(omic1, []).append([filt1, omic1_common_fp, meta_common_fp, pair])
#             songbirds.setdefault(omic2, []).append([filt2, omic2_common_fp, meta_common_fp, pair])
#
#     uni = 0
#
#     all_sh_pbs = {}
#     first_print = 0
#     songbird_outputs = []
#     for dat, filts_tsvs_metas_pair in songbirds.items():
#         if not split:
#             out_sh = '%s/run_songbird_%s%s.sh' % (job_folder2, dat, filt_raref)
#         for (filt, tsv, meta_, pair) in filts_tsvs_metas_pair:
#             if split:
#                 if pair:
#                     out_sh = '%s/run_songbird_%s_%s_%s%s.sh' % (job_folder2, dat, filt, pair, filt_raref)
#                 else:
#                     out_sh = '%s/run_songbird_%s_%s%s.sh' % (job_folder2, dat, filt, filt_raref)
#
#             meta_alphas = '%s_alphas_full.tsv' % splitext(meta_)[0]
#             if isfile(meta_alphas):
#                 meta = meta_alphas
#             else:
#                 meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
#                 if isfile(meta_alphas):
#                     meta = meta_alphas
#                 else:
#                     meta = meta_
#                     if not first_print:
#                         print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
#                               '\t(if you have alpha diversity as a factors in the models)!')
#                         first_print += 1
#             if pair:
#                 dat_pair = '%s_%s' % (dat, pair)
#                 dat_pair_path = '%s/%s' % (dat, pair)
#             else:
#                 dat_pair = '%s' % dat
#                 dat_pair_path = '%s' % dat
#
#             qza = '%s.qza' % splitext(tsv)[0]
#             meta_pd = read_meta_pd(meta)
#             meta_pd = rename_duplicate_columns(meta_pd)
#             meta_pd = meta_pd.set_index('sample_name')
#
#             if dat in songbird_models:
#                 models = check_metadata_models(meta, meta_pd, songbird_models[dat])
#             else:
#                 continue
#
#             cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'songbird')
#             for case_var, case_vals_list in cases_dict.items():
#                 for case_vals in case_vals_list:
#
#                     #####################################################################
#                     # snakemake here: config to organise the inputs/depedencies (joblib)
#                     #####################################################################
#
#                     case = get_case(case_vals, case_var)
#                     for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors,
#                                                                thresh_feats, thresh_samples, trains)):
#                         batch, learn, epoch, diff_prior, thresh_feat, thresh_sample, train = [str(x) for x in it]
#                         params = 'filt_f%s_s%s/%s_%s_%s_%s_%s' % (
#                             thresh_feat, thresh_sample, batch, learn, epoch,
#                             diff_prior.replace('.', ''), train.replace('.', '') )
#                         if uni:
#                             pass
#                             # TO DEVELOP: RUN ALL MODELS BASED ON THE SAME SET OF TESTTRAIN SAMPLES
#                             # uni_meta_vars, unil_meta_var, uni_drop = set(), set(), set()
#                             # uni_datdir = '%s/%s/%s/%s' % (dat_pair_path, filt, case, params)
#                             # uni_datdir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % uni_datdir)
#                             # uni_models = {}
#                             # uni_model_baselines = {}
#                             # for model, formula_meta_var_drop in models.items():
#                             #     uni_model_baselines[model] = {'1': '"1"'}
#                             #     formula, meta_vars, meta_var, drop = formula_meta_var_drop
#                             #     uni_meta_vars.update(set(meta_vars))
#                             #     unil_meta_var.add(meta_var)
#                             #     uni_drop.update(set(drop))
#                             #     uni_models[model] = formula
#                             #     if dat in models_baselines and model in models_baselines[dat]:
#                             #         uni_model_baselines[model].update(models_baselines[dat][model])
#                             #
#                             # print(uni_meta_vars, unil_meta_var, uni_drop)
#                             # print(uni_datdir)
#                             # print(uni_models)
#                             # print(uni_model_baselines)
#                             # uni_new_qza = '%s/uni_tab.qza' % uni_datdir
#                             # uni_new_meta = '%s/uni_metadata.tsv' % uni_datdir
#                             # train_column = get_songbird_metadata_train_test(
#                             #     meta_pd, meta_vars, meta_var, new_meta,
#                             #     train, case, case_var, case_vals, drop)
#                             #
#                             # uni_baselines = {}
#                             # uni_metadata = {}
#                             # print(gfds)
#
#                         for model, formula_meta_var_drop in models.items():
#
#                             datdir = '%s/%s/%s/%s/%s' % (dat_pair_path, filt, case, params, model)
#                             odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % datdir)
#                             new_qza = '%s/tab.qza' % odir
#                             new_meta = '%s/metadata.tsv' % odir
#
#                             formula, meta_vars, meta_var, drop = formula_meta_var_drop
#                             train_column = get_songbird_metadata_train_test(
#                                 meta_pd, meta_vars, meta_var, new_meta,
#                                 train, case, case_var, case_vals, drop)
#                             if not train_column:
#                                 new_meta_invalid = '%s/metadata_invalid' % odir
#                                 with open(new_meta_invalid, 'w') as invalid:
#                                     pass
#                                 continue
#
#                             baselines = {}
#                             metadatas = {}
#                             model_baselines = {'1': '"1"'}
#                             if dat in models_baselines and model in models_baselines[dat]:
#                                 model_baselines = models_baselines[dat][model]
#
#                             for model_baseline, baseline_formula in model_baselines.items():
#                                 odir_base = get_analysis_folder(i_datasets_folder, 'songbird/%s/b-%s' % (datdir, model_baseline))
#
#                                 cur_sh = '%s/run_songbird_%s_%s_%s_%s_%s_%s.sh' % (
#                                     job_folder2, dat_pair, filt, case, model, model_baseline, idx)
#                                 cur_sh = cur_sh.replace(' ', '-')
#                                 all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
#
#                                 diffs, tensor_html = run_single_songbird(
#                                     odir, odir_base, qza, new_qza, new_meta, cur_sh,
#                                     force, batch, learn, epoch, diff_prior, thresh_feat, thresh_sample,
#                                     formula, train_column, metadatas, baselines, model_baseline, baseline_formula
#                                 )
#                                 songbird_outputs.append([dat, filt, '%s_%s' % (params.replace('/', '__'), model), case,
#                                                          diffs, model_baseline, tensor_html, pair])
#
#     job_folder = get_job_folder(i_datasets_folder, 'songbird')
#     main_sh = write_main_sh(job_folder, '2_songbird%s' % filt_raref, all_sh_pbs,
#                             '%s.sngbrd%s' % (prjct_nm, filt_raref),
#                             run_params["time"], run_params["n_nodes"], run_params["n_procs"],
#                             run_params["mem_num"], run_params["mem_dim"],
#                             qiime_env, chmod, noloc)
#     if main_sh:
#         if p_diff_models.startswith('/panfs'):
#             p_diff_models = p_diff_models.replace(os.getcwd(), '')
#         print_message("# Songbird (configs in %s)" % p_diff_models, 'sh', main_sh)
#     return songbird_outputs


def run_songbird(p_diff_models: str, i_datasets_folder: str, datasets: dict,
                 datasets_read: dict, datasets_filt: dict, input_to_filtered: dict,
                 mmvec_outputs: list, force: bool, prjct_nm: str, qiime_env: str,
                 chmod: str, noloc: bool, split: bool, run_params: dict,
                 filt_raref: str, jobs: bool, chunkit: int) -> list:
    """
    Run songbird: Vanilla regression methods for microbiome differential abundance analysis.
    https://github.com/biocore/songbird
    Main per-dataset looper for the songbird datasets.

    :param p_diff_models: Formulas for multinomial regression-based differential abundance ranking.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (default: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'songbird')
    job_folder2 = get_job_folder(i_datasets_folder, 'songbird/chunks')
    songbird_dicts = get_songbird_dicts(p_diff_models)
    songbird_models = songbird_dicts[0]
    songbird_filtering = songbird_dicts[1]
    unique_filtering = get_unique_filterings(songbird_filtering)

    params = songbird_dicts[2]
    models_baselines = songbird_dicts[3]
    songbird_datasets = songbird_dicts[4]
    songbird_subsets = songbird_dicts[5]

    trains = params['train']
    batches = params['batches']
    learns = params['learns']
    epochs = params['epochs']
    thresh_feats = params['thresh_feats']
    thresh_samples = params['thresh_samples']
    diff_priors = params['diff_priors']

    filt_datasets_done, common_datasets_done = check_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_filt, songbird_datasets,
        {}, songbird_filtering, unique_filtering,
        'songbird', input_to_filtered, songbird_subsets)

    already_computed = {}
    filt_datasets, common_datasets = make_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_filt,
        datasets_read, songbird_datasets, {}, songbird_filtering,
        unique_filtering, job_folder, force, prjct_nm, qiime_env,
        chmod, noloc, 'songbird', filt_raref, filt_datasets_done,
        common_datasets_done, input_to_filtered, already_computed,
        songbird_subsets, jobs)

    songbird_models.update(dict((input_to_filtered[x], y)
        for x, y in songbird_models.items() if x in input_to_filtered))

    songbirds = {}
    for dat, filts_files in filt_datasets.items():
        for (case, filts), files in filts_files.items():
            songbirds.setdefault(dat[0], []).append([case, filts, files[0], files[2], ''])

    if mmvec_outputs:
        mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
        mmvec_outputs_pd.to_csv('~/testtest.txt', index=False, sep='\t')
        for r, row in mmvec_outputs_pd.iterrows():
            pair = row['pair']
            case = row['case']
            omic1 = row['omic1']
            omic2 = row['omic2']
            filt1 = row['filt1']
            filt2 = row['filt2']
            omic1_common_fp = row['omic1_common_fp']
            omic2_common_fp = row['omic2_common_fp']
            meta_common_fp = row['meta_common_fp']
            songbirds.setdefault(omic1, []).append([case, filt1, omic1_common_fp, meta_common_fp, pair])
            songbirds.setdefault(omic2, []).append([case, filt2, omic2_common_fp, meta_common_fp, pair])

    uni = 0

    all_sh_pbs = {}
    first_print = 0
    songbird_outputs = []
    for dat, case_filts_tsvs_metas_pair in songbirds.items():
        if not split:
            out_sh = '%s/run_songbird_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
        for (case, filt, tsv, meta_, pair) in case_filts_tsvs_metas_pair:
            if split:
                if pair:
                    out_sh = '%s/run_songbird_%s_%s_%s_%s_%s%s.sh' % (job_folder2,prjct_nm, dat, case,
                                                                      filt, pair, filt_raref)
                else:
                    out_sh = '%s/run_songbird_%s_%s_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, case, filt, filt_raref)

            meta_alphas = '%s_alphas_full.tsv' % splitext(meta_)[0]
            if isfile(meta_alphas):
                meta = meta_alphas
            else:
                meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
                if isfile(meta_alphas):
                    meta = meta_alphas
                else:
                    meta = meta_
                    if not first_print:
                        print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                              '\t(if you have alpha diversity as a factors in the models)!')
                        first_print += 1
            if pair:
                dat_pair = '%s_%s' % (dat, pair)
                dat_pair_path = '%s/%s' % (dat, pair)
            else:
                dat_pair = dat
                dat_pair_path = dat

            qza = '%s.qza' % splitext(tsv)[0]
            meta_pd = read_meta_pd(meta)
            meta_pd = rename_duplicate_columns(meta_pd)
            meta_pd = meta_pd.set_index('sample_name')
            meta_pd.columns = [x.lower() for x in meta_pd.columns]

            if dat in songbird_models:
                models = check_metadata_models(meta, meta_pd, songbird_models[dat])
            else:
                continue

            #####################################################################
            # snakemake here: config to organise the inputs/depedencies (joblib)
            #####################################################################
            for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors,
                                                       thresh_feats, thresh_samples, trains)):
                batch, learn, epoch, diff_prior, thresh_feat, thresh_sample, train = [str(x) for x in it]
                params = 'filt_f%s_s%s/%s_%s_%s_%s_%s' % (
                    thresh_feat, thresh_sample, batch, learn, epoch,
                    diff_prior.replace('.', ''), train.replace('.', '') )
                if uni:
                    pass
                    # TO DEVELOP: RUN ALL MODELS BASED ON THE SAME SET OF TESTTRAIN SAMPLES
                    # uni_meta_vars, unil_meta_var, uni_drop = set(), set(), set()
                    # uni_datdir = '%s/%s/%s/%s' % (dat_pair_path, filt, case, params)
                    # uni_datdir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % uni_datdir)
                    # uni_models = {}
                    # uni_model_baselines = {}
                    # for model, formula_meta_var_drop in models.items():
                    #     uni_model_baselines[model] = {'1': '"1"'}
                    #     formula, meta_vars, drop = formula_meta_var_drop
                    #     uni_meta_vars.update(set(meta_vars))
                    #     unil_meta_var.add(meta_var)
                    #     uni_drop.update(set(drop))
                    #     uni_models[model] = formula
                    #     if dat in models_baselines and model in models_baselines[dat]:
                    #         uni_model_baselines[model].update(models_baselines[dat][model])
                    #
                    # print(uni_meta_vars, unil_meta_var, uni_drop)
                    # print(uni_datdir)
                    # print(uni_models)
                    # print(uni_model_baselines)
                    # uni_new_qza = '%s/uni_tab.qza' % uni_datdir
                    # uni_new_meta = '%s/uni_metadata.tsv' % uni_datdir
                    # train_column = get_songbird_metadata_train_test(
                    #     meta_pd, meta_vars, new_meta,
                    #     train, case, case_var, case_vals, drop)
                    #
                    # uni_baselines = {}
                    # uni_metadata = {}
                    # print(gfds)

                for modx, model in enumerate(models.keys()):

                    formula, meta_vars, drop = models[model]
                    # print("meta_pd.shape")
                    # print(meta_pd.shape)
                    # print("meta_pd.columns")
                    # print(meta_pd.columns)
                    # print("meta_vars")
                    # print(meta_vars)
                    # for meta_v in meta_vars:
                    #     print("meta_v")
                    #     print(meta_v)
                    #     print("meta_pd[meta_v].value_counts()")
                    #     print(meta_pd[meta_v].value_counts())

                    datdir = '%s/%s/%s/%s/%s' % (dat_pair_path, filt, case, params, model)
                    odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % datdir)
                    new_qza = '%s/tab.qza' % odir
                    new_meta = '%s/metadata.tsv' % odir
                    new_meta_ct = '%s/metadata_traintest.tsv' % odir

                    train_column = get_metadata_train_test(
                        meta_pd, meta_vars, new_meta, train, drop, new_meta_ct)
                    if not train_column:
                        new_meta_invalid = '%s/metadata_invalid' % odir
                        with open(new_meta_invalid, 'w') as invalid:
                            pass
                        continue

                    baselines = {}
                    metadatas = {}
                    model_baselines = {'1': '"1"'}
                    if dat in models_baselines and model in models_baselines[dat]:
                        model_baselines = models_baselines[dat][model]

                    for mdx, model_baseline in enumerate(model_baselines.keys()):
                        baseline_formula = model_baselines[model_baseline]
                        odir_base = get_analysis_folder(i_datasets_folder, 'songbird/%s/b-%s' % (datdir, model_baseline))

                        cur_sh = '%s/run_songbird_%s_%s_%s_%s_%s_%s.sh' % (
                            job_folder2, dat_pair, filt, case, modx, mdx, idx)
                        cur_sh = cur_sh.replace(' ', '-')
                        all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)

                        diffs, tensor_html = run_single_songbird(
                            odir, odir_base, qza, new_qza, new_meta, cur_sh,
                            force, batch, learn, epoch, diff_prior, thresh_feat, thresh_sample,
                            formula, train_column, metadatas, baselines, model_baseline, baseline_formula
                        )
                        songbird_outputs.append([dat, filt, '%s_%s' % (params.replace('/', '__'), model), case,
                                                 diffs, model_baseline, tensor_html, pair])

    job_folder = get_job_folder(i_datasets_folder, 'songbird')
    main_sh = write_main_sh(job_folder, '2_songbird_%s%s' % (prjct_nm, filt_raref), all_sh_pbs,
                            '%s.sngbrd%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        if p_diff_models.startswith('/panfs'):
            p_diff_models = p_diff_models.replace(os.getcwd(), '')
        print_message("# Songbird (configs in %s)" % p_diff_models, 'sh', main_sh, jobs)
    return songbird_outputs
