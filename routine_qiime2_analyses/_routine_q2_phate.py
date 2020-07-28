# # ----------------------------------------------------------------------------
# # Copyright (c) 2020, Franck Lejzerowicz.
# #
# # Distributed under the terms of the Modified BSD License.
# #
# # The full license is in the file LICENSE, distributed with this software.
# # ----------------------------------------------------------------------------
#
# import os
# import random
# import itertools
# import pandas as pd
# from os.path import isfile, splitext
#
# from routine_qiime2_analyses._routine_q2_xpbs import print_message
# from routine_qiime2_analyses._routine_q2_io_utils import (
#     get_job_folder,
#     get_analysis_folder,
#     get_songbird_dicts,
#     write_main_sh,
#     read_meta_pd,
# )
# from routine_qiime2_analyses._routine_q2_metadata import (
#     check_metadata_cases_dict,
#     check_metadata_models,
#     rename_duplicate_columns
# )
# from routine_qiime2_analyses._routine_q2_cmds import (
#     get_new_meta_pd, get_case,
#     write_phate_cmd
# )
# from routine_qiime2_analyses._routine_q2_mmvec import (
#     make_filtered_and_common_dataset,
#     check_filtered_and_common_dataset
# )
#
#
# def get_train_column(new_meta_pd, train):
#     if train.isdigit() or train.replace('.', '').isdigit():
#         train_column = 'TrainTest'
#         if train.isdigit():
#             train_samples = random.sample(
#                 new_meta_pd.index.tolist(),
#                 k=int(train))
#         else:
#             train_float = float(train)
#             if 0 < train_float < 1:
#                 train_samples = random.sample(
#                     new_meta_pd.index.tolist(),
#                     k=int(train_float * new_meta_pd.shape[0]))
#             else:
#                 raise IOError('Float passed as percent of samples for'
#                               ' training not valid (must be in range 0-1)')
#         new_meta_pd[train_column] = ['Train' if x in train_samples else
#                                      'Test' for x in new_meta_pd.index]
#     else:
#         if train in new_meta_pd.columns:
#             if {'Train', 'Test'}.issubset(new_meta_pd[train]):
#                 train_column = train
#                 new_meta_pd = new_meta_pd.loc[
#                     new_meta_pd[train].isin(['Train', 'Test'])
#                 ]
#             else:
#                 raise IOError('Columns passed for training do '
#                               'not have "Train" and "Test" factors')
#         else:
#             raise IOError('Columns passed for training not exists')
#     return new_meta_pd, train_column
#
#
# def run_single_songbird(odir: str, qza: str, meta_pd: pd.DataFrame, cur_sh: str,
#                         case: str, formula_meta_var_drop: list, case_var: str,
#                         case_vals: list, force: bool, batch: str, learn: str,
#                         epoch: str, diff_prior: str, thresh_feat: str,
#                         thresh_sample: str, train: str, baselines: dict,
#                         baseline: str, baseline_formula: str) -> (str, str):
#     """
#     Run songbird: Vanilla regression methods for microbiome differential abundance analysis.
#     https://github.com/biocore/songbird
#     (in-loop function).
#     :param odir: output analysis directory.
#     :param qza: features table input.
#     :param meta_pd: metadata table.
#     :param cur_sh: input bash script file.
#     :param case:
#     :param formula:
#     :param case_var:
#     :param case_vals:
#     :param force:
#     :param batch:
#     :param learn:
#     :param epoch:
#     :param diff_prior:
#     :param thresh_feat:
#     :param thresh_sample:
#     :param train:
#     :param force: Force the re-writing of scripts for all commands.
#     """
#     remove = True
#
#     new_meta = '%s/metadata.tsv' % odir
#     new_qza = '%s/tab.qza' % odir
#     diffs = '%s/differentials.tsv' % odir
#     diffs_qza = '%s/differentials.qza' % odir
#     stats = '%s/differentials-stats.qza' % odir
#     plot = '%s/differentials-biplot.qza' % odir
#     if baseline in baselines:
#         base_diff_qza = ''
#         base_stats = baselines[baseline]
#         base_plot = ''
#     else:
#         base_diff_qza = '%s/differentials-baseline.qza' % odir
#         base_stats = '%s/differentials-stats-baseline.qza' % odir
#         base_plot = '%s/differentialsbiplot-baseline.qza' % odir
#         baselines[baseline] = base_stats
#     tensor = '%s/tensorboard.qzv' % odir
#     tensor_html = '%s/tensorboard.html' % odir
#
#     formula, vars, meta_var, drop = formula_meta_var_drop
#     with open(cur_sh, 'w') as cur_sh_o:
#         if force or not isfile(tensor_html):
#             new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
#             new_meta_pd.columns = [x.lower() for x in new_meta_pd.columns]
#             if meta_var:
#                 meta_vars = list(set(list(vars) + [meta_var]))
#             else:
#                 meta_vars = list(vars)
#             new_meta_pd = new_meta_pd[meta_vars]
#             new_meta_pd = rename_duplicate_columns(new_meta_pd)
#             if len(drop):
#                 new_meta_pd = new_meta_pd.loc[(~new_meta_pd[meta_var.lower()].isin(drop)), :]
#             new_meta_pd, train_column = get_train_column(new_meta_pd, train)
#             new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
#             write_songbird_cmd(qza, new_qza, new_meta, formula, epoch, batch, diff_prior,
#                                learn, thresh_sample, thresh_feat, train_column, diffs,
#                                diffs_qza, stats, plot, base_diff_qza, base_stats,
#                                base_plot, baseline_formula, tensor, tensor_html, cur_sh_o)
#             remove = False
#     if remove:
#         os.remove(cur_sh)
#     return diffs, tensor_html
#
#
# def get_unique_filterings(songbird_filtering):
#     unique_filterings = {}
#     for filt_name, dat_d in songbird_filtering[''].items():
#         for dat, (preval, abund) in dat_d.items():
#             unique_filterings.setdefault(dat, set()).add((filt_name, preval, abund))
#     return unique_filterings
#
#
# def run_songbird(p_phate_config: str, i_datasets_folder: str, datasets: dict,
#                  datasets_read: dict, datasets_filt: dict, input_to_filtered: dict,
#                  force: bool, prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
#                  split: bool, run_params: dict, filt_raref: str) -> list:
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
#     songbird_dicts = get_songbird_dicts(p_phate_config)
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
#         'songbird', input_to_filtered)
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
#             meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
#             if isfile(meta_alphas):
#                 meta = meta_alphas
#             else:
#                 meta = meta_
#                 if not first_print:
#                     print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
#                           '\t(if you have alpha diversity as a factors in the models)!')
#                     first_print += 1
#             qza = '%s.qza' % splitext(tsv)[0]
#             meta_pd = read_meta_pd(meta)
#             meta_pd = rename_duplicate_columns(meta_pd)
#             meta_pd = meta_pd.set_index('sample_name')
#             cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'songbird')
#             # print(' --->', dat, end=' : ')
#             if dat in songbird_models:
#                 models = check_metadata_models(meta, meta_pd, songbird_models[dat])
#             else:
#                 continue
#             baselines = {}
#             for model, formula_meta_var_drop in models.items():
#                 model_baselines = {'1': '"1"'}
#                 if dat in models_baselines and model in models_baselines[dat]:
#                     model_baselines = models_baselines[dat][model]
#                 for model_baseline, baseline_formula in model_baselines.items():
#                     # print(" ** model, formula_meta_var_drop")
#                     # print(model, formula_meta_var_drop)
#                     for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors,
#                                                                thresh_feats, thresh_samples, trains)):
#                         batch, learn, epoch, diff_prior, thresh_feat, thresh_sample, train = [str(x) for x in it]
#                         model_rep = model.replace('+', 'PLUS').replace('*', 'COMBI').replace('-', 'MINUS').replace('/', 'DIVIDE')
#                         params = '%s/filt_f%s_s%s/%s_%s_%s_%s_%s' % (
#                             model_rep, thresh_feat, thresh_sample,
#                             batch, learn, epoch, diff_prior.replace('.', ''), train.replace('.', '')
#                         )
#                         for case_var, case_vals_list in cases_dict.items():
#                             for case_vals in case_vals_list:
#                                 case = get_case(case_vals, case_var, str(idx))
#                                 baseline = '%s_%s_%s' % (model_baseline, '_'.join([str(x) for x in it]), case)
#                                 if pair:
#                                     res_dir = '%s/%s/%s/%s/%s' % (dat, pair, filt, case, params)
#                                 else:
#                                     res_dir = '%s/%s/%s/%s' % (dat, filt, case, params)
#                                 odir = get_analysis_folder(i_datasets_folder, 'songbird/%s/%s' % (res_dir, model_baseline))
#                                 cur_sh = '%s/run_songbird_%s_%s_%s_%s_%s%s.sh' % (
#                                     job_folder2, dat, filt, model_rep, case, pair, filt_raref)
#                                 cur_sh = cur_sh.replace(' ', '-')
#                                 all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
#                                 diffs, tensor_html = run_single_songbird(
#                                     odir, qza, meta_pd, cur_sh, case, formula_meta_var_drop,
#                                     case_var, case_vals, force, batch, learn,
#                                     epoch, diff_prior, thresh_feat, thresh_sample,
#                                     train, baselines, baseline, baseline_formula
#                                 )
#                                 songbird_outputs.append([dat, filt, params.replace('/', '__'),
#                                                          case, diffs, pair])
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