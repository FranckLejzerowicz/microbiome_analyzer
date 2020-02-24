# # ----------------------------------------------------------------------------
# # Copyright (c) 2020, Franck Lejzerowicz.
# #
# # Distributed under the terms of the Modified BSD License.
# #
# # The full license is in the file LICENSE, distributed with this software.
# # ----------------------------------------------------------------------------
#
# import os
# import itertools
# import pandas as pd
# from os.path import basename, isfile, splitext, dirname
# import multiprocessing
#
# from routine_qiime2_analyses._routine_q2_io_utils import (
#     get_job_folder,
#     get_analysis_folder,
#     get_songbird_dict,
#     get_mmvec_dict,
#     write_main_sh,
#     read_meta_pd
# )
# from routine_qiime2_analyses._routine_q2_metadata import (
#     check_metadata_cases_dict,
#     check_metadata_models
# )
# from routine_qiime2_analyses._routine_q2_cmds import (
#     get_new_meta_pd, get_case,
#     write_songbird_cmd
# )
#
#
# def run_multi_mmvec() -> None:
#     """
#     Run mmvec: Neural networks for microbe-metabolite interaction analysis.
#     https://github.com/biocore/mmvec
#     (in-loop function).
#     """
#     remove = True
#
#     cur_rad = odir + '/' + basename(qza).replace('.qza', '_%s' % case)
#     new_meta = '%s.meta' % cur_rad
#     new_qza = '%s.qza' % cur_rad
#     diffs = '%s/differentials.tsv' % cur_rad
#     diffs_qza = '%s/differentials.qza' % cur_rad
#     stats = '%s/differentials-stats.qza' % cur_rad
#     plot = '%s/differentials-biplot.qza' % cur_rad
#     base_diff_qza = '%s/differentials-baseline.qza' % cur_rad
#     base_stats = '%s/differentials-baseline-stats.qza' % cur_rad
#     base_plot = '%s/differentials-baseline-biplot.qza' % cur_rad
#     tensor = '%s/differentials-tensorboard.qzv' % cur_rad
#
#     with open(cur_sh, 'w') as cur_sh_o:
#         if force or not isfile(tensor):
#             new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
#             new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
#             write_songbird_cmd(qza, new_qza, new_meta,
#                                formula, epoch, batch, diff_prior, learn, thresh_sample,
#                                thresh_feat, n_random, diffs, diffs_qza, stats, plot,
#                                base_diff_qza, base_stats, base_plot, tensor, cur_sh_o)
#             remove = False
#     if remove:
#         os.remove(cur_sh)
#
#
# # def make_common(datasets: dict, pair_datasets: list, mmvec_filtering: dict):
#
# def get_mb_flag(dataset: str) -> (str, bool):
#     if dataset.endswith('*'):
#         return dataset[:1], True
#     else:
#         return dataset, False
#
#
# def filter_mb_table(preval_filt: int, abund_filt: int,
#                     tsv_pd: pd.DataFrame) -> pd.DataFrame:
#     if preval_filt or abund_filt:
#         new_cols = []
#         cur_index = tsv_pd.index.tolist()
#         cur_columns = tsv_pd.columns.tolist()
#         for r, row in tsv_pd.iterrows():
#             if sum(row == 0):
#                 min_thresh = min([x for x in row if x > 0]) * abund_filt
#                 cur_row = [x if x >= min_thresh else 0 for x in row]
#                 new_cols.append(cur_row)
#             else:
#                 new_cols.append(row.tolist())
#         tsv_filt_pd = pd.DataFrame(
#             new_cols,
#             index=cur_index,
#             columns=cur_columns
#         )
#         tsv_filt_pd = tsv_filt_pd[tsv_filt_pd.sum(1) > 1]
#         n_perc = (preval_filt / tsv_filt_pd.shape[1]) * 100
#         if preval_filt and abund_filt:
#             tsv_filt_pd = tsv_filt_pd.loc[tsv_filt_pd.astype(bool).sum(1) >= n_perc, :]
#         tsv_filt_pd = tsv_filt_pd.loc[:, tsv_filt_pd.sum(0) > 0]
#         return tsv_filt_pd
#     else:
#         return tsv_pd
#
#
# def filter_non_mb_table(preval_filt: int, abund_filt: int, tsv_pd: pd.DataFrame) -> pd.DataFrame:
#     n_perc = (preval_filt / tsv_pd.shape[1]) * 100
#     if preval_filt and abund_filt:
#         tsv_filt_pd = tsv_pd.loc[(tsv_pd.values >= abund_filt).sum(1) >= n_perc, :].copy()
#     elif preval_filt:
#         tsv_filt_pd = tsv_pd.loc[tsv_pd.values.astype(bool).sum(1) >= n_perc, :].copy()
#     elif abund_filt:
#         tsv_filt_pd = tsv_pd.loc[tsv_pd.sum(1) > abund_filt, :].copy()
#     else:
#         tsv_filt_pd = tsv_pd.copy()
#     tsv_filt_pd = tsv_filt_pd.loc[:, tsv_filt_pd.sum(0) > 0]
#     return tsv_filt_pd
#
#
# def get_datasets_filtered(i_datasets_folder: str, datasets: dict, datasets_read: dict,
#                           unique_datasets: list, mmvec_filtering: dict):
#
#
#     datasets_filts = {}
#     for dataset in unique_datasets:
#         dat, mb = get_mb_flag(dataset)
#
#         job_folder = get_job_folder(i_datasets_folder, 'mmvec/%s' % dat)
#         job_folder2 = get_job_folder(i_datasets_folder, 'mmvec/chunks')
#
#         tsv, meta = datasets[dat]
#         tsv_pd_, meta_pd_ = datasets_read[dat]
#
#         for preval_filt in mmvec_filtering['prevalence']:
#             for abund_filt in mmvec_filtering['abundance']:
#
#                 # ake sure there's no empty row / column
#                 tsv_pd = tsv_pd_.loc[tsv_pd_.sum(1) > 0, :].copy()
#                 tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
#
#                 if mb:
#                     abund_filter = int(abund_filt[1])
#                     tsv_pd = filter_mb_table(int(preval_filt), abund_filter, tsv_pd)
#                 else:
#                     abund_filter = int(abund_filt[0])
#                     tsv_pd = filter_non_mb_table(int(preval_filt), abund_filter, tsv_pd)
#
#                 dat_dir = get_analysis_folder(i_datasets_folder, 'mmvec/%s/datasets' % dat)
#
#                 tsv_out = '%s/%s_%s_%s_%ss.tsv' % (dat_dir, dat, preval_filt, abund_filter, tsv_pd.shape[1])
#                 tsv_pd = tsv_pd.reset_index().rename(
#                     columns={tsv_pd.reset_index().columns[0]: 'Feature ID'}
#                 ).set_index('Feature ID')
#
#
#
#                 if dataset_ != 'mb':
#                     fil_pd.index = [feats_taxa_dict[x] for x in fil_pd.index.tolist()]
#                 fil_pd.reset_index().to_csv(tsv, index=False, sep='\t')
#                 biom = tsv.replace('.tsv', '.biom')
#
#                 meta = tsv.replace('/XX_datasets/', '/00_metadata/').replace('.tsv', '_meta.tsv')
#                 mkdr(meta)
#
#                 tax = tsv.replace('XX_datasets', 'TT_taxonomy').replace('.tsv', '_Tax.tsv')
#                 mkdr(tax)
#
#                 if dataset_ == 'mb':
#                     tax_pd = mb_tax_pd.loc[mb_tax_pd['Feature ID'].isin(fil_pd.index.tolist())].copy()
#                 else:
#                     tax_pd = feats_taxa.loc[feats_taxa['Feature ID'].isin(fil_pd.index.tolist())].copy()
#
#                 cur_meta_pd = meta_pd.loc[meta_pd.sample_name.isin(fil_pd.columns)].copy()
#                 cur_meta_pd.to_csv(meta, index=False, sep='\t')
#
#                 # Leaf_bac_1-Leaf_mb_1-100sams/filt_10/cpu_1000_1e-4_2000/latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.85_beta2_0.90_ordination_meta_Leaf_bac_1.tsv
#                 tax_pd.to_csv(tax, index=False, sep='\t')
#
#                 if dataset_ == 'mb':
#                     k = (dataset, str(preval_filt), str(abund_filt[1]))
#                 else:
#                     k = (dataset, str(preval_filt), str(abund_filt[0]))
#
#                 dataTabs[k] = [
#                     tsv,
#                     fil_pd,
#                     biom,
#                     # qza,
#                     tax,
#                     tax_pd,
#                     meta,
#                     cur_meta_pd
#                 ]
#
#
# def run_mmvec(p_mmvec_pairs: str, p_diff_models: str, i_datasets_folder: str,
#               datasets: dict, datasets_read: dict, force: str, gpu: bool,
#               prjct_nm: str, qiime_env: str, chmod: str) -> None:
#     """
#     Run mmvec: Neural networks for microbe-metabolite interaction analysis.
#     https://github.com/biocore/mmvec
#     Main two-datasets looper for the mmvec co-occurrences.
#
#     :param p_mmvec_pairs: Pairs of datasets for which to compute co-occurrences probabilities.
#     :param p_diff_models: Formulas for multinomial regression-based differential abundance ranking.
#     :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
#     :param datasets: list of datasets.
#     :param datasets_read: dataset -> [tsv table, meta table] (here it updates tsv table after features correction)
#     :param force: Force the re-writing of scripts for all commands.
#     :param gpu: Use GPUs instead of CPUs for MMVEC.
#     :param prjct_nm: Nick name for your project.
#     :param qiime_env: qiime2-xxxx.xx conda environment.
#     :param chmod: whether to change permission of output files (defalt: 775).
#     """
#     job_folder2 = get_job_folder(i_datasets_folder, 'mmvec/chunks')
#     if p_diff_models:
#         songbird_models, songbird_cases_dict, songbird_params = get_songbird_dict(p_diff_models)
#
#     mmvec_pairs, mmvec_filtering, mmvec_params = get_mmvec_dict(p_mmvec_pairs)
#
#     unique_datasets = list(set([dat for pair_dats in mmvec_pairs.values() for dat in pair_dats]))
#     datasets_filts = get_datasets_filtered(datasets, datasets_read, unique_datasets, mmvec_filtering)
#
#     jobs = []
#     all_sh_pbs = {}
#     for pair, pair_datasets in mmvec_pairs.items():
#
#         make_common(datasets, pair_datasets)
#         meta_common_fp = row['meta_common_fp']
#
#         tsvs, metas = [()]
#         for dat, tsv_meta_pds in datasets.items():
#
#         tsv, meta = tsv_meta_pds
#         qza = '%s.qza' % splitext(tsv)[0]
#         meta_pd = read_meta_pd(meta)
#
#
#         cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(songbird_cases_dict), 'mmvec')
#         songbird_models = check_metadata_models(meta, meta_pd, songbird_models, 'mmvec')
#
#         for datasets_pair, datasets, formula in songbird_models.items():
#             out_sh = '%s/run_mmvec_%s_%s.sh' % (job_folder2, dat, model)
#             for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors, thresh_feats, thresh_samples)):
#                 batch, learn, epoch, diff_prior, thresh_feat, thresh_sample = it
#                 res_dir = '%s/%s/filt_f%s_s%s/%s_%s_%s_%s' % (
#                     dat, model.replace('+', 'PLUS').replace('*', 'COMBI').replace('-', 'MINUS').replace('/', 'DIVIDE'),
#                     thresh_feat, thresh_sample, batch, learn, epoch, diff_prior.replace('.', '')
#                 )
#                 odir = get_analysis_folder(i_datasets_folder, 'mmvec/%s' % res_dir)
#                 for case_var, case_vals_list in cases_dict.items():
#                     for case_vals in case_vals_list:
#                         case = get_case(case_vals, model, case_var, str(idx))
#                         cur_sh = '%s/run_mmvec_%s_%s.sh' % (
#                             job_folder2, dat, case)
#                         cur_sh = cur_sh.replace(' ', '-')
#                         all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
#                         p = multiprocessing.Process(
#                             target=run_multi_mmvec,
#                             args=(odir, qza, meta_pd, cur_sh, case, formula,
#                                   batch, learn, epoch, diff_prior, thresh_feat, thresh_sample,
#                                   case_var, case_vals, force))
#                         p.start()
#                         jobs.append(p)
#     for j in jobs:
#         j.join()
#
#     job_folder = get_job_folder(i_datasets_folder, 'mmvec')
#     main_sh = write_main_sh(job_folder, '2_mmvec', all_sh_pbs,
#                             '%s.mmvc' % prjct_nm, '2', '1', '48', '2', 'gb',
#                             qiime_env, chmod)
#     if main_sh:
#         print("# MMVEC (datasets pairs in %s)" % p_mmvec_pairs)
#         print('[TO RUN] sh', main_sh)