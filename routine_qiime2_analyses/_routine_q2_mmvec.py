# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, sys
import itertools
import pandas as pd
from os.path import isfile, splitext
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_mmvec_dicts,
    write_main_sh,
    get_raref_tab_meta_pds
)
from routine_qiime2_analyses._routine_q2_cmds import (
    filter_feature_table,
    run_import, run_export,
    write_mmvec_cmd
)


def filter_mb_table(preval_filt: int, abund_filt: int, tsv_pd: pd.DataFrame) -> pd.DataFrame:
    if preval_filt or abund_filt:
        new_cols = []
        cur_index = tsv_pd.index.tolist()
        cur_columns = tsv_pd.columns.tolist()
        for r, row in tsv_pd.iterrows():
            if sum(row == 0):
                min_thresh = min([x for x in row if x > 0]) * abund_filt
                cur_row = [x if x >= min_thresh else 0 for x in row]
                new_cols.append(cur_row)
            else:
                new_cols.append(row.tolist())
        tsv_pd = pd.DataFrame(
            new_cols,
            index=cur_index,
            columns=cur_columns
        )
        tsv_pd = tsv_pd[tsv_pd.sum(1) > 1]
        n_perc = (preval_filt / tsv_pd.shape[1]) * 100
        if preval_filt and abund_filt:
            tsv_pd = tsv_pd.loc[tsv_pd.astype(bool).sum(1) >= n_perc, :]
        tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
    return tsv_pd


def filter_non_mb_table(preval_filt: int, abund_filt: int, tsv_pd: pd.DataFrame) -> pd.DataFrame:

    n_perc = (preval_filt / tsv_pd.shape[1]) * 100
    if preval_filt and abund_filt:
        tsv_filt_pd = tsv_pd.loc[(tsv_pd.values >= abund_filt).sum(1) >= n_perc, :].copy()
    elif preval_filt:
        tsv_filt_pd = tsv_pd.loc[tsv_pd.values.astype(bool).sum(1) >= n_perc, :].copy()
    elif abund_filt:
        tsv_filt_pd = tsv_pd.loc[tsv_pd.sum(1) > abund_filt, :].copy()
    else:
        tsv_filt_pd = tsv_pd.copy()
    tsv_filt_pd = tsv_filt_pd.loc[:, tsv_filt_pd.sum(0) > 0]
    return tsv_filt_pd


def write_filtered_tsv(tsv_out: str, tsv_pd: pd.DataFrame) -> None:
    tsv_sams_col = tsv_pd.reset_index().columns[0]
    tsv_pd = tsv_pd.reset_index().rename(columns={tsv_sams_col: 'Feature ID'}).set_index('Feature ID')
    print('write_filtered_tsv:', tsv_out)
    tsv_pd.reset_index().to_csv(tsv_out, index=False, sep='\t')


def write_filtered_meta(meta_out: str, meta_pd_: pd.DataFrame, tsv_pd: pd.DataFrame) -> pd.DataFrame:
    meta_filt_pd = meta_pd_.loc[meta_pd_.sample_name.isin(tsv_pd.columns),:].copy()
    # print('write_filtered_meta:', meta_out)
    meta_filt_pd.to_csv(meta_out, index=False, sep='\t')
    return meta_filt_pd


def get_datasets_filtered(i_datasets_folder: str, datasets: dict,
                          datasets_read: dict, unique_datasets: list,
                          mmvec_filtering: dict, force: bool) -> (dict, list):
    """
    Filter the datasets for use in mmvec.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of data_sets.
    :param datasets_read: dataset -> [tsv table, meta table] (here it updates tsv table after features correction)
    :param datasets: datasets names from the yaml pairs.
    :param mmvec_filtering: validated filtering thersholds.
    :param force: Force the re-writing of scripts for all commands.
    :return: list of datasets from filtered threshold.
    """
    filt_jobs = []
    filt_datasets = {}
    for (dat, mb) in unique_datasets:
        dat_dir = get_analysis_folder(i_datasets_folder, 'mmvec/datasets/%s' % dat)
        if datasets_read[dat] == 'raref':
            tsv, meta = datasets[dat]
            if not isfile(tsv):
                print('Must have run rarefaction to use it further...\nExiting')
                sys.exit(0)
            tsv_pd_, meta_pd_ = get_raref_tab_meta_pds(meta, tsv)
            datasets_read[dat] = [tsv_pd_, meta_pd_]
        else:
            tsv_pd_, meta_pd_ = datasets_read[dat]

        dat_filts = {}
        for preval_filt in mmvec_filtering['prevalence']:
            for abund_filt in mmvec_filtering['abundance']:
                # make sure there's no empty row / column
                tsv_pd = tsv_pd_.loc[tsv_pd_.sum(1) > 0, :].copy()
                tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
                if mb:
                    abund_filter = int(abund_filt[1])
                    tsv_pd = filter_mb_table(int(preval_filt), abund_filter, tsv_pd)
                else:
                    abund_filter = int(abund_filt[0])
                    tsv_pd = filter_non_mb_table(int(preval_filt), abund_filter, tsv_pd)
                rad_out = '%s_%s_%s_%ss' % (dat, preval_filt, abund_filter, tsv_pd.shape[1])
                tsv_out = '%s/tab_%s.tsv' % (dat_dir, rad_out)
                tsv_qza = '%s.qza' % splitext(tsv_out)[0]
                meta_out = '%s/meta_%s.tsv' % (dat_dir, rad_out)
                meta_pd = write_filtered_meta(meta_out, meta_pd_, tsv_pd)
                if force or not isfile(tsv_out):
                    write_filtered_tsv(tsv_out, tsv_pd)
                if force or not isfile(tsv_qza):
                    cmd = run_import(tsv_out, tsv_qza, 'FeatureTable[Frequency]')
                    filt_jobs.append(cmd)
                dat_filts[(preval_filt, abund_filter)] = [
                    tsv_qza, meta_out, meta_pd, tsv_pd.columns.tolist()
                ]
        filt_datasets[dat] = dat_filts
    return filt_datasets, filt_jobs


def get_meta_common_sorted(meta: pd.DataFrame, common_sams: list) -> pd.DataFrame:
    meta_subset = meta.loc[meta.sample_name.isin(common_sams),:].copy()
    meta_subset.sort_values('sample_name', inplace=True)
    return meta_subset


def merge_and_write_metas(meta_subset1: pd.DataFrame, meta_subset2: pd.DataFrame, meta_fp: str) -> pd.DataFrame:
    """
    :param meta_subset1:
    :param meta_subset2:
    :param meta_fp:
    :return:
    """
    # get the columns present in both metadata
    common_cols = [x for x in list(set(meta_subset1.columns.tolist()) &
                                   set(meta_subset2.columns.tolist())) if x!='sample_name']
    # get these columns that also have different contents
    diff_cols = [c for c in common_cols if meta_subset1[c].tolist() != meta_subset2[c].tolist()]
    # edit these different columns' names
    if len(diff_cols):
        meta_subset1.rename(columns=dict((c, '%s.1' % c) for c in diff_cols), inplace=True)
        meta_subset2.rename(columns=dict((c, '%s.2' % c) for c in diff_cols), inplace=True)
    meta_subset = meta_subset1.merge(meta_subset2,
        on=(['sample_name'] + [c for c in common_cols if c not in diff_cols]))
    sorting_col =['sample_name'] + [x for x in meta_subset.columns.tolist() if x != 'sample_name']
    meta_subset[sorting_col].to_csv(meta_fp, index=False, sep='\t')
    # print('merge_and_write_metas:', meta_fp)
    return meta_subset


def get_common_datasets(i_datasets_folder: str, mmvec_pairs: dict,
                        filt_datasets: dict, force: bool) -> (dict, list):
    """
    :param i_datasets_folder:
    :param mmvec_pairs:
    :param filt_datasets:
    :param force: Force the re-writing of scripts for all commands.
    :return:
    """
    common_jobs = []
    common_datasets = {}
    for pair, pair_datasets in mmvec_pairs.items():
        data_dir = get_analysis_folder(i_datasets_folder, 'mmvec/common/data/%s' % pair)
        meta_dir = get_analysis_folder(i_datasets_folder, 'mmvec/common/metadata/%s' % pair)
        (omic1, bool1), (omic2, bool2) = pair_datasets
        # pair_datasets
        # e.g. [('dataset_number_3', 1), ('dataset_number_4', 0)]
        # ('dataset_number_3' --> dataset , 1 --> metabolomics == special filtering)
        for fdx, (preval_filt, abund_filter) in enumerate(filt_datasets[omic1]):
            filts_1 = list(filt_datasets[omic1].keys())
            filts_2 = list(filt_datasets[omic2].keys())
            qza1, meta1, meta_pd1, sams1 = filt_datasets[omic1][filts_1[fdx]]
            qza2, meta2, meta_pd2, sams2 = filt_datasets[omic2][filts_2[fdx]]
            common_sams = sorted(set(sams1) & set(sams2))
            meta_subset1 = get_meta_common_sorted(meta_pd1, common_sams)
            meta_subset2 = get_meta_common_sorted(meta_pd2, common_sams)
            meta_fp = '%s/meta_%s__%s_%s_%ss.tsv' % (meta_dir, pair, preval_filt, abund_filter, len(common_sams))
            new_qza1 = '%s/tab_%s__%s__%s_%s_%ss.qza' % (data_dir, omic1, pair, preval_filt, abund_filter, len(common_sams))
            new_qza2 = '%s/tab_%s__%s__%s_%s_%ss.qza' % (data_dir, omic2, pair, preval_filt, abund_filter, len(common_sams))
            new_tsv1 = '%s.tsv' % splitext(new_qza1)[0]
            new_tsv2 = '%s.tsv' % splitext(new_qza2)[0]
            merge_and_write_metas(meta_subset1, meta_subset2, meta_fp)
            if force or not isfile(new_qza1):
                cmd = filter_feature_table(qza1, new_qza1, meta_fp)
                common_jobs.append(cmd)
            if force or not isfile(new_tsv1):
                cmd = run_export(new_qza1, new_tsv1, 'FeatureTable')
                common_jobs.append(cmd)
            if force or not isfile(new_qza2):
                cmd = filter_feature_table(qza2, new_qza2, meta_fp)
                common_jobs.append(cmd)
            if force or not isfile(new_tsv2):
                cmd = run_export(new_qza2, new_tsv2, 'FeatureTable')
                common_jobs.append(cmd)
            common_datasets.setdefault((pair, preval_filt, abund_filter), []).append([meta_fp, new_qza1, new_qza2])
    return common_datasets, common_jobs


def run_multi_mmvec(odir: str, pair: str, meta_fp: str, qza1: str, qza2: str, res_dir: str,
                    cur_sh: str, batch: str, learn: str, epoch: str, prior: str,
                    thresh_feat: str, latent_dim: str, train_column: str,
                    n_example: str, gpu: bool, force: bool, standalone: bool) -> None:
    """
    Run mmvec: Neural networks for microbe-metabolite interaction analysis.
    https://github.com/biocore/mmvec
    (in-loop function).

    :param odir:
    :param pair:
    :param meta_fp:
    :param qza1:
    :param qza2:
    :param res_dir:
    :param cur_sh:
    :param batch:
    :param learn:
    :param epoch:
    :param prior:
    :param thresh_feat:
    :param latent_dim:
    :param train_column:
    :param n_example:
    :param gpu:
    :param standalone:
    :return:
    """
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        cur_rad = '%s/%s' % (odir, pair)
        conditionals_tsv = '%s_conditionals.tsv' % cur_rad
        biplot_tsv = '%s_ordination.tsv' % cur_rad
        if force or not isfile(conditionals_tsv):
            write_mmvec_cmd(meta_fp, qza1, qza2, res_dir,
                            conditionals_tsv, biplot_tsv,
                            batch, learn, epoch, prior,
                            thresh_feat, latent_dim, train_column,
                            n_example, gpu, standalone, cur_sh_o)

            remove = False
    if remove:
        os.remove(cur_sh)


def make_filtered_and_common_dataset(i_datasets_folder:str, datasets: dict,
                                     datasets_read: dict, mmvec_pairs: dict,
                                     mmvec_filtering: dict, job_folder: str,
                                     force: bool, prjct_nm: str, qiime_env: str,
                                     chmod: str) -> (dict, dict):
    """
    :param i_datasets_folder:
    :param datasets: list of data_sets.
    :param datasets_read:
    :param mmvec_pairs:
    :param mmvec_filtering:
    :param job_folder:
    :param force:
    :param prjct_nm:
    :param qiime_env:
    :param chmod:
    :return:
    """
    unique_datasets = list(set([dat for pair_dats in mmvec_pairs.values() for dat in pair_dats]))

    print('\t-> Get datasets filtered...', end = ' ')
    filt_datasets, filt_jobs = get_datasets_filtered(i_datasets_folder, datasets, datasets_read,
                                                     unique_datasets, mmvec_filtering, force)
    print('Done.')
    print('\t-> Get common datasets...', end = ' ')
    common_datasets, common_jobs = get_common_datasets(i_datasets_folder, mmvec_pairs,
                                                       filt_datasets, force)
    print('Done.')
    pre_jobs = filt_jobs + common_jobs
    if len(pre_jobs):
        import_sh = '%s/2_run_mmvec_imports.sh' % job_folder
        import_pbs = '%s.pbs' % splitext(import_sh)[0]
        with open(import_sh, 'w') as import_o:
            for cmd in pre_jobs:
                import_o.write('\necho "%s"\n' % cmd)
                import_o.write('%s\n' % cmd)
        run_xpbs(import_sh, import_pbs, '%s.xprt.lph' % prjct_nm,
                 qiime_env, '2', '1', '1', '150', 'mb', chmod, 1,
                 '# Import common datasets for MMVEC')
    return filt_datasets, common_datasets


def run_mmvec(p_mmvec_pairs: str, i_datasets_folder: str, datasets: dict,
              datasets_read: dict, force: bool, gpu: bool, standalone: bool,
              prjct_nm: str, qiime_env: str, chmod: str) -> dict:
    """
    Run mmvec: Neural networks for microbe-metabolite interaction analysis.
    https://github.com/biocore/mmvec
    Main two-datasets looper for the mmvec co-occurrences.

    :param p_mmvec_pairs: Pairs of datasets for which to compute co-occurrences probabilities.
    :param p_diff_models: Formulas for multinomial regression-based differential abundance ranking.
    :param datasets: list of data_sets.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets_read: dataset -> [tsv table, meta table] (here it updates tsv table after features correction)
    :param force: Force the re-writing of scripts for all commands.
    :param gpu: Use GPUs instead of CPUs for MMVEC.
    :param standalone:
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (default: 644).
    """
    mmvec_pairs, mmvec_filtering, mmvec_params = get_mmvec_dicts(p_mmvec_pairs)
    job_folder = get_job_folder(i_datasets_folder, 'mmvec')
    print(' [mmvec] Make filtered and_common datasets:')
    filt_datasets, common_datasets = make_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_read, mmvec_pairs, mmvec_filtering,
        job_folder, force, prjct_nm, qiime_env, chmod)
    mmvec_outputs = {}

    jobs = []
    all_sh_pbs = {}
    for (pair, preval, abund), meta_qzas in common_datasets.items():

        job_folder2 = get_job_folder(i_datasets_folder, 'mmvec/chunks/%s' % pair)
        out_sh = '%s/chunks/run_mmvec_%s.sh' % (job_folder, pair)
        pair_filt = '_'.join([str(x) for x in [pair, preval, abund]])

        for (meta_fp, qza1, qza2) in meta_qzas:

            print('"""""""')
            print('COMMON:')
            print('""""""""""""""""""""""""""""""')
            print('meta_fp', meta_fp)
            print('qza1', qza1)
            print('qza2', qza2)
            print('""""""""""""""""""""""""""""""')

            train_columns = mmvec_params['train_column']
            n_examples = mmvec_params['n_examples']
            batches = mmvec_params['batches']
            learns = mmvec_params['learns']
            epochs = mmvec_params['epochs']
            priors = mmvec_params['priors']
            thresh_feats = mmvec_params['thresh_feats']
            latent_dims = mmvec_params['latent_dims']
            for idx, it in enumerate(itertools.product(train_columns, n_examples, batches, learns,
                                                       epochs, priors, thresh_feats, latent_dims)):
                train_column, n_example, batch, learn, epoch, prior, thresh_feat, latent_dim = [str(x) for x in it]
                res_dir = 'b-%s_l-%s_e-%s_p-%s_f-%s_d-%s_t-%s_n-%s_gpu-%s' % (
                    batch, learn, epoch, prior.replace('.', ''),
                    thresh_feat, latent_dim, train_column,
                    n_example, str(gpu)[0]
                )
                mmvec_outputs.setdefault(pair, []).append([meta_fp, qza1, qza2, train_column, res_dir])
                odir = get_analysis_folder(i_datasets_folder, 'mmvec/paired/%s/%s_%s/%s' % (pair, preval, abund, res_dir))
                cur_sh = '%s/run_mmvec_%s_%s_%s.sh' % (job_folder2, preval, abund, res_dir)
                all_sh_pbs.setdefault((pair, out_sh), []).append(cur_sh)
                print('[', idx, ']', it)
                p = multiprocessing.Process(
                    target=run_multi_mmvec,
                    args=(odir, pair, meta_fp, qza1, qza2, res_dir, cur_sh,
                          batch, learn, epoch, prior, thresh_feat,
                          latent_dim, train_column, n_example,
                          gpu, force, standalone))
                p.start()
                jobs.append(p)
    for j in jobs:
        j.join()

    main_sh = write_main_sh(job_folder, '3_mmvec', all_sh_pbs,
                            '%s.mmvc' % prjct_nm, '150', '1', '1', '2', 'gb',
                            qiime_env, chmod)
    if main_sh:
        if p_mmvec_pairs.startswith('/panfs'):
            p_mmvec_pairs = p_mmvec_pairs.replace(os.getcwd(), '')
        print_message("# MMVEC (datasets pairs in %s)" % p_mmvec_pairs, 'sh', main_sh)
    return mmvec_outputs