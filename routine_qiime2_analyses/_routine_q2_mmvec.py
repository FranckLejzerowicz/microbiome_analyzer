# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import itertools
import pandas as pd
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_mmvec_dicts,
    write_main_sh,
    get_datasets_filtered,
    check_datasets_filtered
)
from routine_qiime2_analyses._routine_q2_metadata import rename_duplicate_columns
from routine_qiime2_analyses._routine_q2_cmds import (
    filter_feature_table,
    run_export,
    write_mmvec_cmd
)


def get_meta_common_sorted(meta: pd.DataFrame, common_sams: list) -> pd.DataFrame:
    meta_subset = meta.loc[meta.sample_name.isin(common_sams),:].copy()
    meta_subset.columns = [x.lower() for x in meta_subset.columns]
    meta_subset.sort_values('sample_name', inplace=True)
    return meta_subset


def merge_and_write_metas(meta_subset1: pd.DataFrame,
                          meta_subset2: pd.DataFrame,
                          meta_fp: str) -> pd.DataFrame:
    """
    :param meta_subset1:
    :param meta_subset2:
    :param meta_fp:
    :return:
    """
    meta_subset1 = rename_duplicate_columns(meta_subset1)
    meta_subset2 = rename_duplicate_columns(meta_subset2)

    # get the columns present in both metadata
    common_cols = set(meta_subset1.columns) & set(meta_subset2.columns)
    common_cols = [x for x in common_cols if x != 'sample_name']
    # get these columns that also have different contents
    diff_cols = []
    for c in common_cols:
        try:
            meta_col1 = meta_subset1[c].tolist()
        except:
            print(meta_subset1[c])
            print(fkkj)

        meta_col2 = meta_subset2[c].tolist()

        if meta_col1 != meta_col2:
            diff_cols.append(c)

    if len(diff_cols):
        meta_subset2.rename(columns=dict((c, '%s.copy' % c) for c in diff_cols), inplace=True)
    meta_subset = meta_subset1.merge(meta_subset2,
        on=(['sample_name'] + [c for c in common_cols if c not in diff_cols]))
    sorting_col =['sample_name'] + [x for x in meta_subset.columns.tolist() if x != 'sample_name']
    meta_subset[sorting_col].to_csv(meta_fp, index=False, sep='\t')
    return meta_subset


def get_common_datasets(i_datasets_folder: str, mmvec_pairs: dict, filtering: dict,
                        filt_datasets: dict, common_datasets_done: dict,
                        input_to_filtered: dict, force: bool) -> (dict, list):
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
        # print("pair, pair_datasets")
        # print(pair, pair_datasets)

        data_dir = get_analysis_folder(i_datasets_folder, 'mmvec/common/data/%s' % pair)
        meta_dir = get_analysis_folder(i_datasets_folder, 'mmvec/common/metadata/%s' % pair)
        (omic1_, bool1), (omic2_, bool2) = pair_datasets
        omic1 = input_to_filtered[omic1_]
        omic2 = input_to_filtered[omic2_]
        # print()
        # print()
        # print('pair:', pair)
        # print(' >', omic1)
        # print(' >', omic2)
        if (omic1, bool1) not in filt_datasets or (omic2, bool2) not in filt_datasets:
            # print('[NOT IN]', omic1)
            # print('[NOT IN]', omic2)
            # print('[NOT IN]', filt_datasets.keys())
            # print('[NOT IN]', mmvec_pairs_filt)
            continue

        pair_filtering = filtering[pair]
        for preval_abund, preval_abund_dats in pair_filtering.items():
            preval_filt1, abund_filter1 = preval_abund_dats[(omic1, bool1)]
            preval_filt2, abund_filter2 = preval_abund_dats[(omic2, bool2)]
            filt1 = '%s_%s' % (preval_filt1, abund_filter1)
            filt2 = '%s_%s' % (preval_filt2, abund_filter2)
            tsv1, qza1, meta1, meta_pd1, sams1 = filt_datasets[(omic1, bool1)][preval_abund]
            tsv2, qza2, meta2, meta_pd2, sams2 = filt_datasets[(omic2, bool2)][preval_abund]
            # print()
            # print()
            # print('omic_filt1', omic_filt1)
            # print('omic_filt2', omic_filt2)
            # print('tsv1, qza1, meta1')
            # print(' -', tsv1)
            # print(' -', qza1)
            # print(' -', meta1)
            # print('tsv2, qza2, meta2')
            # print(' -', tsv2)
            # print(' -', qza2)
            # print(' -', meta2)
            common_sams = sorted(set(sams1) & set(sams2))
            len_common_sams = len(common_sams)
            if len_common_sams < 10:
                print('Not enough samples: %s (%s) vs %s (%s) -> skipping' % (omic1, filt1, omic2, filt2))
                continue

            meta_fp = '%s/meta_%s_%s_%s__%s_%s_%s__%s_%ss.tsv' % (
                meta_dir, omic1, preval_filt1, abund_filter1,
                omic2, preval_filt2, abund_filter2,
                pair, len_common_sams
            )
            new_tsv1 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
                data_dir, omic1, preval_filt1,
                abund_filter1, pair, len_common_sams
            )
            new_qza1 = '%s.qza' % splitext(new_tsv1)[0]
            new_tsv2 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
                data_dir, omic2, preval_filt2,
                abund_filter2, pair, len_common_sams
            )
            new_qza2 = '%s.qza' % splitext(new_tsv2)[0]

            common_datasets.setdefault(pair, []).append(
                [meta_fp, omic1, omic2, filt1, filt2,
                 new_tsv1, new_tsv2, new_qza1,
                 new_qza2, len_common_sams]
            )
            if meta_fp in common_datasets_done[pair]:
                print('\t\t\t* [DONE]', pair, ':', omic1, filt1, omic2, filt2)
                continue
            meta_subset1 = get_meta_common_sorted(meta_pd1, common_sams)
            meta_subset2 = get_meta_common_sorted(meta_pd2, common_sams)
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
            print(
                '\t\t\t* [TODO]', pair, ':',
                omic1, '[%s: %s]' % (filt1, meta_subset1.shape[0]),
                omic2, '[%s: %s]' % (filt2, meta_subset2.shape[0]))
    return common_datasets, common_jobs


def check_common_datasets(i_datasets_folder: str, mmvec_pairs: dict,
                          mmvec_filtering: dict, filt_datasets_pass: dict,
                          input_to_filtered: dict) -> (dict, list):
    """
    :param i_datasets_folder:
    :param mmvec_pairs:
    :param force: Force the re-writing of scripts for all commands.
    :return:
    """
    common_datasets_pass = {}
    for pair, pair_datasets in mmvec_pairs.items():
        pair_filtering = mmvec_filtering[pair]
        common_datasets_pass[pair] = []
        data_dir = get_analysis_folder(i_datasets_folder, 'mmvec/common/data/%s' % pair)
        meta_dir = get_analysis_folder(i_datasets_folder, 'mmvec/common/metadata/%s' % pair)
        (omic1_, bool1), (omic2_, bool2) = pair_datasets
        omic1 = input_to_filtered[omic1_]
        omic2 = input_to_filtered[omic2_]
        if (omic1, bool1) not in filt_datasets_pass or (omic2, bool2) not in filt_datasets_pass:
            continue
        for preval_abund in pair_filtering:
            preval_filt1, abund_filter1 = pair_filtering[preval_abund][(omic1, bool1)]
            preval_filt2, abund_filter2 = pair_filtering[preval_abund][(omic2, bool2)]
            if not filt_datasets_pass[(omic1, bool1)][preval_abund] or not filt_datasets_pass[(omic2, bool2)][preval_abund]:
                continue
            filt1 = '_'.join([preval_filt1, abund_filter1])
            filt2 = '_'.join([preval_filt2, abund_filter2])
            tsv1, qza1, meta1, meta_pd1, sams1 = filt_datasets_pass[(omic1, bool1)][preval_abund]
            tsv2, qza2, meta2, meta_pd2, sams2 = filt_datasets_pass[(omic2, bool2)][preval_abund]
            common_sams = sorted(set(sams1) & set(sams2))
            if len(common_sams) < 10:
                print('Not enough samples: %s (%s) vs %s (%s) -> skipping' % (omic1, filt1, omic2, filt2))
                continue
            meta_fp = '%s/meta_%s_%s_%s__%s_%s_%s__%s_%ss.tsv' % (
                meta_dir, omic1, preval_filt1, abund_filter1,
                omic2, preval_filt2, abund_filter2, pair, len(common_sams)
            )
            new_tsv1 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
                data_dir, omic1, preval_filt1, abund_filter1, pair, len(common_sams)
            )
            new_qza1 = '%s.qza' % splitext(new_tsv1)[0]
            new_tsv2 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
                data_dir, omic2, preval_filt2, abund_filter2, pair, len(common_sams)
            )
            new_qza2 = '%s.qza' % splitext(new_tsv2)[0]
            if isfile(meta_fp) and isfile(new_qza1) and isfile(new_qza2):
                common_datasets_pass[pair].append(meta_fp)
    return common_datasets_pass


def run_single_mmvec(odir: str, pair: str, meta_fp: str, qza1: str,
                     qza2: str,  res_dir: str, cur_sh: str, batch: str,
                     learn: str, epoch: str,  prior: str, thresh_feat: str,
                     latent_dim: str, train_column: str, n_example: str,
                     gpu: bool, force: bool, standalone: bool, qiime_env: str) -> None:
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
        # ranks_tsv = '%s/%s_ranks.tsv' % (pair, odir)
        # ordination_tsv = '%s/%s_ordination.txt' % (pair, odir)
        ranks_tsv = '%s/ranks.tsv' % odir
        ordination_tsv = '%s/ordination.txt' % odir
        if force or not isfile(ordination_tsv) or not isfile(ranks_tsv):
            write_mmvec_cmd(meta_fp, qza1, qza2, res_dir, odir,
                            ranks_tsv, ordination_tsv,
                            batch, learn, epoch, prior,
                            thresh_feat, latent_dim, train_column,
                            n_example, gpu, standalone, cur_sh_o, qiime_env)

            remove = False
    if remove:
        os.remove(cur_sh)


def check_filtered_and_common_dataset(
        i_datasets_folder:str, datasets: dict, datasets_filt: dict,
        unique_datasets: list, mmvec_pairs: dict, mmvec_filtering: dict,
        unique_filterings: dict, analysis: str, input_to_filtered: dict) -> (dict, dict):
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

    print('\t-> [%s] Check datasets filtered...' % analysis)
    filt_datasets_pass = check_datasets_filtered(
        i_datasets_folder, datasets, datasets_filt,
        unique_datasets, unique_filterings, analysis,
        input_to_filtered
    )
    filt_datasets_todo = [x for x, y in filt_datasets_pass.items() if not len(y)]
    if len(filt_datasets_todo):
        print('\t\t--> %s dataset(s) to prepare:' % len(filt_datasets_todo))
        for filt_datasets_td in filt_datasets_todo:
            print('\t\t\t*', filt_datasets_td)

    common_datasets_pass = {}
    if analysis == 'mmvec':
        print('\t-> [mmvec] Check common datasets...')
        common_datasets_pass = check_common_datasets(
            i_datasets_folder, mmvec_pairs, mmvec_filtering,
            filt_datasets_pass, input_to_filtered
        )
        common_datasets_todo = [x for x, y in common_datasets_pass.items() if not len(y)]
        if len(common_datasets_todo):
            print('\t\t--> %s common dataset(s) to prepare:' % len(common_datasets_todo))
            for common_datasets_td in common_datasets_todo:
                print('\t\t\t*', common_datasets_td)
    return filt_datasets_pass, common_datasets_pass


def make_filtered_and_common_dataset(
        i_datasets_folder:str, datasets: dict, datasets_filt: dict,
        datasets_read: dict, unique_datasets: list,
        mmvec_pairs: dict, filtering: dict, unique_filterings: dict,
        job_folder: str, force: bool,
        prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
        analysis: str, filt_raref: str, filt_datasets_done: dict,
        common_datasets_done: dict, input_to_filtered: dict,
        already_computed: dict) -> (dict, dict):
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

    print('\t-> [%s] Get datasets filtered...' % analysis)
    filt_datasets, filt_jobs = get_datasets_filtered(
        i_datasets_folder, datasets, datasets_read, datasets_filt,
        unique_datasets, unique_filterings, force, analysis,
        filt_datasets_done, input_to_filtered, already_computed)

    common_jobs = []
    common_datasets = {}
    if analysis == 'mmvec':
        print('\t-> [mmvec] Get common datasets...')
        common_datasets, common_jobs = get_common_datasets(
            i_datasets_folder, mmvec_pairs, filtering, filt_datasets,
            common_datasets_done, input_to_filtered, force
        )

    pre_jobs = filt_jobs + common_jobs
    if len(pre_jobs):
        import_sh = '%s/2_run_%s_imports%s.sh' % (job_folder, analysis, filt_raref)
        import_pbs = '%s.pbs' % splitext(import_sh)[0]
        with open(import_sh, 'w') as import_o:
            for cmd in pre_jobs:
                import_o.write('\necho "%s"\n' % cmd)
                import_o.write('%s\n' % cmd)
        run_xpbs(import_sh, import_pbs, '%s.mprt.mmsb.%s%s' % (prjct_nm, analysis, filt_raref),
                 qiime_env, '2', '1', '1', '150', 'mb', chmod, 1,
                 '# Import datasets for %s' % analysis, None, noloc)

    return filt_datasets, common_datasets


def get_unique_mmvec_filtering(mmvec_filtering):
    unique_filterings = {}
    for pair, filt_name_d in mmvec_filtering.items():
        for filt_name, dat_d in filt_name_d.items():
            for dat, prev_abund in dat_d.items():
                unique_filterings.setdefault(dat, []).append((filt_name, prev_abund))
    return unique_filterings


def run_mmvec(p_mmvec_pairs: str, i_datasets_folder: str, datasets: dict,
              datasets_filt: dict, datasets_read: dict, force: bool,
              gpu: bool, standalone: bool, prjct_nm: str, qiime_env: str,
              chmod: str, noloc: bool, split: bool, filt_raref: str,
              run_params: dict, input_to_filtered: dict) -> list:
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

    mmvec_dicts = get_mmvec_dicts(p_mmvec_pairs)
    mmvec_pairs = mmvec_dicts[0]
    mmvec_filtering = mmvec_dicts[1]
    mmvec_params = mmvec_dicts[2]
    mmvec_subsets = mmvec_dicts[3]
    unique_datasets = list(set([dat for pair_dats in mmvec_pairs.values() for dat in pair_dats]))
    unique_filterings = get_unique_mmvec_filtering(mmvec_filtering)


    print("mmvec_filtering")
    print(mmvec_filtering)

    print("unique_filterings.keys()")
    print(unique_filterings.keys())

    print("unique_datasets")
    print(unique_datasets)
    print(datasets_filtds)


    filt_datasets_done, common_datasets_done = check_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_filt,
        unique_datasets, mmvec_pairs, mmvec_filtering,
        unique_filterings, 'mmvec', input_to_filtered)

    already_computed = {}
    job_folder = get_job_folder(i_datasets_folder, 'mmvec')
    filt_datasets, common_datasets = make_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_filt, datasets_read,
        unique_datasets, mmvec_pairs, mmvec_filtering, unique_filterings, job_folder,
        force, prjct_nm, qiime_env, chmod, noloc, 'mmvec',
        filt_raref, filt_datasets_done, common_datasets_done,
        input_to_filtered, already_computed)

    jobs = []
    all_sh_pbs = {}
    mmvec_outputs = []

    for pair, pair_data in common_datasets.items():

        # print()
        # print("pair")
        # print(pair)
        # for idx, i in enumerate(pair_data):
        #     print(' -', idx)
        #     for j in i:
        #         print('     -', j)

        job_folder2 = get_job_folder(i_datasets_folder, 'mmvec/chunks/%s' % pair)
        if not split:
            out_sh = '%s/chunks/run_mmvec_%s%s.sh' % (job_folder, pair, filt_raref)

        for (meta_fp, omic1, omic2, filt1, filt2, tsv1, tsv2, qza1, qza2, ncommon) in pair_data:
            train_columns = mmvec_params['train_column']
            n_examples = mmvec_params['n_examples']
            batches = mmvec_params['batches']
            learns = mmvec_params['learns']
            epochs = mmvec_params['epochs']
            priors = mmvec_params['priors']
            thresh_feats = mmvec_params['thresh_feats']
            latent_dims = mmvec_params['latent_dims']
            if split:
                out_sh = '%s/chunks/run_mmvec_%s_%s_%s_%s_%s%s.sh' % (job_folder, pair, omic1, filt1,
                                                                      omic2, filt2, filt_raref)
            for idx, it in enumerate(itertools.product(train_columns, n_examples, batches, learns,
                                                       epochs, priors, thresh_feats, latent_dims)):
                train_column, n_example, batch, learn, epoch, prior, thresh_feat, latent_dim = [str(x) for x in it]
                res_dir = 'b-%s_l-%s_e-%s_p-%s_f-%s_d-%s_t-%s_n-%s_gpu-%s' % (
                    batch, learn, epoch, prior.replace('.', ''),
                    thresh_feat, latent_dim, train_column,
                    n_example, str(gpu)[0]
                )
                odir = get_analysis_folder(i_datasets_folder, 'mmvec/paired/%s/%s_%s__%s_%s/%s' % (
                    pair, omic1, filt1, omic2, filt2, res_dir
                ))
                mmvec_outputs.append([
                    pair, omic1, omic2, filt1, filt2,
                    ncommon, meta_fp, tsv1, tsv2, qza1, qza2,
                    'mmvec_out__%s' % res_dir, odir
                ])
                cur_sh = '%s/run_mmvec_%s_%s_%s_%s%s.sh' % (job_folder2, pair, filt1, filt2, res_dir, filt_raref)
                all_sh_pbs.setdefault((pair, out_sh), []).append(cur_sh)
                run_single_mmvec(
                    odir, pair, meta_fp,
                    qza1, qza2, res_dir, cur_sh,
                    batch, learn, epoch, prior, thresh_feat,
                    latent_dim, train_column, n_example,
                    gpu, force, standalone, qiime_env
                )
    if standalone:
        qiime_env = 'mmvec2'

    main_sh = write_main_sh(job_folder, '3_mmvec%s' % filt_raref, all_sh_pbs,
                            '%s.mmvc%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_mmvec_pairs.startswith('/panfs'):
            p_mmvec_pairs = p_mmvec_pairs.replace(os.getcwd(), '')
        print_message("# MMVEC (datasets pairs in %s)" % p_mmvec_pairs, 'sh', main_sh)

    return mmvec_outputs