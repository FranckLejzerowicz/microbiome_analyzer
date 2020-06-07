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
from routine_qiime2_analyses._routine_q2_mmvec import make_filtered_and_common_dataset


def run_single_songbird(odir: str, qza: str, meta_pd: pd.DataFrame, cur_sh: str,
                        case: str, formula_meta_var_drop: list, case_var: str, case_vals: list, force: bool,
                        batch: str, learn: str, epoch: str, diff_prior: str,
                        thresh_feat: str, thresh_sample: str, n_random: str,
                        baselines: dict, baseline: str) -> (str, str):
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
    :param n_random:
    :param force: Force the re-writing of scripts for all commands.
    """
    remove = True

    new_meta = '%s/metadata.tsv' % odir
    new_qza = '%s/tab.qza' % odir
    diffs = '%s/differentials.tsv' % odir
    diffs_qza = '%s/differentials.qza' % odir
    stats = '%s/differentials-stats.qza' % odir
    plot = '%s/differentials-biplot.qza' % odir
    if baseline in baselines:
        base_diff_qza = ''
        base_stats = baselines[baseline]
        base_plot = ''
    else:
        base_diff_qza = '%s/differentials-baseline.qza' % odir
        base_stats = '%s/differentials-baseline-stats.qza' % odir
        base_plot = '%s/differentials-baseline-biplot.qza' % odir
        baselines[baseline] = base_stats
    tensor = '%s/differentials-tensorboard.qzv' % odir
    tensor_html = '%s/differentials-tensorboard.html' % odir

    formula, vars, meta_var, drop = formula_meta_var_drop
    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(tensor_html):
            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
            new_meta_pd.columns = [x.lower() for x in new_meta_pd.columns]
            meta_vars = list(set(list(vars) + list(meta_var)))
            new_meta_pd = new_meta_pd[meta_vars]
            new_meta_pd = rename_duplicate_columns(new_meta_pd)
            if len(drop):
                new_meta_pd = new_meta_pd.loc[(~new_meta_pd[meta_var.lower()].isin(drop)), :]
            new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
            write_songbird_cmd(qza, new_qza, new_meta, formula, epoch,
                               batch, diff_prior, learn, thresh_sample,
                               thresh_feat, n_random, diffs, diffs_qza, stats, plot,
                               base_diff_qza, base_stats, base_plot, tensor, tensor_html,
                               cur_sh_o)
            remove = False
    if remove:
        os.remove(cur_sh)
    return diffs, tensor_html, baselines


def run_songbird(p_diff_models: str, i_datasets_folder: str, datasets: dict,
                 datasets_read: dict, datasets_filt_map: dict, mmvec_outputs: list, force: bool,
                 prjct_nm: str, qiime_env: str, chmod: str, noloc: bool,
                 split: bool, filt_raref: str) -> list:
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
    params = songbird_dicts[2]
    songbird_datasets = songbird_dicts[3]
    main_cases_dict = songbird_dicts[4]

    batches = params['batches']
    learns = params['learns']
    epochs = params['epochs']
    thresh_feats = params['thresh_feats']
    thresh_samples = params['thresh_samples']
    diff_priors = params['diff_priors']
    n_randoms = params['n_randoms']

    filt_datasets, common_datasets = make_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_filt_map,
        datasets_read, songbird_datasets, {}, songbird_filtering,
        job_folder, force, prjct_nm, qiime_env,
        chmod, noloc, 'songbird', filt_raref)

    songbirds = {}
    for dat, filts_files in filt_datasets.items():
        for filts, files in filts_files.items():
            songbirds.setdefault(dat, []).append(['_'.join(filts), files[0], files[2], ''])

    if mmvec_outputs:
        mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
        for r, row in mmvec_outputs_pd.iterrows():
            pair = row['pair']
            omic1 = row['omic1']
            omic2 = row['omic2']
            filt1 = row['filt1']
            filt2 = row['filt2']
            omic1_common_fp = row['omic1_common_fp']
            omic2_common_fp = row['omic2_common_fp']
            meta_common_fp = row['meta_common_fp']
            songbirds.setdefault(omic1, []).append([filt1, omic1_common_fp, meta_common_fp, pair])
            songbirds.setdefault(omic2, []).append([filt2, omic2_common_fp, meta_common_fp, pair])

    all_sh_pbs = {}
    first_print = 0
    songbird_outputs = []
    for dat, filts_tsvs_metas_pair in songbirds.items():

        if not split:
            out_sh = '%s/run_songbird_%s%s.sh' % (job_folder2, dat, filt_raref)
        for (filt, tsv, meta_, pair) in filts_tsvs_metas_pair:

            if split:
                out_sh = '%s/run_songbird_%s_%s_%s%s.sh' % (job_folder2, dat, filt, pair, filt_raref)

            meta_alphas = '%s_alphas.tsv' % splitext(meta_)[0]
            if isfile(meta_alphas):
                meta = meta_alphas
            else:
                meta = meta_
                if not first_print:
                    print('\nWarning: Make sure you first run alpha -> alpha merge -> alpha export\n'
                          '\t(if you have alpha diversity as a factors in the models)!')
                    first_print += 1

            qza = '%s.qza' % splitext(tsv)[0]
            meta_pd = read_meta_pd(meta)
            meta_pd = rename_duplicate_columns(meta_pd)
            meta_pd = meta_pd.set_index('sample_name')
            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'songbird')
            if dat in songbird_models:
                models = check_metadata_models(meta, meta_pd, songbird_models[dat], 'songbird')
            else:
                continue

            baselines = {}
            for model, formula_meta_var_drop in models.items():
                for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors,
                                                           thresh_feats, thresh_samples, n_randoms)):
                    batch, learn, epoch, diff_prior, thresh_feat, thresh_sample, n_random = [str(x) for x in it]
                    model_rep = model.replace('+', 'PLUS').replace('*', 'COMBI').replace('-', 'MINUS').replace('/', 'DIVIDE')

                    params = '%s/filt_f%s_s%s/%s_%s_%s_%s' % (
                        model_rep, thresh_feat, thresh_sample,
                        batch, learn, epoch, diff_prior.replace('.', '')
                    )
                    for case_var, case_vals_list in cases_dict.items():
                        for case_vals in case_vals_list:
                            case = get_case(case_vals, case_var, str(idx))
                            baseline = '%s_%s' % ('_'.join([str(x) for x in it]), case)
                            if pair:
                                res_dir = '%s/%s/%s/%s/%s' % (dat, pair, filt, case, params)
                            else:
                                res_dir = '%s/%s/%s/%s' % (dat, filt, case, params)
                            odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % res_dir)
                            cur_sh = '%s/run_songbird_%s_%s_%s_%s_%s%s.sh' % (
                                job_folder2, dat, filt, model_rep, case, pair, filt_raref)
                            cur_sh = cur_sh.replace(' ', '-')
                            all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                            diffs, tensor_html = run_single_songbird(
                                odir, qza, meta_pd, cur_sh, case, formula_meta_var_drop,
                                case_var, case_vals, force, batch, learn,
                                epoch, diff_prior, thresh_feat,
                                thresh_sample, n_random, baselines, baseline
                            )
                            if pair:
                                songbird_outputs.append([dat, filt, params.replace('/', '__'),
                                                         case, diffs, pair])

    job_folder = get_job_folder(i_datasets_folder, 'songbird')
    main_sh = write_main_sh(job_folder, '2_songbird%s' % filt_raref, all_sh_pbs,
                            '%s.sngbrd%s' % (prjct_nm, filt_raref), '48', '1', '1', '2', 'gb',
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_diff_models.startswith('/panfs'):
            p_diff_models = p_diff_models.replace(os.getcwd(), '')
        print_message("# Songbird (configs in %s)" % p_diff_models, 'sh', main_sh)
    return songbird_outputs
