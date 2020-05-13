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
from os.path import basename, isfile, splitext
import multiprocessing

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
    check_metadata_models
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_songbird_cmd
)
from routine_qiime2_analyses._routine_q2_mmbird import get_mmvec_outputs
from routine_qiime2_analyses._routine_q2_mmvec import make_filtered_and_common_dataset


def run_multi_songbird(odir: str, qza: str, meta_pd: pd.DataFrame, cur_sh: str,
                       case: str, formula: str, case_var: str, case_vals: list, force: bool,
                       batch: str, learn: str, epoch: str, diff_prior: str,
                       thresh_feat: str, thresh_sample: str, n_random: str) -> None:
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

    cur_rad = odir + '/' + basename(qza).replace('.qza', '_%s' % case)
    new_meta = '%s.meta' % cur_rad
    new_qza = '%s.qza' % cur_rad
    diffs = '%s/differentials.tsv' % cur_rad
    diffs_qza = '%s/differentials.qza' % cur_rad
    stats = '%s/differentials-stats.qza' % cur_rad
    plot = '%s/differentials-biplot.qza' % cur_rad
    base_diff_qza = '%s/differentials-baseline.qza' % cur_rad
    base_stats = '%s/differentials-baseline-stats.qza' % cur_rad
    base_plot = '%s/differentials-baseline-biplot.qza' % cur_rad
    tensor = '%s/differentials-tensorboard.qzv' % cur_rad
    tensor_html = '%s/differentials-tensorboard.html' % cur_rad

    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(tensor):
            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
            new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
            write_songbird_cmd(qza, new_qza, new_meta,
                               formula, epoch, batch, diff_prior, learn, thresh_sample,
                               thresh_feat, n_random, diffs, diffs_qza, stats, plot,
                               base_diff_qza, base_stats, base_plot, tensor, tensor_html,
                               cur_sh_o)
            remove = False
    if remove:
        os.remove(cur_sh)


def run_single_songbird(odir: str, qza: str, meta_pd: pd.DataFrame, cur_sh: str,
                        case: str, formula: str, case_var: str, case_vals: list, force: bool,
                        batch: str, learn: str, epoch: str, diff_prior: str,
                        thresh_feat: str, thresh_sample: str, n_random: str) -> (str, str):
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

    cur_rad = odir + '/' + basename(qza).replace('.qza', '_%s' % case)
    new_meta = '%s.meta' % cur_rad
    new_qza = '%s.qza' % cur_rad
    diffs = '%s_differentials.tsv' % cur_rad
    diffs_qza = '%s_differentials.qza' % cur_rad
    stats = '%s_differentials-stats.qza' % cur_rad
    plot = '%s_differentials-biplot.qza' % cur_rad
    base_diff_qza = '%s_differentials-baseline.qza' % cur_rad
    base_stats = '%s_differentials-baseline-stats.qza' % cur_rad
    base_plot = '%s_differentials-baseline-biplot.qza' % cur_rad
    tensor = '%s_differentials-tensorboard.qzv' % cur_rad
    tensor_html = '%s_differentials-tensorboard.html' % cur_rad

    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(tensor_html):
            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
            new_meta_pd.columns = [x.lower() for x in new_meta_pd.columns]
            new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
            write_songbird_cmd(qza, new_qza, new_meta,
                               formula, epoch, batch, diff_prior, learn, thresh_sample,
                               thresh_feat, n_random, diffs, diffs_qza, stats, plot,
                               base_diff_qza, base_stats, base_plot, tensor, tensor_html,
                               cur_sh_o)
            remove = False
    if remove:
        os.remove(cur_sh)
    return diffs, tensor_html


def run_songbird(p_diff_models: str, i_datasets_folder: str, datasets: dict,
                 datasets_read: dict, mmvec_outputs: list, force: bool,
                 prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> list:
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
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'songbird')
    job_folder2 = get_job_folder(i_datasets_folder, 'songbird/chunks')
    songbird_dicts = get_songbird_dicts(p_diff_models)
    songbird_models = songbird_dicts[0]
    songbird_filtering = songbird_dicts[1]
    params = songbird_dicts[2]
    songbird_datasets = songbird_dicts[3]
    main_cases_dict = songbird_dicts[4]

    # print("songbird_models")
    # print(songbird_models)
    # {
    #   'vioscreen_foods_consumed_grams_per_day_1800s_noLiquids':
    #       {'age': 'age_years', 'bmi': 'bmi'},
    #   'vioscreen_micromacro_qemistree_1800s':
    #       {'age': 'age_years', 'bmi': 'bmi'}}

    # print("songbird_filtering")
    # print(songbird_filtering)
    # {'prevalence': ['0', '10'],
    #  'abundance': [['0', '0'], ['1', '3']]}

    # print("main_cases_dict")
    # print(main_cases_dict)
    # {'ALL': [[]]}

    # print("params")
    # print(params)
    # {'batches': ['2'], 'learns': ['1e-3'], 'epochs': ['20'],
    #  'thresh_feats': ['0'], 'thresh_samples': ['0'], 'diff_priors': ['0.5'], 'n_randoms': ['250']}

    batches = params['batches']
    learns = params['learns']
    epochs = params['epochs']
    thresh_feats = params['thresh_feats']
    thresh_samples = params['thresh_samples']
    diff_priors = params['diff_priors']
    n_randoms = params['n_randoms']

    filt_datasets, common_datasets = make_filtered_and_common_dataset(
        i_datasets_folder, datasets, datasets_read,
        songbird_datasets, {}, songbird_filtering,
        job_folder, force, prjct_nm, qiime_env,
        chmod, noloc, 'songbird')

    songbirds = {}
    for dat, filts_files in filt_datasets.items():
        for filts, files in filts_files.items():
            songbirds.setdefault(dat, []).append(['_'.join(filts), files[0], files[2]])

    if mmvec_outputs:
        mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
        for r, row in mmvec_outputs_pd.iterrows():
            dat1 = row['omic1']
            dat2 = row['omic2']
            dat_filt1 = row['omic_filt1']
            dat_filt2 = row['omic_filt2']
            omic1_common_fp = row['omic1_common_fp']
            omic2_common_fp = row['omic2_common_fp']
            meta_common_fp = row['meta_common_fp']
            songbirds.setdefault(dat1, []).append([dat_filt1, omic1_common_fp, meta_common_fp])
            songbirds.setdefault(dat2, []).append([dat_filt2, omic2_common_fp, meta_common_fp])

    jobs = []
    all_sh_pbs = {}
    first_print = 0
    songbird_outputs = []
    for dat, filts_tsvs_metas in songbirds.items():

        out_sh = '%s/run_songbird_%s.sh' % (job_folder2, dat)
        for (filt, tsv, meta_) in filts_tsvs_metas:

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
            meta_pd = meta_pd.set_index('sample_name')

            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'songbird')
            if dat in songbird_models:
                print()
                print(dat)
                print(meta)
                models = check_metadata_models(meta, meta_pd, songbird_models[dat], 'songbird')
            else:
                continue

            for model, formula in models.items():
                for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors,
                                                           thresh_feats, thresh_samples, n_randoms)):
                    batch, learn, epoch, diff_prior, thresh_feat, thresh_sample, n_random = [str(x) for x in it]
                    model_rep = model.replace('+', 'PLUS').replace('*', 'COMBI').replace('-', 'MINUS').replace('/', 'DIVIDE')
                    params = '%s/%s/filt_f%s_s%s/%s_%s_%s_%s' % (
                        filt, model_rep, thresh_feat, thresh_sample,
                        batch, learn, epoch, diff_prior.replace('.', '')
                    )
                    res_dir = '%s/%s' % (dat, params)
                    odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % res_dir)
                    for case_var, case_vals_list in cases_dict.items():
                        for case_vals in case_vals_list:
                            case = get_case(case_vals, case_var, str(idx))
                            cur_sh = '%s/run_songbird_%s_%s_%s_%s.sh' % (
                                job_folder2, dat, filt, model_rep, case)
                            cur_sh = cur_sh.replace(' ', '-')
                            all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                            diffs, tensor_html = run_single_songbird(
                                odir, qza, meta_pd, cur_sh, case, formula,
                                case_var, case_vals, force, batch, learn,
                                epoch, diff_prior, thresh_feat,
                                thresh_sample, n_random
                            )
                            if '/mmvec/' in tsv:
                                songbird_outputs.append([dat, filt, params.replace('/', '__'),
                                                         case_var, case, diffs, tensor_html])
                            # p = multiprocessing.Process(
                            #     target=run_multi_songbird,
                            #     args=(odir, qza, meta_pd, cur_sh, case, formula, case_var, case_vals, force,
                            #           batch, learn, epoch, diff_prior, thresh_feat, thresh_sample, n_random
                            #           ))
                            # p.start()
                            # jobs.append(p)
    # for j in jobs:
    #     j.join()

    job_folder = get_job_folder(i_datasets_folder, 'songbird')
    main_sh = write_main_sh(job_folder, '2_songbird', all_sh_pbs,
                            '%s.sngbrd' % prjct_nm, '150', '1', '1', '25', 'gb',
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_diff_models.startswith('/panfs'):
            p_diff_models = p_diff_models.replace(os.getcwd(), '')
        print_message("# Songbird (configs in %s)" % p_diff_models, 'sh', main_sh)
    return songbird_outputs