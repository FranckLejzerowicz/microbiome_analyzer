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

from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_songbird_dict,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict,
    check_metadata_models
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_songbird_cmd
)


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

    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(tensor):
            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
            new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
            write_songbird_cmd(qza, new_qza, new_meta,
                               formula, epoch, batch, diff_prior, learn, thresh_sample,
                               thresh_feat, n_random, diffs, diffs_qza, stats, plot,
                               base_diff_qza, base_stats, base_plot, tensor, cur_sh_o)
            remove = False
    if remove:
        os.remove(cur_sh)


def run_songbird(p_diff_models: str, i_datasets_folder: str, datasets: dict,
                 force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> dict:
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
    job_folder2 = get_job_folder(i_datasets_folder, 'songbird/chunks')
    main_models, main_cases_dict, params = get_songbird_dict(p_diff_models)

    batches = params['batches']
    learns = params['learns']
    epochs = params['epochs']
    thresh_feats = params['thresh_feats']
    thresh_samples = params['thresh_samples']
    diff_priors = params['diff_priors']

    songbird_outputs = {}

    jobs = []
    all_sh_pbs = {}
    for dat, tsv_meta_pds in datasets.items():

        tsv, meta = tsv_meta_pds
        qza = '%s.qza' % splitext(tsv)[0]
        meta_pd = read_meta_pd(meta)

        cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'songbird')
        models = check_metadata_models(meta, meta_pd, main_models, 'songbird')

        for model, formula in models.items():
            out_sh = '%s/run_songbird_%s_%s.sh' % (job_folder2, dat, model)
            for idx, it in enumerate(itertools.product(batches, learns, epochs, diff_priors,
                                                       thresh_feats, thresh_samples)):
                batch, learn, epoch, diff_prior, thresh_feat, thresh_sample = [str(x) for x in it]
                res_dir = '%s/%s/filt_f%s_s%s/%s_%s_%s_%s' % (
                    dat, model.replace('+', 'PLUS').replace('*', 'COMBI').replace('-', 'MINUS').replace('/', 'DIVIDE'),
                    thresh_feat, thresh_sample, batch, learn, epoch, diff_prior.replace('.', '')
                )
                songbird_outputs.setdefault(dat, []).append([dat, meta, qza, model, res_dir])
                odir = get_analysis_folder(i_datasets_folder, 'songbird/%s' % res_dir)
                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        case = get_case(case_vals, model, case_var, str(idx))
                        cur_sh = '%s/run_songbird_%s_%s.sh' % (
                            job_folder2, dat, case)
                        cur_sh = cur_sh.replace(' ', '-')
                        all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                        p = multiprocessing.Process(
                            target=run_multi_songbird,
                            args=(odir, qza, meta_pd, cur_sh, case, formula,
                                  batch, learn, epoch, diff_prior, thresh_feat, thresh_sample,
                                  case_var, case_vals, force))
                        p.start()
                        jobs.append(p)
    for j in jobs:
        j.join()

    job_folder = get_job_folder(i_datasets_folder, 'songbird')
    main_sh = write_main_sh(job_folder, '2_songbird', all_sh_pbs,
                            '%s.sngbrd' % prjct_nm, '2', '1', '24', '2', 'gb',
                            qiime_env, chmod)
    if main_sh:
        print("# Songbird (configs in %s)" % p_diff_models)
        print('[TO RUN] sh', main_sh)

    return songbird_outputs