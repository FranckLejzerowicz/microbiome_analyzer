# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import skew
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    read_yaml_file, get_job_folder, get_analysis_folder, simple_chunks)
from routine_qiime2_analyses._routine_q2_cmds import write_rarefy, run_export
np.set_printoptions(precision=2, suppress=True)


def get_digit_depth(
        depth_: str,
        tsv_sums: pd.Series) -> int:
    """Get the rarefaction depth integer

    Parameters
    ----------
    depth_ : str
        Rarefaction depth
    tsv_sums : pd.Series
        Sum of reads per sample

    Returns
    -------
    depth : int
        Rarefaction depth
    """
    if depth_.isdigit():
        depth = int(depth_)
    else:
        depth = int(np.floor(min(tsv_sums)))
    return depth


def run_rarefy(
        i_datasets_folder: str, datasets: dict, datasets_read: dict,
        datasets_phylo: dict, datasets_filt_map: dict, datasets_rarefs: dict,
        p_raref_depths: str, eval_rarefs: bool, force: bool, prjct_nm: str,
        qiime_env: str, chmod: str, noloc: bool, run_params: dict,
        filt_raref: str, filt_only: bool, jobs: bool, slurm: bool,
        chunkit: int) -> dict:
    """
    Run rarefy: Rarefy table.
    https://docs.qiime2.org/2019.10/plugins/available/feature-table/rarefy/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: deta divesity matrices.
    """

    evaluation = ''
    eval_depths = {}
    datasets_raref_depths, datasets_raref_evals = check_rarefy_need(
        i_datasets_folder, datasets_read, p_raref_depths)
    if eval_rarefs:
        evaluation = '_eval'

    set_filt_rarefy(datasets_raref_depths, datasets_filt_map)

    datasets_update = {}
    datasets_read_update = {}
    datasets_phylo_update = {}
    datasets_append = {}

    main_written = 0
    job_folder = get_job_folder(
        i_datasets_folder, 'rarefy%s' % evaluation)
    job_folder2 = get_job_folder(
        i_datasets_folder, 'rarefy%s/chunks' % evaluation)
    to_chunk = []
    run_pbs = '%s/1_run_rarefy_%s%s%s.sh' % (
        job_folder, prjct_nm, evaluation, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds_ in datasets.items():

            written = 0
            if dat not in datasets_raref_depths:
                continue
            if filt_only and dat not in datasets_filt_map:
                continue

            odir = get_analysis_folder(
                i_datasets_folder, 'rarefy%s/%s' % (evaluation, dat))
            out_sh = '%s/run_rarefy_%s%s_%s.sh' % (
                job_folder2, prjct_nm, evaluation, dat)
            if slurm:
                out_pbs = '%s.slm' % splitext(out_sh)[0]
            else:
                out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:

                depths = datasets_raref_depths[dat][1]
                if eval_rarefs:
                    depths = datasets_raref_evals[dat]

                tsv_pd, meta_pd = datasets_read[dat][0]
                tsv_sums = tsv_pd.sum()
                for tsv_meta_pds in tsv_meta_pds_:
                    tsv, meta = tsv_meta_pds
                    for depth_ in depths:
                        depth = get_digit_depth(depth_, tsv_sums)
                        dat_raref = '%s_raref%s%s' % (
                            dat, evaluation, str(depth))
                        meta_out = '%s/meta_%s.tsv' % (odir, dat_raref)
                        remaining_samples = tsv_sums[
                            tsv_sums >= depth].index.tolist()
                        meta_raref_pd = meta_pd.loc[
                            meta_pd.sample_name.isin(remaining_samples), :]
                        meta_raref_pd.to_csv(meta_out, index=False, sep='\t')

                        qza = tsv.replace('.tsv', '.qza')
                        qza_out = '%s/tab_%s.qza' % (odir, dat_raref)
                        tsv_out = '%s.tsv' % splitext(qza_out)[0]
                        if force or not os.path.isfile(tsv_out):
                            cmd = write_rarefy(qza, qza_out, depth)
                            cur_sh.write('echo "%s"\n' % cmd)
                            cur_sh.write('%s\n\n' % cmd)
                            cmd = run_export(
                                qza_out, tsv_out, 'FeatureTable[Frequency]')
                            cur_sh.write('echo "%s"\n' % cmd)
                            cur_sh.write('%s\n\n' % cmd)
                            main_written += 1
                            written += 1

                        if eval_rarefs:
                            eval_depths.setdefault(dat, []).append(
                                '%s_%s' % (dat, str(depth)))
                            datasets_update['%s_%s' % (
                                dat, str(depth))] = [[tsv_out, meta_out]]
                            datasets_read_update['%s_%s' % (
                                dat, str(depth))] = ('raref', str(depth))
                            datasets_phylo_update['%s_%s' % (
                                dat, str(depth))] = datasets_phylo[dat]
                        else:
                            datasets_append.setdefault(dat, []).append(
                                [tsv_out, meta_out])
                            if isfile(tsv_out) and isfile(meta_out):
                                tab_filt_pd = pd.read_csv(
                                    tsv_out, index_col=0, header=0, sep='\t')
                                with open(meta_out) as f:
                                    for line in f:
                                        break
                                meta_filt_pd = pd.read_csv(
                                    meta_out, header=0, sep='\t',
                                    dtype={line.split('\t')[0]: str},
                                    low_memory=False)
                                datasets_read[dat].append([tab_filt_pd,
                                                           meta_filt_pd])
                            else:
                                datasets_read[dat].append(('raref', str(depth)))
                            datasets_rarefs.setdefault(dat, []).append(
                                '_raref%s%s' % (evaluation, str(depth)))

            to_chunk.append(out_sh)
            if not chunkit:
                run_xpbs(
                    out_sh, out_pbs, '%s.bt%s.%s%s' % (prjct_nm, evaluation,
                                                       dat, filt_raref),
                    qiime_env, run_params["time"], run_params["n_nodes"],
                    run_params["n_procs"], run_params["mem_num"],
                    run_params["mem_dim"], chmod, written, 'single', o,
                    noloc, slurm, jobs)

    if to_chunk and chunkit:
        simple_chunks(
            run_pbs, job_folder2, to_chunk, 'rarefy%s' % evaluation,
            prjct_nm, run_params["time"], run_params["n_nodes"],
            run_params["n_procs"], run_params["mem_num"], run_params["mem_dim"],
            qiime_env, chmod, noloc, slurm, jobs, chunkit, None)

    if main_written:
        print_message('# Get rarefied datasets', 'sh', run_pbs, jobs)

    if eval_rarefs:
        datasets.update(datasets_update)
        datasets_read.update(datasets_read_update)
        datasets_phylo.update(datasets_phylo_update)
    else:
        for dat, fps in datasets_append.items():
            datasets[dat].extend(fps)

    return eval_depths


def get_datasets_raref_evals(
        sam_sum: pd.Series) -> set:
    """

    Parameters
    ----------
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    datasets_raref_evals : set
    """
    datasets_raref_evals = set([
        int(x) for x in sam_sum.describe(
            percentiles=[
                x / 100 for x in range(10, 101, 10)
            ])[4:-1]
    ])
    return datasets_raref_evals


def get_depths(
        depths_yml: list,
        dat: str,
        sam_sum: pd.Series) -> tuple:
    """

    Parameters
    ----------
    depths_yml : list
        Depths at which to rarefy
    dat : str
        Dataset name
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    skip : bool
        Whether to skip a rarefaction
    depths : list
        Depths at which to rarefy
    """
    skip = False
    depths = []
    for depth in depths_yml:
        if depth == 'min' or sum(sam_sum >= int(depth)) > 10:
            depths.append(depth)
    if not depths:
        print('[%s] Min. proposed rarefaction depths would leave <10 samples: '
              '%s (not rarefaction)' % (dat, ', '.join(depths_yml)))
        skip = True
    elif len(depths) != len(depths_yml):
        print('[%s] Proposed rarefaction depths would leave <10 samples: %s ('
              'not rarefied)' % (dat, ', '.join([x for x in depths_yml
                                                 if x not in depths])))
    return skip, depths


def print_not_rarefying(
        dat: str,
        tsv_sam_sum: pd.Series) -> None:
    """
    Parameters
    ----------
    dat: str
        Dataset name
    tsv_sam_sum : pd.Series
        Sum of reads per sample
    """
    print('[%s] Second quantile of the reads-per-sample '
          'distribution is <1000' % dat)
    print('- The sequencing might have failed! Analyze with caution')
    print('- reads-per-sample distribution described:')
    for x, y in tsv_sam_sum.describe().to_dict().items():
        print('\t%s: %s' % (x, round(y, 3)))
    print('!!! NOT RAREFYING %s !!!' % dat)


def print_skew(
        dat: str,
        tsv_sam_sum: pd.Series) -> None:
    """
    Parameters
    ----------
    dat : str
        Dataset name
    tsv_sam_sum : pd.Series
        Sum of reads per sample
    """
    count, division = np.histogram(tsv_sam_sum)
    skw = skew(count)
    if abs(skw) > 1:
        print()
        print(' ==> Consider rarefying <==')
        print('[%s] Reads-per-sample distribution [skewness=%s] (>1!)' % (
            dat, round(abs(float(skw)), 3)))
        division_std = np.interp(
            count, (min(count), max(count)), (0, 20))
        print('\treadsbin\tsamples\thistogram')
        for ddx, div in enumerate(division_std):
            if div > 1:
                print('\t%s\t%s\t%s' % (
                format(division[ddx], '6.3E'), count[ddx], '-' * int(div)))
            elif div == 0:
                print('\t%s\t%s\t%s' % (
                format(division[ddx], '6.3E'), count[ddx], ''))
            else:
                print('\t%s\t%s\t%s' % (
                format(division[ddx], '6.3E'), count[ddx], '-'))


def get_default_raref_depth(
        i_datasets_folder: str,
        dat: str,
        sam_sum: pd.Series):
    """
    Parameters
    ----------
    i_datasets_folder : str
        Path to the folder containing the data/metadata sub-folders
    dat : str
        Dataset name
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    depths : list
        Rarefaction depths
    """
    raref_files = glob.glob('%s/qiime/rarefy/%s/tab_raref*.qza' % (
        i_datasets_folder, dat))
    if len(raref_files):
        depths = [
            x.split('_raref')[-1].split('.tsv')[0] for x in raref_files
        ]
    else:
        second_quantile = sam_sum.quantile(0.2)
        print_skew(dat, sam_sum)
        if second_quantile < 1000:
            depths = []
            print_not_rarefying(dat, sam_sum)
        else:
            nfigure = len(str(int(second_quantile)))
            second_quantile_to_round = second_quantile / (10 ** (nfigure - 2))
            second_quantile_rounded = round(second_quantile_to_round) * (
                        10 ** (nfigure - 2))
            depths = [str(int(second_quantile_rounded))]
            print('[%s] Proposed rarefaction depth: %s '
                  '(second quantile)' % (dat, depths[0]))
    return depths


def get_dat_depths(
        dat: str,
        i_datasets_folder: str,
        depths_yml: dict,
        sam_sum: pd.Series) -> tuple:
    """
    Parameters
    ----------
    dat : str
        Dataset name
    i_datasets_folder : str
        Path to the folder containing the data/metadata sub-folders
    depths_yml : dist
        Mapping Dataset nanme -> Depths at which to rarefy
    sam_sum : pd.Series
        Sum of reads per sample

    Returns
    -------
    skip : bool
        Whether to skip a rarefaction
    depths_tuple : tuple
        (boolean, Rarefaction depths)
    """
    skip = False
    if not depths_yml:
        depths = get_default_raref_depth(i_datasets_folder, dat, sam_sum)
        depths_tuple = (0, depths)
    elif dat in depths_yml:
        skip, depths = get_depths(depths_yml[dat], dat, sam_sum)
        depths_tuple = (1, depths)
    else:
        skip = True
        depths_tuple = []
    return skip, depths_tuple


def check_rarefy_need(i_datasets_folder: str,
                      datasets_read: dict,
                      p_raref_depths: str) -> (dict, dict):
    """Check the distribution of reads per sample and its skewness to
    warn user for the need for rarefaction of the feature tables.

    Parameters
    ----------
    i_datasets_folder : str
        Path to the folder containing the data/metadata sub-folders
    datasets_read : dict
    p_raref_depths : str
        Path to a rarefaction config yaml file.

    Returns
    -------

    """
    depths_yml = read_yaml_file(p_raref_depths)
    datasets_raref_depths = {}
    datasets_raref_evals = {}
    for dat, tsv_meta_pds in datasets_read.items():
        for (tsv_pd, meta_pd) in tsv_meta_pds:
            tsv_sam_sum = tsv_pd.sum()
            # not in the class code yet
            datasets_raref_evals[dat] = get_datasets_raref_evals(tsv_sam_sum)
            # if depths_yml:
            #     if dat in depths_yml:
            #     continue
            skip, depths = get_dat_depths(
                dat, i_datasets_folder, depths_yml, tsv_sam_sum)
            if skip:
                continue
            datasets_raref_depths[dat] = depths
            if depths[0]:
                # not in the class code yet
                datasets_raref_evals[dat].update([
                    int(x) if str(x).isdigit()
                    else np.floor(min(tsv_sam_sum))
                    for x in depths[1]
                ])
    return datasets_raref_depths, datasets_raref_evals


def set_filt_rarefy(
        datasets_raref_depths: dict,
        datasets_filt_map: dict) -> None:
    """
    Parameters
    ----------
    datasets_raref_depths : dict
    datasets_filt_map : dict
    """
    for dat_filt, dat in datasets_filt_map.items():
        if datasets_raref_depths[dat][0]:
            datasets_raref_depths[dat_filt] = datasets_raref_depths[dat]

