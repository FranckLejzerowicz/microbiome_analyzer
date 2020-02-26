# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import numpy as np

from scipy.stats import skew
from os.path import basename, dirname, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics,
    get_job_folder,
    get_analysis_folder
)
from routine_qiime2_analyses._routine_q2_cmds import write_rarefy, run_export
np.set_printoptions(precision=2, suppress=True)


def run_rarefy(i_datasets_folder: str, datasets: dict, datasets_read: dict,
               datasets_features: dict, datasets_phylo: dict, datasets_raref_depths: dict,
               force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    Run rarefy: Rarefy table.
    https://docs.qiime2.org/2019.10/plugins/available/feature-table/rarefy/

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param datasets_raref_depths: The rarefaction depths.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: deta divesity matrices.
    """

    job_folder = get_job_folder(i_datasets_folder, 'rarefy')
    job_folder2 = get_job_folder(i_datasets_folder, 'rarefy/chunks')

    datasets_update = {}
    datasets_read_update = {}
    datasets_features_update = {}
    datasets_phylo_update = {}

    written = 0
    run_pbs = '%s/1_run_rarefy.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():

            if dat not in datasets_raref_depths:
                continue

            odir = get_analysis_folder(i_datasets_folder, 'rarefy/%s' % dat)
            depth = datasets_raref_depths[dat]
            dat_raref = '%s_raref%s' % (dat, depth)

            tsv, meta = tsv_meta_pds
            meta_out = '%s/meta_%s.tsv' % (odir, dat_raref)
            subprocess.call(['cp', meta, meta_out])

            out_sh = '%s/run_rarefy_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                qza = tsv.replace('.tsv', '.qza')
                qza_out = '%s/tab_%s.qza' % (odir, dat_raref)
                tsv_out = splitext(qza_out)[0]
                if force or not os.path.isfile(qza_out):
                    write_rarefy(qza, qza_out, depth, cur_sh)
                    written += 1
                if force or not os.path.isfile(tsv_out):
                    cmd = run_export(qza_out, tsv_out, 'FeatureTable[Frequency]')
                    cur_sh.write('echo "%s"\n' % cmd)
                    cur_sh.write('%s\n\n' % cmd)
                    written += 1

                datasets_update[dat_raref] = [tsv_out, meta_out]
                datasets_read_update[dat_raref] = [None, None]
                datasets_phylo_update[dat_raref] = datasets_phylo[dat]
                datasets_features_update[dat_raref] = None

            run_xpbs(out_sh, out_pbs, '%s.bt.%s' % (prjct_nm, dat),
                     qiime_env, '24', '1', '1', '10', 'gb',
                     chmod, written, 'single', o)
    if written:
        print('# Calculate beta diversity indices')
        print('[TO RUN] sh', run_pbs)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_features.update(datasets_features_update)
    datasets_phylo.update(datasets_phylo_update)


def check_rarefy_need(datasets_read: dict) -> dict:
    """
    Check the distribution of reads per sample and its skewness to
    warn user for the need for rarefaction of the feature tables.

    :param datasets_read: dataset -> [tsv table, meta table]
    :return datasets_raref_depths: Rarefaction depths for eac dataset.
    """
    datasets_raref_depths = {}
    for dat, (tsv_pd, meta_pd) in datasets_read.items():
        tsv_sam_sum = tsv_pd.sum()
        count, division = np.histogram(tsv_sam_sum)
        skw = skew(count)
        if abs(skw) > 1:
            print('[%s] Reads-per-sample distribution [skewness=%s] (>1!)' % (dat, skw))
            division_std = np.interp(count, (count.min(), count.max()), (0, 20))
            print('readsbin\tnsamples\thist_representation')
            for ddx, div in enumerate(division_std):
                if div > 1:
                    print('\t%s\t%s\t%s' % (format(division[ddx], '6.3E'), count[ddx], '-' * int(div)))
                elif div == 0:
                    print('\t%s\t%s\t%s' % (format(division[ddx], '6.3E'), count[ddx], ''))
                else:
                    print('\t%s\t%s\t%s' % (format(division[ddx], '6.3E'), count[ddx], '-'))
            print(' ==> Consider rarefying <==')
        second_quantile = tsv_sam_sum.quantile(0.2)
        if second_quantile < 1000:
            print('[%s] Second quantile of the reads-per-sample distribution is <1000' % dat)
            print('- The sequencing might have failed! Analyze with caution')
            print('- reads-per-sample distribution described:')
            for x,y in tsv_sam_sum.describe().to_dict():
                print('\t%s: %s' % (x, round(y, 3)))
        else:
            nfigure = len(str(int(second_quantile)))
            second_quantile_to_round = second_quantile / ( 10 ** (nfigure - 2) )
            second_quantile_rounded = round(second_quantile_to_round) * ( 10 ** (nfigure - 2) )
            print('[%s] Proposed rarefaction depth: %s (second quantile)' % (dat, second_quantile_rounded))
            datasets_raref_depths[dat] = str(second_quantile_rounded)
    return datasets_raref_depths