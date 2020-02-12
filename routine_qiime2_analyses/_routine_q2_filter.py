# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
from os.path import splitext

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, run_import


def import_datasets(i_folder: str, datasets: dict, force: bool,
                    prjct_nm: str, qiime_env: str) -> None:

    print('# Import tables to qiime2')
    job_folder = get_job_folder(i_folder, 'import_tables')

    out_sh = '%s/0_run_import.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    written = 0
    with open(out_sh, 'w') as sh:
        for dataset, tsvs_metas_lists in datasets.items():
            for tsv_meta in tsvs_metas_lists:
                cmd = run_import(tsv_meta[0], '%s.qza' % splitext(tsv_meta[0])[0], 'FeatureTable[Frequency]')
                sh.write('%s\n' % cmd)
                written += 1

    if force or written:
        xpbs_call(out_sh, out_pbs, '%s.mprt' % prjct_nm, qiime_env,
                  '1', '1', '1', '100', 'mb')
        print('[TO RUN] qsub', out_pbs)
    else:
        print('\nNothing written in', out_sh, '--> removed')
        os.remove(out_sh)


def filter_rare_samples(i_folder: str, datasets: dict, datasets_read: dict, datasets_features: dict,
                        force: bool, prjct_nm: str, qiime_env: str, thresh: int, gid: bool) -> (dict, dict, dict):
    """
    Filter the rare features, keep samples with enough reads/features and import to Qiime2.

    :param i_folder: Path to the folder containing the .tsv datasets.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param thresh: min number of reads per sample to keep it.
    :param gid: If feature names have the genome ID (to use the Web of Life tree).
    :return:
    """

    print('# Filter samples for a min number of %s reads' % thresh)
    job_folder = get_job_folder(i_folder, 'import_filtered')

    written = 0
    datasets_update = {}
    datasets_read_update = {}
    datasets_features_update = {}
    out_sh = '%s/1_run_import_filtered.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    with open(out_sh, 'w') as sh:
        for dat, tabs_metas_pds in datasets_read.items():
            for tab_pd, meta_pd in tabs_metas_pds:

                meta_pd = meta_pd.set_index(meta_pd.columns.tolist()[0])
                tab_filt_pd = tab_pd.loc[:, tab_pd.sum(0) >= thresh].copy()
                tab_filt_pd = tab_filt_pd.loc[tab_filt_pd.sum(1) > 0, :]
                tab_filt_pd = tab_filt_pd.loc[:, tab_filt_pd.sum(0) > 0]
                meta_filt_pd = meta_pd.loc[tab_filt_pd.columns.tolist()].copy()
                dat_filt = '%s_min%s_%ss' % (dat, thresh, meta_filt_pd.shape[0])
                tab_filt_fp = '%s/data/tab_%s.tsv' % (i_folder, dat_filt)
                meta_filt_fp = tab_filt_fp.replace('/data/', '/metadata/').replace('tab_', 'meta_')
                datasets_update.setdefault(dat_filt, []).append([tab_filt_fp, meta_filt_fp])

                tab_filt_pd.reset_index().to_csv(tab_filt_fp, index=False, sep='\t')
                meta_filt_pd.reset_index().to_csv(meta_filt_fp, index=False, sep='\t')
                datasets_read_update.setdefault(dat_filt, []).append([tab_filt_pd, meta_filt_pd.reset_index()])

                if gid and str(tab_filt_pd.index.dtype) == 'object':
                    features_names = tab_filt_pd.index.tolist()
                    found_gids = {}
                    for features_name in features_names:
                        if re.search('G\d{9}', features_name):
                            found_gids[re.search('G\d{9}', features_name).group(0)] = features_name
                    if len(found_gids) == len(features_names):
                        datasets_features[dat_filt] = found_gids

                qza = tab_filt_fp.replace('.tsv', '.qza')

                cmd = run_import(tab_filt_fp, qza, "FeatureTable[Frequency]")
                sh.write('%s\n' % ' '.join(cmd))
                written += 1

    if force or written:
        xpbs_call(out_sh, out_pbs, '%s.fltr' % prjct_nm, qiime_env,
                  '4', '4', '1', '100', 'mb')
        print('[TO RUN] qsub', out_pbs)
    else:
        print('\nNothing written in', out_sh, '--> removed')
        os.remove(out_sh)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_features.update(datasets_features_update)
