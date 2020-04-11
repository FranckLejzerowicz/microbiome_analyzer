# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_raref_tab_meta_pds
from routine_qiime2_analyses._routine_q2_cmds import run_import


def import_datasets(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
                    force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    Initial import of the .tsv datasets in to Qiime2 Artefact.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'import_tables')

    out_sh = '%s/0_run_import.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    written = 0
    with open(out_sh, 'w') as sh:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            qza = '%s.qza' % splitext(tsv)[0]
            if datasets_phylo[dat][1]:
                cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n' % cmd)
                written += 1
            elif force or not isfile(qza):
                cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n' % cmd)
                written += 1
    run_xpbs(out_sh, out_pbs, '%s.mprt' % prjct_nm,
             qiime_env,'1', '1', '1', '100', 'mb', chmod, written,
             '# Import tables to qiime2')


def filter_rare_samples(i_datasets_folder: str, datasets: dict, datasets_read: dict,
                        datasets_features: dict, datasets_phylo: dict,
                        prjct_nm: str, qiime_env: str, thresh: int, chmod: str) -> None:
    """
    Filter the rare features, keep samples with enough reads/features and import to Qiime2.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param thresh: min number of reads per sample to keep it.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    written = 0
    datasets_update = {}
    datasets_read_update = {}
    datasets_features_update = {}
    datasets_phylo_update = {}
    job_folder = get_job_folder(i_datasets_folder, 'import_filtered')
    out_sh = '%s/1_run_import_filtered.sh' % job_folder
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    with open(out_sh, 'w') as sh:
        for dat, tab_meta_pds in datasets_read.items():

            if tab_meta_pds == 'raref':
                tsv, meta = datasets[dat]
                if not isfile(tsv):
                    print('Must have run rarefcation to use it further...\nExiting')
                    sys.exit(0)
                tab_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                datasets_read[dat] = [tab_pd, meta_pd]
            else:
                tab_pd, meta_pd = tab_meta_pds

            if datasets_features[dat] == 'raref':
                dat_no_raref = '_raref'.join(dat.split('_raref')[:-1])
                datasets_features[dat] = dict(
                    gid_feat for gid_feat in datasets_features[dat_no_raref].items() if gid_feat[1] in tab_pd.index
                )

            meta_pd = meta_pd.set_index(meta_pd.columns.tolist()[0])
            tab_filt_pd = tab_pd.loc[:, tab_pd.sum(0) >= thresh].copy()
            tab_filt_pd = tab_filt_pd.loc[tab_filt_pd.sum(1) > 0, :]
            tab_filt_pd = tab_filt_pd.loc[:, tab_filt_pd.sum(0) > 0]
            meta_filt_pd = meta_pd.loc[tab_filt_pd.columns.tolist()].copy()
            dat_filt = '%s_min%s_%ss' % (dat, thresh, meta_filt_pd.shape[0])
            tab_filt_fp = '%s/data/tab_%s.tsv' % (i_datasets_folder, dat_filt)
            meta_filt_fp = tab_filt_fp.replace('/data/', '/metadata/').replace('tab_', 'meta_')

            tab_filt_pd.reset_index().to_csv(tab_filt_fp, index=False, sep='\t')
            meta_filt_pd.reset_index().to_csv(meta_filt_fp, index=False, sep='\t')

            datasets_update[dat_filt] = [tab_filt_fp, meta_filt_fp]
            datasets_read_update[dat_filt] = [tab_filt_pd, meta_filt_pd.reset_index()]
            datasets_phylo_update[dat_filt] = datasets_phylo[dat]

            datasets_features_update[dat_filt] = dict(
                gid_feat for gid_feat in datasets_features[dat].items() if gid_feat[1] in tab_filt_pd.index
            )
            qza = tab_filt_fp.replace('.tsv', '.qza')
            cmd = run_import(tab_filt_fp, qza, "FeatureTable[Frequency]")
            sh.write('echo "%s"\n' % cmd)
            sh.write('%s\n' % cmd)
            written += 1

    run_xpbs(out_sh, out_pbs, '%s.fltr' % prjct_nm,
             qiime_env, '4', '4', '1', '100', 'mb', chmod, written,
             '# Filter samples for a min number of %s reads' % thresh)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_features.update(datasets_features_update)
    datasets_phylo.update(datasets_phylo_update)
