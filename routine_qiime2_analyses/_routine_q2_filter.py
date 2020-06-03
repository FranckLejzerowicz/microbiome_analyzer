# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import yaml
import pandas as pd
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_raref_tab_meta_pds
from routine_qiime2_analyses._routine_q2_cmds import run_import


def import_datasets(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
                    force: bool, prjct_nm: str, qiime_env: str,  chmod: str,
                    noloc: bool, run_params: dict) -> None:
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
    run_xpbs(out_sh, out_pbs, '%s.mprt' % prjct_nm, qiime_env,
             run_params["time"], run_params["n_nodes"], run_params["n_procs"],
             run_params["mem_num"], run_params["mem_dim"],
             chmod, written, '# Import tables to qiime2', None, noloc)


def get_threshs(p_filt_threshs):
    if not isfile(p_filt_threshs):
        print('yaml file for filtering thresholds does not exist:\n%s\nExiting...' % p_filt_threshs)
        sys.exit(0)
    with open(p_filt_threshs) as handle:
        threshs_d = yaml.load(handle, Loader=yaml.FullLoader)
        return threshs_d


def filter_rare_samples(i_datasets_folder: str, datasets: dict, datasets_read: dict,
                        datasets_features: dict, datasets_filt: dict, datasets_phylo: dict,
                        prjct_nm: str, qiime_env: str, p_filt_threshs: str,
                        chmod: str, noloc: bool, run_params: dict) -> None:
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
    threshs_dats = get_threshs(p_filt_threshs)

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
            if dat not in threshs_dats:
                continue
            threshs_d = threshs_dats[dat]
            names = []
            if 'names' in threshs_d:
                names = threshs_d['names']
            thresh_sam = 0
            if 'samples' in threshs_d:
                thresh_sam = threshs_d['samples']
            thresh_feat = 0
            if 'features' in threshs_d:
                thresh_feat = threshs_d['features']

            if not thresh_sam and not thresh_feat:
                print('Filtering threshold(s) of 0 do nothing: skipping...')
                continue
            if not isinstance(thresh_sam, (float, int)) or not isinstance(thresh_feat, (float, int)):
                print('Filtering threshold for %s not a integer/float: skipping...' % dat)
                continue
            if thresh_sam < 0 or thresh_feat < 0:
                print('Filtering threshold must be positive: skipping...')
                continue

            if tab_meta_pds == 'raref':
                tsv, meta = datasets[dat]
                if not isfile(tsv):
                    print('Must have run rarefaction to use it further...\nExiting')
                    sys.exit(0)
                tab_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                # --> datasets_read <--
                # path_pd : indexed with feature name
                # meta_pd : not indexed -> "sample_name" as first column
                datasets_read[dat] = [tab_pd, meta_pd]
            else:
                tab_pd, meta_pd = tab_meta_pds

            if datasets_features[dat] == 'raref':
                datasets_features[dat] = dict(
                    gid_feat for gid_feat in datasets_features[dat].items() if gid_feat[1] in tab_pd.index
                )

            dat_filt = []
            if names:
                dat_filt.append('%srm' % len(names))
            if thresh_sam:
                if thresh_sam > 1:
                    dat_filt.append('minSam%s' % thresh_sam)
                else:
                    dat_filt.append('minSam%s' % str(thresh_sam).replace('.', ''))

            if thresh_feat:
                if thresh_feat > 1:
                    dat_filt.append('minFeat%s' % thresh_feat)
                else:
                    dat_filt.append('minFeat%s' % str(thresh_feat).replace('.', ''))
            dat_filt = '%s_%s' % (dat, '-'.join(dat_filt))
            datasets_filt[dat] = dat_filt
            datasets_filt[dat_filt] = dat
            tab_filt_fp = '%s/data/tab_%s.tsv' % (i_datasets_folder, dat_filt)
            qza = tab_filt_fp.replace('.tsv', '.qza')
            meta_filt_fp = tab_filt_fp.replace(
                '%s/data/' % i_datasets_folder,
                '%s/metadata/' % i_datasets_folder
            ).replace('tab_', 'meta_')
            if isfile(qza) and isfile(meta_filt_fp):
                datasets_update[dat_filt] = [tab_filt_fp, meta_filt_fp]
                tab_filt_pd = pd.read_csv(tab_filt_fp, index_col=0, header=0, sep='\t')
                meta_filt_pd = pd.read_csv(meta_filt_fp, header=0, sep='\t', dtype={'sample_name': str})
                datasets_read_update[dat_filt] = [tab_filt_pd, meta_filt_pd]
                datasets_phylo_update[dat_filt] = datasets_phylo[dat]
                datasets_features_update[dat_filt] = dict(
                    gid_feat for gid_feat in datasets_features[dat].items() if gid_feat[1] in tab_filt_pd.index
                )
                continue

            meta_pd = meta_pd.set_index('sample_name')

            dat_filt = []
            if names:
                dat_filt.append('%srm' % len(names))
                tab_filt_pd = tab_pd[[x for x in tab_pd.columns if x not in names]].copy()
            else:
                tab_filt_pd = tab_pd.copy()

            if thresh_sam:
                if thresh_sam > 1:
                    tab_filt_pd = tab_filt_pd.loc[:, tab_filt_pd.sum(0) >= thresh_sam]
                    dat_filt.append('minSam%s' % thresh_sam)
                else:
                    tab_perc_min = tab_filt_pd.sum(0).mean() * thresh_sam
                    tab_filt_pd = tab_filt_pd.loc[:, tab_filt_pd.sum(0) >= tab_perc_min]
                    dat_filt.append('minSam%s' % str(thresh_sam).replace('.', ''))

            if thresh_feat:
                if thresh_feat > 1:
                    tab_filt_rm = tab_filt_pd < thresh_feat
                    dat_filt.append('minFeat%s' % thresh_feat)
                else:
                    tab_perc = tab_filt_pd/tab_filt_pd.sum(0)
                    tab_filt_rm = tab_perc < thresh_feat
                    dat_filt.append('minFeat%s' % str(thresh_feat).replace('.', ''))
                tab_filt_pd[tab_filt_rm] = 0

            tab_filt_pd = tab_filt_pd.loc[tab_filt_pd.sum(1) > 0, :]
            tab_filt_pd = tab_filt_pd.loc[:, tab_filt_pd.sum(0) > 0]

            dat_filt = '%s_%s' % (dat, '-'.join(dat_filt))
            if tab_filt_pd.shape[0] < 2 or tab_filt_pd.shape[1] < 2:
                print('Filtering too harsh (no more data for %s): skipping...' % dat_filt)
                continue

            meta_filt_pd = meta_pd.loc[tab_filt_pd.columns.tolist()].copy()
            tab_filt_pd.reset_index().to_csv(tab_filt_fp, index=False, sep='\t')
            meta_filt_pd.reset_index().to_csv(meta_filt_fp, index=False, sep='\t')

            datasets_update[dat_filt] = [tab_filt_fp, meta_filt_fp]
            datasets_read_update[dat_filt] = [tab_filt_pd, meta_filt_pd.reset_index()]
            datasets_phylo_update[dat_filt] = datasets_phylo[dat]
            datasets_features_update[dat_filt] = dict(
                gid_feat for gid_feat in datasets_features[dat].items() if gid_feat[1] in tab_filt_pd.index
            )
            cmd = run_import(tab_filt_fp, qza, "FeatureTable[Frequency]")
            sh.write('echo "%s"\n' % cmd)
            sh.write('%s\n' % cmd)
            written += 1
    if written:
        run_xpbs(out_sh, out_pbs, '%s.fltr' % prjct_nm, qiime_env,
                 run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                 run_params["mem_num"], run_params["mem_dim"], chmod, written,
                 '# Filter samples for a min number of %s reads' % p_filt_threshs, None, noloc)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_features.update(datasets_features_update)
    datasets_phylo.update(datasets_phylo_update)
