# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import yaml
import itertools
import numpy as np
import pandas as pd
from os.path import isfile, splitext

import plotly
import plotly.graph_objs as go

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder, get_raref_tab_meta_pds, get_raref_table, simple_chunks,
    get_analysis_folder, filter_mb_table, filter_non_mb_table)
from routine_qiime2_analyses._routine_q2_cmds import run_import
from routine_qiime2_analyses._routine_q2_mmvec import get_mmvec_dicts
from routine_qiime2_analyses._routine_q2_songbird import get_songbird_dicts


def import_datasets(i_datasets_folder: str, datasets: dict, datasets_phylo: dict,
                    force: bool, prjct_nm: str, qiime_env: str,  chmod: str,
                    noloc: bool, run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
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
    job_folder2 = get_job_folder(i_datasets_folder, 'import_tables/chunks')

    to_chunk = []
    main_written = 0
    run_pbs = '%s/0_run_import_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds_ in datasets.items():
            written = 0
            out_sh = '%s/0_run_import_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for tsv_meta_pds in tsv_meta_pds_:  # REMOVE IF FIXED NOT KEPT
                    tsv, meta = tsv_meta_pds
                    qza = '%s.qza' % splitext(tsv)[0]
                    if datasets_phylo[dat][1]:
                        cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                        cur_sh.write('echo "%s"\n' % cmd)
                        cur_sh.write('%s\n' % cmd)
                        written += 1
                    elif force or not isfile(qza):
                        cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                        cur_sh.write('echo "%s"\n' % cmd)
                        cur_sh.write('%s\n' % cmd)
                        written += 1
            if written:
                main_written += 1
                to_chunk.append(out_sh)
                if not chunkit:
                    run_xpbs(out_sh, out_pbs, '%s.mprt.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                             run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                             run_params["mem_num"], run_params["mem_dim"],
                             chmod, written, 'single', o, noloc, jobs)
    if to_chunk:
        simple_chunks(run_pbs, job_folder2, to_chunk, 'imports',
                      prjct_nm, run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                      run_params["mem_num"], run_params["mem_dim"],
                      qiime_env, chmod, noloc, jobs, chunkit, None)

    if main_written:
        print_message('# Import tables to qiime2', 'sh', run_pbs, jobs)


def get_threshs(p_filt_threshs):
    if not isfile(p_filt_threshs):
        print('yaml file for filtering thresholds does not exist:\n%s\nExiting...' % p_filt_threshs)
        sys.exit(0)
    with open(p_filt_threshs) as handle:
        try:
            threshs_d = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            threshs_d = yaml.load(handle)
        return threshs_d


def deleted_non_filt(datasets: dict, datasets_read: dict, datasets_features: dict,
                     datasets_phylo: dict, datasets_rarefs: dict, taxonomies: dict,
                     datasets_filt: dict, datasets_filt_map: dict):
    for d in [datasets, datasets_read, datasets_features,
              datasets_phylo, datasets_rarefs, taxonomies]:
        to_delete = []
        for dat in d:
            if dat not in datasets_filt_map and dat in datasets_filt:
                to_delete.append(dat)
        for delete in to_delete:
            d.pop(delete)
        break


def filter_rare_samples(i_datasets_folder: str, datasets: dict, datasets_read: dict,
                        datasets_features: dict, datasets_rarefs: dict, datasets_filt: dict,
                        datasets_filt_map: dict, datasets_phylo: dict, prjct_nm: str,
                        qiime_env: str, p_filt_threshs: str, chmod: str, noloc: bool,
                        run_params: dict, filt_raref: str, jobs: bool, chunkit: int) -> None:
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
    out_sh = '%s/1_run_import_filtered_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
    out_pbs = '%s.pbs' % splitext(out_sh)[0]
    to_chunk = []
    with open(out_sh, 'w') as sh:
        for dat, tab_meta_pds_ in datasets_read.items():
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
            datasets_filt_map[dat_filt] = dat
            datasets_rarefs[dat_filt] = ['']
            tab_filt_fp = '%s/data/tab_%s.tsv' % (i_datasets_folder, dat_filt)
            qza = tab_filt_fp.replace('.tsv', '.qza')
            meta_filt_fp = tab_filt_fp.replace(
                '%s/data/' % i_datasets_folder,
                '%s/metadata/' % i_datasets_folder
            ).replace('tab_', 'meta_')
            if isfile(qza) and isfile(meta_filt_fp):
                # datasets_update[dat_filt] = [tab_filt_fp, meta_filt_fp]
                datasets_update[dat_filt] = [[tab_filt_fp, meta_filt_fp]]
                tab_filt_pd = pd.read_csv(tab_filt_fp, index_col=0, header=0, sep='\t')
                with open(meta_filt_fp) as f:
                    for line in f:
                        break
                meta_filt_pd = pd.read_csv(meta_filt_fp, header=0, sep='\t',
                                           dtype={line.split('\t')[0]: str},
                                           low_memory=False)
                # datasets_read_update[dat_filt] = [tab_filt_pd, meta_filt_pd]
                datasets_read_update[dat_filt] = [[tab_filt_pd, meta_filt_pd]]
                datasets_phylo_update[dat_filt] = datasets_phylo[dat]
                datasets_features_update[dat_filt] = dict(
                    gid_feat for gid_feat in datasets_features[dat].items() if gid_feat[1] in tab_filt_pd.index
                )
                continue

            for (tab_pd, meta_pd) in tab_meta_pds_:
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

                # datasets_update[dat_filt] = [tab_filt_fp, meta_filt_fp]
                datasets_update[dat_filt] = [[tab_filt_fp, meta_filt_fp]]
                # datasets_read_update[dat_filt] = [tab_filt_pd, meta_filt_pd.reset_index()]
                datasets_read_update[dat_filt] = [[tab_filt_pd, meta_filt_pd.reset_index()]]
                datasets_phylo_update[dat_filt] = datasets_phylo[dat]
                datasets_features_update[dat_filt] = dict(
                    gid_feat for gid_feat in datasets_features[dat].items() if gid_feat[1] in tab_filt_pd.index
                )
                cmd = run_import(tab_filt_fp, qza, "FeatureTable[Frequency]")
                sh.write('echo "%s"\n' % cmd)
                sh.write('%s\n' % cmd)
                written += 1
    if written:
        run_xpbs(out_sh, out_pbs, '%s.fltr%s' % (prjct_nm, filt_raref), qiime_env,
                 run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                 run_params["mem_num"], run_params["mem_dim"], chmod, written,
                 '# Filter samples for a min number of %s reads' % p_filt_threshs,
                 None, noloc, jobs)

    datasets.update(datasets_update)
    datasets_read.update(datasets_read_update)
    datasets_features.update(datasets_features_update)
    datasets_phylo.update(datasets_phylo_update)


def get_filt3d_params(p_config, analysis):
    if analysis == 'mmvec':
        mmvec_dicts = get_mmvec_dicts(p_config)
        filtering = mmvec_dicts[1]
    else:
        songbird_dicts = get_songbird_dicts(p_config)
        filtering = songbird_dicts[1]
    return filtering


def explore_filtering(i_datasets_folder, datasets, datasets_read,
                      datasets_filt, filtering, p_filt3d_config):

    if p_filt3d_config and isfile(p_filt3d_config):
        with open(p_filt3d_config) as handle:
            try:
                explor = yaml.load(handle, Loader=yaml.FullLoader)
            except AttributeError:
                explor = yaml.load(handle)

        defaults = (
            {'prevalCount': set(explor['prevalCount']),
             'prevalPercent': set(explor['prevalPercent'])},
            {'abundCount': set(explor['abundCount']),
             'abundPercent': set(explor['abundPercent'])})
    else:
        defaults = (
            {'prevalCount': {0, 1, 2, 3, 5, 10, 20, 30},
             'prevalPercent': {0, 0.01, 0.02, 0.05, 0.1}},
            {'abundCount': {0, 1, 2, 5, 10, 100, 1000},
             'abundPercent': {0, .001, .01, .02, .03, .05, .1}})

    scales = {}
    currents = {}
    for pair, filt_d in filtering.items():
        for filt, dat_preval_abund in filt_d.items():
            for (dat, mb), preval_abund in dat_preval_abund.items():

                preval, abund = map(float, preval_abund)

                if preval == 0 and abund == 0:
                    continue
                if preval < 1:
                    preval_label = 'prevalPercent'
                else:
                    preval_label = 'prevalCount'

                if abund < 1:
                    abund_label = 'abundPercent'
                else:
                    abund_label = 'abundCount'

                if (preval_label, abund_label) not in scales:
                    scales[(preval_label, abund_label)] = {}
                if (dat, mb) not in scales[(preval_label, abund_label)]:
                    currents.setdefault((dat, mb), []).append((preval, abund))
                    scales[(preval_label, abund_label)][(dat, mb)] = (
                        defaults[0][preval_label], defaults[1][abund_label])
                scales[(preval_label, abund_label)][(dat, mb)][0].add(preval)
                scales[(preval_label, abund_label)][(dat, mb)][1].add(abund)

    print("scales")
    print(scales)
    for (preval_label, abund_label), dats_d in scales.items():
        out_dir = get_analysis_folder(
            i_datasets_folder, 'filter3D/scale_%s_%s' % (preval_label, abund_label))
        for (dat_, mb), prevals_abunds in dats_d.items():
            if dat_ in datasets_filt:
                dat = datasets_filt[dat_]
            else:
                dat = dat_
            if dat not in datasets:
                if '__raref' in dat:
                    split = dat.split('__raref')
                    dat = '__raref'.join(split[:-1])
                    raref = '_raref%s' % '__raref'.join(split[-1:])
                    if dat in datasets_filt:
                        dat = datasets_filt[dat]
                    tsv_pd_, meta_pd_, meta = get_raref_table(dat, raref, i_datasets_folder, 'filter3D')
                    if not tsv_pd_.shape[0]:
                        continue
                    dat = '%s_%s' % (dat, raref)
                else:
                    print('dataset "%s" not found...' % dat)
                    continue
            elif not isinstance(datasets_read[dat][0], pd.DataFrame) and datasets_read[dat][0] == 'raref':
                tsv, meta = datasets[dat]
                if not isfile(tsv):
                    print('Must have run rarefaction to use it further...\nExiting')
                    sys.exit(0)
                tsv_pd_, meta_pd_ = get_raref_tab_meta_pds(meta, tsv)
                datasets_read[dat] = [[tsv_pd_, meta_pd_]]
            else:
                tsv_pd_, meta_pd_ = datasets_read[dat][0]

            res = []
            rdx = 0
            prevals, abunds = prevals_abunds
            for (preval, abund) in itertools.product(*[sorted(prevals), sorted(abunds)]):
                rdx += 1
                tsv_pd = tsv_pd_.loc[tsv_pd_.sum(1) > 0, :].copy()
                tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
                if mb:
                    tsv_pd, cur_res = filter_mb_table(preval, abund, tsv_pd, True)
                else:
                    tsv_pd, cur_res = filter_non_mb_table(preval, abund, tsv_pd, True)
                if (preval, abund) in currents[(dat, mb)]:
                    cur_res.append(1)
                else:
                    cur_res.append(0)
                res.append(cur_res)
            res_pd = pd.DataFrame(res, columns=['preval_filt', 'abund_filt', 'features',
                                                'samples', 'data'])
            res_pd['features'] = np.log10(res_pd['features']+1)
            x = res_pd.preval_filt.unique()
            y = res_pd.abund_filt.unique()
            X, Y = np.meshgrid(x, y)
            Z = res_pd.features.values.reshape(X.shape, order='f')

            layout = go.Layout(
                scene=dict(
                    xaxis=dict(title=abund_label),
                    yaxis=dict(title=preval_label),
                    zaxis=dict(title='log10(features)')),
                autosize=True,
                width=700, height=700,
                title="Filtering process",
                margin=dict(l=65, r=50, b=65, t=90))
            fig = go.Figure(
                data=[
                    go.Surface(
                        x=Y, y=X, z=Z,
                        colorscale='Viridis',
                        reversescale=True)
                ],
                layout=layout
            )
            fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                              highlightcolor="limegreen", project_z=True))
            fig.add_scatter3d(
                y=X.flatten(), x=Y.flatten(), z=Z.flatten(),
                mode='markers', marker=dict(size=4, color='black'))
            res_data_pd = res_pd.loc[(res_pd.data == 1)].copy()
            x = res_data_pd.preval_filt.unique()
            y = res_data_pd.abund_filt.unique()
            X, Y = np.meshgrid(x, y)
            Z = res_data_pd.features.values.reshape(X.shape, order='f')
            fig.add_scatter3d(
                y=X.flatten(), x=Y.flatten(), z=Z.flatten(),
                mode='markers', marker=dict(size=6, color='red'))
            html_fo = '%s/%s_%s.html' % (out_dir, dat, mb)
            print(' -> Written:', html_fo)
            plotly.offline.plot(fig, filename=html_fo, auto_open=False)
