# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import itertools
import numpy as np
import pandas as pd

import plotly
import plotly.graph_objs as go

from microbiome_analyzer._io_utils import convert_to_biom


def no_filtering(
        dat: str,
        thresh_sam: int,
        thresh_feat: int) -> bool:
    """Checks whether to skip filtering or not.

    Parameters
    ----------
    dat : str
        Dataset name
    thresh_sam : int
        Samples threshold
    thresh_feat : int
        Features threshold

    Returns
    -------
    skip : bool
        Whether to skip filtering or not
    """
    skip = False
    if not thresh_sam and not thresh_feat:
        print('Filtering threshold(s) of 0 do nothing: skipping...')
        skip = True
    thresh_sam_is_numeric = isinstance(thresh_sam, (float, int))
    thresh_feat_is_numeric = isinstance(thresh_feat, (float, int))
    if not thresh_sam_is_numeric or not thresh_feat_is_numeric:
        print('Filtering threshold for %s not a '
              'integer/float: skipping...' % dat)
        skip = True
    if thresh_sam < 0 or thresh_feat < 0:
        print('Filtering threshold must be positive: skipping...')
        skip = True
    return skip


def get_dat_filt(
        dat: str,
        names: list,
        thresh_sam: int,
        thresh_feat: int) -> str:
    """Get a build-up new name for
    the filtered version of a dataset.

    Parameters
    ----------
    dat : str
        Dataset name
    names : list
        Name for the threshold
    thresh_sam : int
        Samples threshold
    thresh_feat : int
        Features threshold

    Returns
    -------
    dat_filt : str
        New dataset name for the filtered version
    """
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
    return dat_filt


def filtering_names(
        names: list,
        biom_filt: biom.Table) -> str:
    """
    Parameters
    ----------
    names : list
        Name for the threshold
    biom_filt : biom.Table
        BIOM table

    Returns
    -------
    filt_cmds : str
        Text about what was filtered.
    """
    filt_cmds = ''
    if names:
        ids = set(biom_filt.ids(axis='sample'))
        names_in = list(ids & set(names))
        names_keep = list(ids - set(names))
        filt_cmds += '- [sample names] Removed %s samples:\n' % len(names_in)
        for name_in in names_in:
            filt_cmds += '  * %s\n' % name_in
        biom_filt.filter(ids_to_keep=names_keep)
    return filt_cmds


def filtering_samples(
        thresh_sam: int,
        biom_filt: biom.Table) -> str:
    """
    Parameters
    ----------
    thresh_sam : int
        Samples threshold
    biom_filt : biom.Table
        BIOM table

    Returns
    -------
    filt_cmds : str
        Text about what was filtered.
    """
    filt_cmds = ''
    if thresh_sam:
        samples = biom_filt.ids(axis='sample')
        biom_filt_sum = pd.Series(biom_filt.sum(axis='sample'), index=samples)
        if thresh_sam > 1:
            thresh = thresh_sam
        else:
            thresh = biom_filt_sum.mean() * thresh_sam
        to_drop = biom_filt_sum[biom_filt_sum < thresh]
        if to_drop.size:
            filt_cmds += '- Removed %s samples (min %s reads):\n' % (
                len(to_drop), str(thresh))
            for to_d in to_drop.index:
                filt_cmds += '  * %s (%s reads)\n' % (to_d, biom_filt_sum[to_d])
            to_keep = list(set(samples) - set(to_drop.index))
            biom_filt.filter(ids_to_keep=to_keep)
    return filt_cmds


def filtering_features(
        thresh_feat: int,
        biom_filt: biom.Table):
    """

    Parameters
    ----------
    thresh_feat : int
        Features threshold
    biom_filt : biom.Table
        BIOM table

    Returns
    -------
    filt_cmds : str
        Text about what was filtered.
    biom_filt : biom.Table
        Feature filtered BIOM table
    """
    filt_cmds = ''
    biom_pd = biom_filt.to_dataframe(dense=True)
    if thresh_feat:
        if thresh_feat > 1:
            biom_rm = biom_pd < thresh_feat
            rm_thresh_feat = thresh_feat
        else:
            tab_perc = biom_pd / biom_pd.sum(0)
            biom_rm = tab_perc < thresh_feat
            rm_thresh_feat = '%s%s' % ((thresh_feat * 100), '%')
        filt_cmds += '- Removed %s features (min %s reads):\n' % (
            biom_rm.any(1).sum(), rm_thresh_feat)
        tab_filt_rmd = biom_rm.sum(1).to_dict()
        filt_cmds += '- Number of samples from which each feature if removed:\n'
        for feat, n_rm in tab_filt_rmd.items():
            filt_cmds += '  * %s\t%s\n' % (feat, n_rm)
        biom_pd[biom_rm] = 0
    biom_filt = convert_to_biom(biom_pd)
    return filt_cmds, biom_filt


def filtering_thresholds(
        names: list,
        thresh_sam: int,
        thresh_feat: int,
        biom_table: biom.Table) -> tuple:
    """

    Parameters
    ----------
    names : list
        Name for the threshold
    thresh_sam : int
        Samples threshold
    thresh_feat : int
        Features threshold
    biom_table : biom.Table
        BIOM table

    Returns
    -------
    biom_filt : biom.Table
        Output filtered feature table
    filt_cmds : str
        Command used to filter
    """
    filt_cmds = ''
    biom_filt = biom_table.copy()
    filt_cmds += filtering_names(names, biom_filt)
    filt_cmds += filtering_samples(thresh_sam, biom_filt)
    filt_cmd, biom_filt = filtering_features(thresh_feat, biom_filt)
    filt_cmds += filt_cmd
    biom_filt.remove_empty()
    return biom_filt, filt_cmds


def harsh_filtering(
        dat_filt: str,
        data_filt_biom: biom.Table) -> bool:
    """

    Parameters
    ----------
    dat_filt : str
        New dataset name for the filtered version
    data_filt_biom : biom.Table
        Filtered feature table

    Returns
    -------
    skip : bool
        Whether to skip a too harsh filtering
    """
    skip = False
    if data_filt_biom.shape[0] < 10 or data_filt_biom.shape[1] < 2:
        print('Filtering too harsh (no more data for %s): '
              'skipping...' % dat_filt)
        skip = True
    return skip


def filter_mb_table(preval: str, abund: str, tsv_pd: pd.DataFrame,
                    do_res: bool = False) -> (pd.DataFrame, list):
    preval = float(preval)
    abund = float(abund)
    if abund:
        new_cols = {}
        cur_index = tsv_pd.index.tolist()
        cur_columns = tsv_pd.columns.tolist()
        for col in cur_columns:
            cur_col = tsv_pd[col]
            min_thresh = min([x for x in cur_col if x > 0]) * abund
            new_col = [x if x > min_thresh else 0 for x in cur_col]
            new_cols[col] = new_col
        tsv_pd = pd.DataFrame(new_cols, index=cur_index, columns=cur_columns)
        tsv_pd = tsv_pd[tsv_pd.sum(1) > 1]
    if preval:
        if preval < 1:
            n_perc = tsv_pd.shape[1] * preval
        else:
            n_perc = preval
        tsv_pd = tsv_pd.loc[tsv_pd.astype(bool).sum(1) >= n_perc, :]
    tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
    res = []
    if do_res:
        res = [preval, abund, tsv_pd.shape[0], tsv_pd.shape[1]]
    return tsv_pd, res


def filter_non_mb_table(preval: str, abund: str, tsv_pd: pd.DataFrame,
                        do_res: bool = False) -> (pd.DataFrame, list):
    res = []
    preval = float(preval)
    abund = float(abund)
    if preval + abund == 0:
        if do_res:
            res = [0, 0, tsv_pd.shape[0], tsv_pd.shape[1]]
        return tsv_pd, res
    tsv_filt_pd = tsv_pd.copy()
    # get the min number of samples based on prevalence percent
    if preval < 1:
        n_perc = tsv_pd.shape[1] * preval
    else:
        n_perc = preval
    # abundance filter in terms of min reads counts
    if abund < 1:
        # change counts to percent per sample
        tsv_pd_perc = tsv_filt_pd / tsv_filt_pd.sum()
        # get percent of features across sample
        tsv_pd_perc_sum = tsv_filt_pd.sum(1) / tsv_filt_pd.sum(1).sum()
    else:
        # do not change counts
        tsv_pd_perc = tsv_filt_pd.copy()
        # get sum of features counts across samples
        tsv_pd_perc_sum = tsv_filt_pd.sum(1)
    abund_mode = 'sample'
    if abund_mode == 'sample':
        # remove features from feature table that are not present in enough
        # samples with the minimum number/percent of reads in these samples
        # - get the table's cells which percent value is > [0-1] of "abund"
        # - counts the number of sample occurrences for these cells
        # - only keep the features§ rows that have > "n_perc" such occurrences
        tsv_filt_pd = tsv_filt_pd.loc[(tsv_pd_perc > abund).sum(1) > n_perc, :]
    elif abund_mode == 'dataset':
        tsv_filt_pd = tsv_filt_pd.loc[tsv_pd_perc_sum > abund, :]
    elif abund_mode == 'both':
        tsv_filt_pd = tsv_filt_pd.loc[(tsv_pd_perc > abund).sum(1) > n_perc, :]
        if abund < 1:
            fil_pd_perc_sum = tsv_filt_pd.sum(1) / tsv_filt_pd.sum(1).sum()
        else:
            fil_pd_perc_sum = tsv_filt_pd.sum(1)
        tsv_filt_pd = tsv_filt_pd.loc[fil_pd_perc_sum > abund, :]
    else:
        raise Exception('"%s" mode not recognized' % abund_mode)
    # - remove empty rows and empty columns
    tsv_filt_pd = tsv_filt_pd.loc[tsv_filt_pd.sum(1) > 0,
                                  tsv_filt_pd.sum(0) > 0]
    if do_res:
        res = [preval, abund, tsv_filt_pd.shape[0], tsv_filt_pd.shape[1]]
    return tsv_filt_pd, res


def filter_3d(dat, pv, ab, defaults, biom_, targeted, out_dir):
    prevals = defaults['preval'][pv]
    abunds = defaults['abund'][ab]

    tab_pd = biom_.to_dataframe(dense=True)
    res = []
    for (preval, abund) in itertools.product(*[sorted(prevals),
                                               sorted(abunds)]):
        # if mb:
        #     tsv_pd, cur_res = filter_mb_table(
        #         preval, abund, tsv_pd, True)
        # else:

        if float(preval) < 1:
            pv_ = 'fraction'
        else:
            pv_ = 'count'
        if float(abund) < 1:
            ab_ = 'fraction'
        else:
            ab_ = 'count'

        tsv_pd, cur_res = filter_non_mb_table(preval, abund, tab_pd, True)
        if (preval, abund) in targeted.get('global', []):
            cur_res.append(1)
        elif (preval, abund) in targeted.get(dat, []):
            cur_res.append(1)
        else:
            cur_res.append(0)
        res.append(cur_res)

    res_pd = pd.DataFrame(res, columns=[
        'preval_filt', 'abund_filt', 'features', 'samples', 'data'])
    res_pd.preval_filt = res_pd.preval_filt.astype(float)
    res_pd.abund_filt = res_pd.abund_filt.astype(float)
    res_pd['features'] = np.log10(res_pd['features'] + 1)

    x = res_pd.preval_filt.unique()
    y = res_pd.abund_filt.unique()
    X, Y = np.meshgrid(x, y)
    Z = res_pd.features.values.reshape(X.shape, order='f')

    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='Min abundance (%s)' % ab_),
            yaxis=dict(title='Min prevalence (%s)' % pv_),
            zaxis=dict(title='log10(features)')),
        autosize=True,
        width=700, height=700,
        title="Filtering process: %s" % dat,
        margin=dict(l=65, r=50, b=65, t=90))

    fig = go.Figure(
        data=[go.Surface(x=Y, y=X, z=Z,
                         colorscale='Viridis',
                         reversescale=True)],
        layout=layout)
    fig.update_traces(
        contours_z=dict(show=True, usecolormap=True,
                        highlightcolor="limegreen", project_z=True))
    fig.add_scatter3d(y=X.flatten(), x=Y.flatten(), z=Z.flatten(),
                      mode='markers', marker=dict(size=4, color='black'))
    res_data_pd = res_pd.loc[(res_pd.data == 1)].copy()

    x = res_data_pd.preval_filt.unique()
    y = res_data_pd.abund_filt.unique()
    X, Y = np.meshgrid(x, y)
    Z = res_data_pd.features.values.reshape(X.shape, order='f')
    fig.add_scatter3d(
        y=X.flatten(), x=Y.flatten(), z=Z.flatten(),
        mode='markers', marker=dict(size=6, color='red'))
    html_fo = '%s/prevalence-as-%s_abundance-as-%s.html' % (out_dir, pv_, ab_)
    print(' -> Written:', html_fo)
    plotly.offline.plot(fig, filename=html_fo, auto_open=False)
