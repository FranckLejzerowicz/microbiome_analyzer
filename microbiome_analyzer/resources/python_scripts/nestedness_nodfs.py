# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np
from os.path import isfile

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgba
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import seaborn as sns


dat = '<DAT>'
odir = '<ODIR>'
nodfs = '<NODFS>'
coll = '<COLLAPSED>'

arefiles = [x for x in nodfs if isfile(x)]
if arefiles:
    nodfs_pd = pd.concat([pd.read_csv(x, sep='\t')
                          for x in arefiles], sort=False)
    nodfs = nodfs_pd.drop(columns=['GRAPH_ID'])
    nodfs_cols = ['NODF_OBSERVED', 'NODF_NULL_MEAN']
    nodfs = nodfs.set_index([x for x in nodfs.columns if x not in nodfs_cols])
    nodfs = nodfs.stack().reset_index().rename(
        columns={'level_10': 'NODF_TYPE', 0: 'NODF'})
    nodfs_obs = nodfs.loc[nodfs.NODF_TYPE == 'NODF_OBSERVED'].drop(
        columns=['NODF_NULL_STDEV'])
    nodfs_null = nodfs.loc[nodfs.NODF_TYPE == 'NODF_NULL_MEAN'].copy()

    nodfs_null_up = nodfs_null.copy()
    nodfs_null_down = nodfs_null.copy()
    nodfs_null_up["NODF"] = nodfs_null_up.NODF + nodfs_null_up.NODF_NULL_STDEV
    nodfs_null_down["NODF"] = nodfs_null_up.NODF - nodfs_null_up.NODF_NULL_STDEV
    nodfs_and_null = pd.concat([nodfs_obs,
                                nodfs_null_up,
                                nodfs_null_down], sort=False)
    nodfs_metas = nodfs_and_null.loc[nodfs_and_null.COMPARISON != 'overall']
    nodfs_leg = dict(zip(nodfs_metas.COMPARISON.unique(), sns.color_palette(
                palette='Set1', n_colors=nodfs_metas.COMPARISON.unique().size)))
    nodfs_leg['overall'] = (0, 0, 0)

    comparison_raref = 'COMPARISON'
    if sum(nodfs_and_null.RAREF.fillna('NaN').str.contains('raref')):
        comparison_raref = 'COMPARISON_RAREF'
        nodfs_and_null[comparison_raref] = nodfs_and_null[
            ['COMPARISON', 'RAREF']
        ].fillna('').apply(func=lambda x: ''.join(x), axis=1)
        rarefs = sorted(set([
            x.split('_raref')[-1] for x in nodfs_and_null.RAREF.unique()]))
        alphas = np.linspace(0.4, 1, len(
            nodfs_and_null.RAREF.fillna('').unique()))
        nodfs_leg = dict(
            ('%s_raref%s' % (x, raref), to_rgba(y, alphas[rdx]))
            for x, y in nodfs_leg.items()
            for rdx, raref in enumerate(rarefs))

    nodfs_and_null['MODE_METADATA'] = nodfs_and_null[
        ['MODE', 'METADATA']
    ].fillna('').apply(func=lambda x: ' - '.join([y for y in x if y]), axis=1)
    levels = dict(
        (y[0], x) for x, y in enumerate(sorted(coll.items(), key=lambda i: i[1])
        ) if y[0] in nodfs_and_null.LEVEL.unique())
    levels['feature'] = 1
    if levels:
        levels['feature'] = max(levels.values()) + 1
    nodfs_and_null['LEVEL_SORT'] = [levels[x] for x in nodfs_and_null['LEVEL']]
    nodfs_and_null = nodfs_and_null.sort_values(
        ['LEVEL_SORT', 'NODF_TYPE'], ascending=True)
    nodfs_and_null = nodfs_and_null.drop(columns=['RAREF', 'MODE', 'METADATA'])
    table_fpo = '%s/nodfs.tsv' % odir
    nodfs_and_null.to_csv(table_fpo, index=False, sep='\t')

    nodfs_pdf = '%s/nodfs.pdf' % odir
    with PdfPages(nodfs_pdf) as pdf:
        for (mode_meta, null), nodfs_ in nodfs_and_null.groupby(
                ['MODE_METADATA', 'NULL']):
            nodfs = nodfs_.drop(columns=['MODE_METADATA'])
            nodfs['SORT'] = nodfs['LEVEL_SORT'].astype(str) + nodfs['LEVEL']
            SUBSETS = sorted([
                x for x in nodfs.SAMPLE_SUBSET.unique().tolist()if x != 'ALL'])
            if 'ALL' in nodfs.SAMPLE_SUBSET.unique():
                SUBSETS = ['ALL'] + SUBSETS
            nodfs.index = range(nodfs.shape[0])
            g = sns.relplot(
                data=nodfs, x="SORT", y="NODF", col="SAMPLE_SUBSET", col_wrap=2,
                col_order=SUBSETS, hue=comparison_raref, palette=nodfs_leg,
                style='NODF_TYPE', kind="line", err_style="bars", ci='sd',
                style_order=['NODF_OBSERVED', 'NODF_NULL_MEAN'],
                facet_kws={'sharey': False, 'sharex': False})
            axes = g.axes
            for sdx, SUBSET in enumerate(SUBSETS):
                ticks = axes[sdx].get_xticks()
                axes[sdx].xaxis.set_major_locator(mticker.FixedLocator(ticks))
                axes[sdx].set_xticklabels([x[0] for x in sorted(
                    levels.items(), key=lambda item: item[1])], rotation=45)
                axes[sdx].set_xlabel(SUBSET)
                SUBSET_pd = nodfs.loc[(nodfs.SAMPLE_SUBSET == SUBSET) &
                                      (nodfs.NODF_TYPE == 'NODF_OBSERVED')
                ].drop(columns=['SAMPLE_SUBSET'])
                for (LEVEL_SORT, NODF, PVALUE, COMP) in SUBSET_pd[
                    ['LEVEL_SORT', 'NODF', 'PVALUE', comparison_raref]
                ].values:
                    if PVALUE:
                        axes[sdx].text(LEVEL_SORT, NODF+0.005, PVALUE,
                                       color=nodfs_leg[COMP], fontsize=5)
            suptitle = '%s - %s (%s)' % (dat, mode_meta, null)
            plt.suptitle(suptitle, fontsize=12)
            plt.subplots_adjust(top=.90, hspace=0.3)
            pdf.savefig(bbox_inches='tight', dpi=300)
            plt.close()

    print('[Nestedness] Written:', nodfs_pdf)
