import pandas as pd
from os.path import isfile

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns


dat = 'DAT'
odir = 'ODIR'
nodfs = NODFS
collapsed = COLLAPSED

arefiles = [x for x in nodfs if isfile(x)]
if arefiles and dat in collapsed:
    nodfs_pd = pd.concat([pd.read_table(x) for x in arefiles])
    levels = collapsed[dat]
    levels['feature'] = max(levels.values()) + 1
    nodfs = nodfs_pd.drop(columns=['GRAPH_ID'])
    nodfs_cols = ['NODF_OBSERVED', 'NODF_NULL_MEAN']
    nodfs = nodfs.set_index([x for x in nodfs.columns if x not in nodfs_cols])
    nodfs = nodfs.stack().reset_index().rename(columns={'level_11': 'NODF_TYPE', 0: 'NODF'})
    nodfs_obs = nodfs.loc[nodfs.NODF_TYPE == 'NODF_OBSERVED'].drop(columns=['NODF_NULL_STDEV'])
    nodfs_null = nodfs.loc[nodfs.NODF_TYPE == 'NODF_NULL_MEAN'].copy()

    nodfs_null_up = nodfs_null.copy()
    nodfs_null_down = nodfs_null.copy()
    nodfs_null_up["NODF"] = nodfs_null_up.NODF + nodfs_null_up.NODF_NULL_STDEV
    nodfs_null_down["NODF"] = nodfs_null_up.NODF - nodfs_null_up.NODF_NULL_STDEV
    nodfs_null = pd.concat([nodfs_obs, nodfs_null_up, nodfs_null_down])

    nodfs_metas = nodfs_null.loc[nodfs_null.COMPARISON != 'overall']
    nodfs_leg = dict(zip(nodfs_metas.COMPARISON.unique(), sns.color_palette(
                palette='Set1', n_colors=nodfs_metas.COMPARISON.unique().size).as_hex()))
    nodfs_leg['overall'] = '#000000'

    nodfs_null['SUBSET'] = nodfs_null[
        ['FEATURE_SUBSET', 'SAMPLE_SUBSET']].apply(func=lambda x: ''.join(x), axis=1)
    comparison_raref = 'COMPARISON_RAREF'
    nodfs_null[comparison_raref] = nodfs_null[
        ['COMPARISON', 'RAREF']].apply(func=lambda x: ''.join(x), axis=1)
    if not sum(nodfs_null.COMPARISON_RAREF.str.contains('raref')):
        comparison_raref = 'COMPARISON'

    nodfs_null['MODE_METADATA'] = nodfs_null[
        ['MODE', 'METADATA']].fillna('').apply(func=lambda x: ' - '.join([y for y in x if y]), axis=1)
    nodfs_null['LEVEL_SORT'] = [levels[x] for x in nodfs_null['LEVEL']]
    nodfs_null = nodfs_null.sort_values(['LEVEL_SORT', 'NODF_TYPE'], ascending=True)

    nodfs_null = nodfs_null.drop(columns=[
        'FEATURE_SUBSET', 'SAMPLE_SUBSET', 'RAREF', 'MODE', 'METADATA'
    ])

    table_fpo = '%s/nodfs.tsv' % odir
    nodfs_null.to_csv(table_fpo, index=False, sep='\t')
    nodfs_pdf = '%s/nodfs.pdf' % odir
    with PdfPages(nodfs_pdf) as pdf:
        for (mode_metadata), nodfs_null_mode_ in nodfs_null.groupby('MODE_METADATA'):
            nodfs_null_mode = nodfs_null_mode_.drop(columns=['MODE_METADATA'])
            SUBSETS = ['ALL'] + sorted([x for x in nodfs_null_mode.SUBSET.unique().tolist() if x != 'ALL'])
            g = sns.relplot(
                data=nodfs_null_mode, x="LEVEL", y="NODF",
                col="SUBSET", col_wrap=2, col_order=SUBSETS,
                hue=comparison_raref, palette=nodfs_leg,
                style='NODF_TYPE', style_order=['NODF_OBSERVED', 'NODF_NULL_MEAN'],
                kind="line", err_style="bars", ci='sd')

            axes = g.axes
            X = 0
            for SUBSET in SUBSETS:
                SUBSET_pd = nodfs_null_mode.loc[
                    (nodfs_null_mode.SUBSET == SUBSET) &
                    (nodfs_null_mode.NODF_TYPE == 'NODF_OBSERVED')
                ].drop(columns=['SUBSET'])

                for (LEVEL_SORT, NODF, PVALUE, COMP) in SUBSET_pd[
                    ['LEVEL_SORT', 'NODF', 'PVALUE', comparison_raref]].values:
                    if PVALUE:
                        axes[X].text(
                            (LEVEL_SORT-1), NODF+0.005, PVALUE,
                            color=nodfs_leg[COMP], fontsize=5
                        )
                X += 1

            suptitle = '%s - %s' % (dat, mode_metadata)
            plt.suptitle(suptitle, fontsize=12)
            plt.subplots_adjust(top=.90, hspace=0.3)
            pdf.savefig(bbox_inches='tight', dpi=300)
            plt.close()
