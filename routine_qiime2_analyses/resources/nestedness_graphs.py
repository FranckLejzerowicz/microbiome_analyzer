import sys
import pandas as pd
import numpy as np
from os.path import isfile, splitext

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


dat = '<DAT>'
cur_raref = '<CUR_RAREF>'
tab_pd = pd.read_table('<TAB_FP>', index_col=0)
if tab_pd.shape[0] > 1000:
    print('No plotting all %s features...' % tab_pd.shape[0])
else:
    meta_pd = pd.read_table('<META_FP>', dtype={'sample_name': str})
    colors_sample = '<COLORS_SAMPLE>'
    colors_feature = '<COLORS_FEATURE>'
    stats_tax_dat = '<STATS_TAX_DAT>'
    split_taxa_fp = '<SPLIT_TAXA_FP>'
    level = '<LEVEL>'
    collapsed = '<COLLAPSED>'
    nestedness_raref = '<NESTEDNESS_RAREF>'

    if colors_sample:
        meta_pd = meta_pd.rename(columns={meta_pd.columns[0]: 'SAMPLE_ID'})
        meta_pd = meta_pd.set_index('SAMPLE_ID')
        if set(colors_sample).issubset(meta_pd.columns):
            colors_in = list(set(colors_sample) & set(meta_pd.columns))

    tab_pd = np.log10(tab_pd + 1).stack().reset_index().rename(
        columns={'level_1': 'SAMPLE_ID', 0: 'log10_reads'})
    tab_pd = tab_pd.rename(columns={tab_pd.columns[0]: 'OBSERVATION_ID'})

    tax_cols = []
    if colors_feature and split_taxa_fp:
        tax_pd = pd.read_table(split_taxa_fp, index_col=0)
        if level != 'feature':
            if level not in collapsed[stats_tax_dat]:
                sys.exit(0)
            tax_pd = tax_pd.iloc[
                     :, :collapsed[stats_tax_dat][level]
                     ].drop_duplicates()
            tax_pd['OBSERVATION_ID'] = tax_pd.apply(func=lambda x: ';'.join([y for y in x if str(y) != 'nan']), axis=1)
        else:
            tax_pd = tax_pd.reset_index()
            tax_pd = tax_pd.rename(columns={tax_pd.columns[0]: 'OBSERVATION_ID'})

        tax_pd = tax_pd.set_index('OBSERVATION_ID')
        tax_cols = [tax_pd.columns.tolist()[( x -1)] if isinstance(x, int)
                                                      and x not in tax_pd.columns and tax_pd.shape[1] >= x
                    else x for x in colors_feature]
        tax_cols = [x for x in tax_cols if x in tax_pd.columns]
        if tax_cols:
            tax_pd = tax_pd[tax_cols]

    for (group, case), res in nestedness_raref.items():

        # fields_fp = res['fields']
        graphs_fp = res['graph']
        # if not isfile(fields_fp) or not isfile(graphs_fp):
        if not isfile(graphs_fp):
            continue

        graphs = pd.read_csv(graphs_fp, header=0, sep=',', dtype={'SAMPLE_ID': str})
        samples_order = [y for x, y in sorted(dict(graphs[['SAMPLE_RANK', 'SAMPLE_ID']].values).items())][::-1]
        features_order = [y for x, y in sorted(dict(graphs[['OBSERVATION_RANK', 'OBSERVATION_ID']].values).items())][::-1]
        graphs = graphs.merge(tab_pd, on=['SAMPLE_ID', 'OBSERVATION_ID'], how='left')

        matrix = graphs[['OBSERVATION_ID', 'SAMPLE_ID', 'log10_reads']].pivot_table(
            values='log10_reads', index='OBSERVATION_ID', columns='SAMPLE_ID')
        matrix = matrix.loc[features_order, samples_order]

        # fields = [x.strip() for x in open(fields_fp).readlines()]
        cur_meta_pd = meta_pd.loc[matrix.columns.astype(str).tolist(), colors_in].copy()

        # new_meta = 'metadata'
        # if new_meta in cur_meta_pd.columns:
        #     new_meta = 'passed_metadata'
        # if len(fields) > 1:
        #     cur_meta_pd[new_meta] = cur_meta_pd[fields].apply(func=lambda x: '-'.join(x), axis=1)
        cur_meta_leg_hex = cur_meta_pd.apply(func=lambda x: dict(
            zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size).as_hex())
        )).to_dict()
        cur_meta_leg_rgb = cur_meta_pd.apply(func=lambda x: dict(
            zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size))
        )).to_dict()
        cur_meta_pd_hex = cur_meta_pd.apply(func=lambda x: x.map(dict(
            zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size).as_hex())
        )))
        cur_meta_pd_rgb = cur_meta_pd.apply(func=lambda x: x.map(dict(
            zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size))
        )))

        if tax_cols:

            # print()
            # print(tax_cols)
            # print()
            # print(tax_pd)
            # print(tax_pd.index)
            # print()
            # print(matrix)
            # print(matrix.index)

            cur_tax_pd = tax_pd.loc[matrix.index.astype(str).tolist(), :].copy()
            cur_tax_leg_hex = cur_tax_pd.apply(func=lambda x: dict(
                zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size).as_hex())
            )).to_dict()
            cur_tax_leg_rgb = cur_tax_pd.apply(func=lambda x: dict(
                zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size))
            )).to_dict()
            cur_tax_pd_hex = cur_tax_pd.apply(func=lambda x: x.map(dict(
                zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size).as_hex())
            )))
            cur_tax_pd_rgb = cur_tax_pd.apply(func=lambda x: x.map(dict(
                zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size))
            )))

        X, Y = matrix.shape
        graphs_pdf = res['graph_pdf']
        graphs_pdf_complex = '%s_complex.pdf' % splitext(graphs_pdf)[0]
        with PdfPages(graphs_pdf_complex) as pdf:
            fig, ax = plt.subplots(figsize=(20, (20 * Y/X)))
            if tax_cols:
                g = sns.clustermap(
                    matrix, col_cluster=False, row_cluster=False,
                    linewidths=0.1, cmap='coolwarm',
                    row_colors=cur_tax_pd_hex, col_colors=cur_meta_pd_hex,
                    yticklabels=False, xticklabels=False
                )
                simples = [('feature', cur_tax_pd_rgb, cur_tax_leg_rgb),
                           ('sample', cur_meta_pd_rgb, cur_meta_leg_rgb)]
            else:
                g = sns.clustermap(
                    matrix, col_cluster=False, row_cluster=False,
                    linewidths=0.1, cmap='coolwarm',
                    col_colors=cur_meta_pd_hex,
                    yticklabels=False, xticklabels=False
                )
                simples = [('sample', cur_meta_pd_rgb, cur_meta_leg_rgb)]

            n = 0
            N = 1 / cur_meta_pd_hex.columns.size
            for cdx, col in enumerate(cur_meta_pd_hex.columns):
                n += N
                if not cdx:
                    first_leg = g.ax_col_dendrogram.legend(handles=[
                        mpatches.Patch(label=x, color=y) for x, y in cur_meta_leg_hex[col].items()
                    ], loc="upper left")
                    g.ax_col_dendrogram.add_artist(first_leg)
                else:
                    g.ax_col_dendrogram.legend(handles=[
                        mpatches.Patch(label=x, color=y) for x, y in cur_meta_leg_hex[col].items()
                    ], loc="upper left", bbox_to_anchor=(n, 1))

            if tax_cols:
                n = 0
                N = 1 / cur_tax_pd_hex.columns.size
                for cdx, col in enumerate(cur_tax_pd_hex.columns):
                    n += N
                    if not cdx:
                        first_leg = g.ax_row_dendrogram.legend(handles=[
                            mpatches.Patch(label=x, color=y) for x, y in cur_tax_leg_hex[col].items()
                        ], loc="lower left")
                        g.ax_row_dendrogram.add_artist(first_leg)
                    else:
                        g.ax_row_dendrogram.legend(handles=[
                            mpatches.Patch(label=x, color=y) for x, y in cur_tax_leg_hex[col].items()
                        ], loc="lower left", bbox_to_anchor=(0, n))

            g.ax_heatmap.set_xlabel('Samples (sorted by richness)')
            g.ax_heatmap.set_ylabel('Taxon (sorted by prevalence)')

            suptitle = '%s%s' % (dat, cur_raref)
            if case != 'ALL':
                suptitle = suptitle + '\nsamples subset: %s' % case
            if group != '':
                suptitle = suptitle + '\nfeatures subset: %s' % group
            plt.suptitle(suptitle, fontsize=20)
            plt.subplots_adjust(top=.95, hspace=0.3)
            pdf.savefig(bbox_inches='tight', dpi=300)
            plt.close()

        print('[Nestedness] Written:', graphs_pdf_complex)

        graphs_pdf_simples = '%s_simples.pdf' % splitext(graphs_pdf)[0]
        with PdfPages(graphs_pdf_simples) as pdf:
            for (typ, tab, leg) in simples:
                num = tab.columns.size
                fig, axes = plt.subplots(num, 1, figsize=(9, (tab.columns.size * 4)))
                for cdx, col in enumerate(tab.columns):
                    leg_col2rgb = dict((x, leg[col][x]) for idx, x in enumerate(leg[col]))
                    cols_d = tab[col].to_dict()
                    if typ == 'sample':
                        mat = [
                            [
                                [np.nan, np.nan, np.nan, 0.] if str(x[1]) == 'nan' else
                                ([y for y in cols_d[x[0]]] + [1.]) for x in row.reset_index().values
                            ] for r, row in matrix.iterrows()
                        ]
                    if typ == 'feature':
                        mat = [
                            [
                                [np.nan, np.nan, np.nan, 0.] if str(x) == 'nan' else
                                ([y for y in cols_d[r]] + [1.]) for x in row
                            ] for r, row in matrix.iterrows()
                        ]
                    mat = np.array(mat)
                    hands = [mpatches.Patch(label=x, color=y) for x, y in leg_col2rgb.items()]
                    if num == 1:
                        axes.imshow(mat, aspect="auto")
                        axes.legend(handles=hands, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                        axes.set_title(col)
                        axes.set_xlabel('Samples (sorted by richness)')
                        axes.set_ylabel('%s (sorted by prevalence)' % level)
                    else:
                        axes[cdx].imshow(mat, aspect="auto")
                        axes[cdx].legend(handles=hands, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                        axes[cdx].set_title(col)
                        axes[cdx].set_xlabel('Samples (sorted by richness)')
                        axes[cdx].set_ylabel('%s (sorted by prevalence)' % level)
                suptitle = '%s%s (%s coloring)' % (dat, cur_raref, typ)
                top = .93
                if case != 'ALL':
                    top -= 0.03
                    suptitle = suptitle + '\nsamples subset: %s' % case
                if group != '':
                    top -= 0.03
                    suptitle = suptitle + '\nfeatures subset: %s' % group
                plt.suptitle(suptitle, fontsize=15)
                plt.subplots_adjust(top=top, hspace=0.3)
                pdf.savefig(bbox_inches='tight', dpi=300)
                plt.close()

        print('[Nestedness] Written:', graphs_pdf_simples)
