import os
from os.path import isfile, splitext
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

metrics = METRICS

perms = []
rt = 'ROUTINE_FOLDER/qiime/permanova/DATASET'
for root, dirs, files in os.walk(rt):
    for fil in files:
        if fil.endswith('.html'):
            path = root + '/' + fil
            for metric in metrics:
                if metric in fil:
                    break
            dataset = fil.split('_%s' % metric)[0].split('tab_')[-1]
            if fil.endswith('permanova.html'):
                method = 'permanova'
                subset_test = fil.split('_%s__' % metric)[-1].split('_permanova.html')[0].split('__')
            elif fil.endswith('permdisp.html'):
                method = 'permdisp'
                subset_test = fil.split('_%s__' % metric)[-1].split('_permdisp.html')[0].split('__')
            subset = subset_test[0]
            test = '__'.join(subset_test[1:])
            lines = [x.strip() for x in open(path).readlines()]
            fil_data = []
            for ldx, line in enumerate(lines):
                if '<th>sample size</th>' in line:
                    size = lines[ldx+1].split('<td>')[-1].split('</td>')[0]
                elif '<th>number of groups</th>' in line:
                    number_of_groups = lines[ldx+1].split('<td>')[-1].split('</td>')[0]
                elif '<th>test statistic</th>' in line:
                    test_statistic = lines[ldx+1].split('<td>')[-1].split('</td>')[0]
                elif '<th>p-value</th>' in line:
                    p_value = lines[ldx+1].split('<td>')[-1].split('</td>')[0]
            cv_fp = '%s.cv' % splitext(path)[0]
            cv = {}
            if isfile(cv_fp):
                cv = dict(x.strip().split('\t') for x in open(cv_fp).readlines())
            perms.append([dataset, metric, subset, method, test, size,
                          number_of_groups, test_statistic, p_value, cv])
            if dataset == 'index.html':
                print(path)
perms_pd = pd.DataFrame(perms, columns=['dataset', 'metric', 'subset', 'method', 'test', 'size',
                                        'number_of_groups', 'test_statistic', 'p_value', 'cv'])
perms_pd = perms_pd.replace(
    {'test': dict(
        (x, x.replace('__vioscreen_component_percents', ''))
        for x in perms_pd.test.unique().tolist() if '__vioscreen_component_percents' in x)})
perms_pd = perms_pd.replace(
    {'test': dict(
        (x, x.replace('__vioscreen_micromacro', ''))
        for x in perms_pd.test.unique().tolist() if '__vioscreen_micromacro' in x)})


mats = {}
for dataset, perms_dat_pd in perms_pd.groupby('dataset'):
    pdf_fp = 'ROUTINE_FOLDER/permanova/DATASET/permdisp_matrices.pdf'
    with PdfPages(pdf_fp) as pdf:
        for subset, perms_filt_pd in perms_dat_pd.groupby('subset'):
            print(dataset, subset)

            perms_test_cols = ['method', 'metric', 'test', 'test_statistic']
            perms_filt_test = perms_filt_pd[perms_test_cols].sort_values(
                perms_test_cols[:-1]).set_index(perms_test_cols[:-1]).unstack()
            perms_filt_test.columns = perms_filt_test.columns.droplevel()
            perms_filt_test = perms_filt_test.reset_index()

            perms_pval_cols = ['method', 'metric', 'test', 'p_value']
            perms_filt_pval = perms_filt_pd[perms_pval_cols].copy()
            perms_filt_pval.p_value = pd.cut(
                x=perms_filt_pval['p_value'].astype(float),
                bins=[-1, 0.01, 0.03, 0.05, 1],
                labels=['***', '**', '*', '']
            )
            perms_filt_pval = perms_filt_pval.sort_values(
                perms_pval_cols[:-1]).set_index(perms_pval_cols[:-1]).unstack()
            perms_filt_pval.columns = perms_filt_pval.columns.droplevel()
            perms_filt_pval = perms_filt_pval.reset_index()

            perms_pval_cols = ['method', 'metric', 'test', 'size']
            perms_filt_size = perms_filt_pd[perms_pval_cols].copy()
            perms_filt_size = perms_filt_size.sort_values(
                perms_pval_cols[:-1]).set_index(perms_pval_cols[:-1]).unstack()
            perms_filt_size.columns = perms_filt_size.columns.droplevel()
            perms_filt_size = perms_filt_size.reset_index()

            perms_filt_cvs = dict(x for x in perms_filt_pd[['test', 'cv']].values)

            ntests = perms_filt_test['test'].unique().size

            fig, axes = plt.subplots(2, 1, figsize=((2 + ntests), 4))
            for mdx, method in enumerate(['permanova', 'permdisp']):
                method_perms_filt_test = perms_filt_test.loc[perms_filt_test.method == method].set_index(
                    'metric').drop(columns='method')
                method_perms_filt_pval = perms_filt_pval.loc[perms_filt_pval.method == method].set_index(
                    'metric').drop(columns='method')
                method_perms_filt_size = perms_filt_size.loc[perms_filt_size.method == method].set_index(
                    'metric').drop(columns='method')

                method_perms_filt_pval_size = (
                    method_perms_filt_pval.astype(str) + '\n' + method_perms_filt_size.astype(str)
                ).replace(
                    dict((x, {'\nnan': '', '***\nnan': '***', '**\nnan': '**', '*\nnan': '*'})
                         for x in method_perms_filt_size.columns)
                )
                g = sns.heatmap(
                    method_perms_filt_test.astype(float),
                    annot=method_perms_filt_pval,
                    fmt='', cmap="viridis",
                    annot_kws={"size": 10},
                    ax=axes[mdx]
                )
                g.set_title(method, loc='left')
                if not mdx:
                    g.set_xticklabels(['' for x in g.get_xticklabels()])
                    g.set_xlabel('')
                else:
                    g.set_xticklabels(['%s [%s]' % (
                        x.get_text(), ', '.join(
                            ['%s:%s' % (y, z) for y, z in perms_filt_cvs[x.get_text()].items()]
                        )
                    ) for x in g.get_xticklabels()])

            plt.suptitle('%s - %s (%s samples)' % (
                dataset, subset, method_perms_filt_size.stack().astype(int).max()))
            pdf.savefig(bbox_inches='tight')
            plt.close()
