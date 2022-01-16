# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from os.path import isfile

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

fig_o = '<FIG_O>'
dat = '<DATASET>'
dat_fp = '<DAT_FP>'
dat_pd = pd.read_csv(dat_fp, sep='\t')

rarefs = set()
cohorts = set()
modes = set()
tabs = []
for raref_, raref_pd in dat_pd.groupby('raref'):
    raref = 'raw'
    if raref_:
        raref = raref_.split('_raref')[-1]
    for (me, co), meco_pd in raref_pd.groupby(['metric', 'cohort']):
        cohorts.add(co)
        for (mode, mode_group, fp) in meco_pd[['mode', 'mode_group',
                                               'tsv']].values:
            if not isfile(fp):
                continue
            cur_mode = mode
            if mode_group:
                cur_mode = '%s (%s)' % (mode, mode_group)
            modes.add(cur_mode)
            decay = pd.read_csv(fp)
            decay['rarefaction'] = raref
            decay['metric'] = me
            decay['metric / rarefaction'] = decay.apply(
                lambda x: ' / '.join(x[['metric', 'rarefaction']]),
                axis=1)
            decay['samples subset'] = co
            decay['analysis mode'] = cur_mode
            tabs.append(decay)

tab_pd = pd.concat(tabs)
if 'aitchison' in set(tab_pd.metric):
    pds = {'_aitchison': tab_pd.loc[tab_pd.metric == 'aitchison'],
           '': tab_pd.loc[tab_pd.metric != 'aitchison']}
else:
    pds = {'': tab_pd.copy()}

with PdfPages(fig_o) as pdf:
    for aitchison, tab_pd in pds.items():
        col = 'samples subset'
        if cohorts == {'ALL'}:
            col = ''
        title = dat
        style = 'analysis mode'
        hue = 'metric / rarefaction'
        if rarefs == {'raw'}:
            title += ' - no rarefaction'
            hue = 'metric'
        if modes == {'analysis mode'}:
            if col:
                style = col
                col = ''
        height = 4
        if col:
            g = sns.relplot(
                data=tab_pd, x='step', y='min', hue=hue, style=style,
                col=col, kind='line', col_wrap=3, facet_kws={
                    'sharex': False}, height=height)
            g.set_titles('{col_name}')
        else:
            g = sns.relplot(
                data=tab_pd, x='step', y='min', hue=hue, style=style,
                kind='line', height=height)
        plt.suptitle(title, fontsize=12)
        plt.subplots_adjust(top=0.93)
        pdf.savefig(bbox_inches='tight')
print('[dm_decay] Written figure: %s' % fig_o)
