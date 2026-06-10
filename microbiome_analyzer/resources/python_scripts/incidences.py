#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2026, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import argparse
import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True, help='Input feature table (.tsv)')
    parser.add_argument(
        '-o', nargs=1, required=True,help='Output incidence table')
    parser.add_argument(
        '-d', nargs=1, required=True, help='Dataset name')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_incidences(pbiom, mins, minp, mina):
    tab = pbiom.copy()
    relab = (tab / tab.sum()).fillna(0)
    tab = tab.where(relab >= mina, other=0)
    tab = tab.loc[tab.sum(1) > 0, tab.sum() > 0]
    pdata = round((tab > 0).sum(1).describe(), 3)
    pdata.index = ['features'] + [
        'features_prevalence_%s' % x for x in pdata.index[1:]]
    ndata = round(tab.sum().describe(), 3)
    ndata.index = ['samples'] + [
        'samples_reads_%s' % x for x in ndata.index[1:]]
    pres = ndata.to_dict()
    pres['min_sample_reads'] = mins
    pres.update(pdata)
    if minp >= 1:
        pres['min_n_prevalence'] = minp
    else:
        pres['min_%s_prevalence' % "%"] = 100 * minp
    pres['min_%s_abundance' % "%"] = 100 * mina
    return pres


def get_res(data, sams, smax):
    res = []
    min_sams = [0, 1000, 2000, 5000, 10000, 50000]
    min_prevs = [x / 1000 for x in range(1, 21)]
    min_abunds = [round(x / 1000, 3) for x in np.logspace(0, 3, 10)][:-1]
    for mins in min_sams:
        if mins > smax:
            continue
        sbiom = data.loc[:, sams >= mins]
        prevs = (sbiom > 0).sum(axis=1)
        for minp in sorted(prevs.unique())[::-1][1:]:
            pbiom = sbiom.loc[prevs >= minp, :]
            for mina in min_abunds:
                res.append(get_incidences(pbiom, mins, minp, mina))
        prevs = (sbiom > 0).sum(axis=1) / sbiom.shape[1]
        for minp in min_prevs:
            pbiom = sbiom.loc[prevs >= minp, :]
            for mina in min_abunds:
                res.append(get_incidences(pbiom, mins, minp, mina))
    return res


def incidences(filin, filou, dat):
    data = pd.read_table(filin, index_col=0)
    sams = data.sum()
    smax = sams.max()
    r = pd.DataFrame(get_res(data, sams, smax))
    r['features'] = r['features'].astype(int)
    r['samples'] = r['samples'].astype(int)
    r['features_prevalence_min'] = r['features_prevalence_min'].astype(int)
    r['features_prevalence_max'] = r['features_prevalence_max'].astype(int)
    r['samples_reads_min'] = r['samples_reads_min'].astype(int)
    r['samples_reads_max'] = r['samples_reads_max'].astype(int)
    r['dataset'] = dat
    ord = ['dataset', 'min_sample_reads', 'min_n_prevalence',
           'min_%s_prevalence' % "%", 'min_%s_abundance' % "%",
           'samples', 'features']
    r = r[ord + [x for x in r.columns if x not in ord]]
    r.to_csv(filou, index=False, sep='\t')


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    dat = args['d'][0]
    incidences(filin, filou, dat)
