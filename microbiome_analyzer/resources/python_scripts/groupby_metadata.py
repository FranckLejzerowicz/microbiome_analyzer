#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2024, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True,
        help='Input metadata (tab-separated)')
    parser.add_argument(
        '-o', nargs=1, required=True,
        help='Output metadata (tab-separated)')
    parser.add_argument(
        '-c', nargs="*", required=True,
        help='Column(s) by which to group metadata samples into new samples')
    parser.add_argument(
        '-g', choices=['first', 'mode'], required=False, default='first',
        help='Group values are "first" of sorted samples, or "mode" of values)')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def add_groups(gb, cols):
    for column in cols:
        gb[column] = gb.index.get_level_values(column)


def reformat_meta(filin, cols, group):
    tab = pd.read_table(filin)
    tab = tab.sort_values(tab.columns.tolist()[0])
    tab = tab.iloc[:, 1:]
    if group == 'mode':
        gb = tab.groupby(cols).agg(lambda x: pd.Series.mode(x)[0])
        add_groups(gb, cols)
    else:
        gb = tab.groupby(cols).agg({x: 'first' for x in tab.columns})
    if len(cols) > 1:
        gb.index = ['_'.join(map(str, x)) for x in gb.index]
    return gb


def write_meta(gb, filou):
    gb.index.name = 'sample_name'
    gb.reset_index(inplace=True)
    gb.to_csv(filou, index=False, sep='\t')


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    cols = args['c']
    group = args['g'][0]
    gb = reformat_meta(filin, cols, group)
    write_meta(gb, filou)
