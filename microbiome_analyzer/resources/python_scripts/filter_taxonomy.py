#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2024, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import sys
import argparse
import pandas as pd
from qiime2 import Artifact


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True,
        help='Input taxonomy table in `.qza` format')
    parser.add_argument(
        '-o', nargs=1, required=True,
        help='Output taxonomy table in `.qza` format')
    parser.add_argument(
        '-t', nargs="*", required=True,
        help='Feature Table in `.qza` format')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def filter_taxonomy(taxin, taxou, feats):
    tab_pd = Artifact.load(feats).view(pd.DataFrame)
    tax_pd = Artifact.load(taxin).view(pd.DataFrame)
    tab_feats = set(tab_pd.columns)
    tax_feats = set(tax_pd.index)
    only_tab = tab_feats.difference(tax_feats)
    if only_tab:
        print("%s features only present in feature table:" % len(only_tab))
        if len(only_tab) > 10:
            print('Showing 10 first: %s ...' % ', '.join(sorted(only_tab)[:10]))
        else:
            print(', '.join(sorted(only_tab)))
        sys.exit(1)
    tax_pd = tax_pd.loc[list(tab_feats), :]
    tax_tmp = '%s.temp.qza' % taxou
    tax_pd.to_csv(tax_tmp, sep='\t')
    tax_qza = Artifact.import_data('FeatureData[Taxonomy]', tax_tmp)
    tax_qza.save(taxou)
    os.remove(tax_tmp)


if __name__ == '__main__':
    args = get_args()
    taxin = args['i'][0]
    taxou = args['o'][0]
    feats = args['t'][0]
    filter_taxonomy(taxin, taxou, feats)
