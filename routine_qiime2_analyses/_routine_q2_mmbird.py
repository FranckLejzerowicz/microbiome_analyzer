# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, sys
import pandas as pd
from os.path import basename, isfile, splitext


def get_mmvec_outputs(mmvec_outputs: list):
    mmvec_outputs_pd = pd.DataFrame(
        mmvec_outputs,
        columns=[
             'pair',
             'omic1',
             'omic2',
             'filt1',
             'filt2',
             'n_common',
             'meta_common_fp',
             'omic1_common_fp',
             'omic2_common_fp',
             'omic1_common_qza',
             'omic2_common_qza',
             'mmvec_parameters',
             'mmvec_out'
         ])
    mmvec_outputs_pd = mmvec_outputs_pd.set_index(
        mmvec_outputs_pd.columns.tolist()[:-1]).unstack()
    mmvec_outputs_pd.columns = mmvec_outputs_pd.columns.droplevel()
    mmvec_outputs_pd.reset_index(inplace=True)
    return mmvec_outputs_pd


def get_songbird_outputs(songbird_outputs: list):
    songbird_outputs_pd = pd.DataFrame(
        songbird_outputs,
        columns=[
            'songbird_dat',
            'songbird_filt',
            'songbird_parameters',
            'songbird_case',
            'songbird_fp'
            'pair'
         ])
    # songbird_outputs_pd = mmvec_outputs_pd.set_index(
    #     mmvec_outputs_pd.columns.tolist()[:-1]).unstack()
    # mmvec_outputs_pd.columns = mmvec_outputs_pd.columns.droplevel()
    # mmvec_outputs_pd.reset_index(inplace=True)
    return songbird_outputs_pd


def run_mmbird(i_datasets_folder: str, datasets: dict, taxonomies: dict,
               songbird_outputs: list, mmvec_outputs: list,
               prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> None:

    # print("songbird_outputs")
    # print(songbird_outputs)

    if not mmvec_outputs:
        print('No mmvec output detected...')
        sys.exit(0)
    if not songbird_outputs:
        print('No songbird output detected...')
        sys.exit(0)

    mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
    mmvec_outputs_pd.to_csv('/Users/franck/Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/nut/mmvec_outputs_pd.tsv', index=False, sep='\t')

    songbird_outputs_pd = get_songbird_outputs(songbird_outputs)
    songbird_outputs_pd.to_csv('/Users/franck/Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/nut/songbird_outputs_pd.tsv', index=False, sep='\t')
    pass
