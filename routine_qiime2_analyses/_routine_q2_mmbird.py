# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# import os
# import itertools
import pandas as pd
from os.path import basename, isfile, splitext
# import multiprocessing


def get_mmvec_outputs(mmvec_outputs: list):
    mmvec_outputs_pd = pd.DataFrame(
        mmvec_outputs,
        columns=[
             'pair',
             'omic_filt1',
             'omic_filt2',
             'omic1',
             'omic2',
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
    mmvec_outputs_pd.to_csv('/Users/franck/Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/nut/mmvec_outputs_pd.tsv', index=False, sep='\t')
    return mmvec_outputs_pd


def run_mmbird(i_datasets_folder: str, datasets: dict, taxonomies: dict,
               songbird_outputs: dict, mmvec_outputs: dict,
               prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> None:

    # print("songbird_outputs")
    # print(songbird_outputs)

    songbird_table = []
    for omic, params in songbird_outputs.items():
        songbird_table.append(([omic] + params))
    songbird_table_pd = pd.DataFrame(
        songbird_table,
        columns=[
             'omic',
             'songbird_parameters',
             'songbird_fp',
             'songbird_q2'
         ])
    songbird_table_pd = songbird_table_pd.set_index(songbird_table_pd.columns.tolist()[:2]).unstack()
    songbird_table_pd.columns = songbird_table_pd.columns.droplevel()
    songbird_table_pd.reset_index(inplace=True)
    songbird_table_pd.to_csv('/Users/franck/Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/nut/songbird_table_pd.tsv', index=False, sep='\t')
    pass
