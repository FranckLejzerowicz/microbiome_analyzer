# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# import os
# import itertools
# import pandas as pd
# from os.path import basename, isfile, splitext
# import multiprocessing


def run_mmbird(i_datasets_folder: str, datasets: dict, taxonomies: dict,
               songbird_outputs: dict, mmvec_outputs: dict,
               prjct_nm: str, qiime_env: str, chmod: str) -> None:
    print(i_datasets_folder)
    print(songbird_outputs)
    print(mmvec_outputs)

