# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import itertools
import pandas as pd
from os.path import basename, isfile, splitext
import multiprocessing

# def run_mmbird(i_datasets_folder: str, songbird_outputs: dict, mmvec_outputs: dict) -> None:
#     pass