# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import yaml
import pandas as pd

from os.path import isfile


def read_yaml_file(file_path: str) -> dict:
    """Simply reads a yaml and return
    its contents in a dictionary structure.

    Parameters
    ----------
    file_path: str
        Path to a yaml file.

    Returns
    -------
    yaml_dict : dict
        Dictionary object returned
        by reading the yaml file.
        (could be an empty dict)
    """
    yaml_dict = {}
    if file_path and isfile(file_path):
        with open(file_path) as yaml_handle:
            try:
                yaml_dict = yaml.load(
                    yaml_handle, Loader=yaml.FullLoader)
            except AttributeError:
                yaml_dict = yaml.load(yaml_handle)
    return yaml_dict


def read_meta_pd(
        meta_tab: str,
        rep_col: str = 'sample_name') -> pd.DataFrame:
    """Read metadata with first column as index.

    Parameters
    ----------
    meta_tab : str
        Path to the data table file
    rep_col : str
        Columns to use for index name

    Returns
    -------
    meta_tab_pd : pd.DataFrame
        Metadata table
    """
    meta_tab_sam_col = get_feature_sample_col(meta_tab)
    meta_tab_pd = pd.read_csv(
        meta_tab,
        header=0,
        sep='\t',
        dtype={meta_tab_sam_col: str},
        low_memory=False
    )
    meta_tab_pd.rename(
        columns={meta_tab_sam_col: rep_col},
        inplace=True
    )
    return meta_tab_pd


def get_feature_sample_col(meta_tab: str) -> str:
    """Get the first column of the metadata file.

    Parameters
    ----------
    meta_tab : str
        Data of metadata file path

    Returns
    -------
    first_column : str
        name of the first column in the table
    """
    n = 0
    with open(meta_tab) as f:
        for line in f:
            n += 1
            break
    if n:
        first_column = line.split('\t')[0]
        return first_column
    else:
        sys.exit('Empty: %s (possibly being written elsewhere..)' % meta_tab)


