# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys


def check_sample_subset_format(yaml):
    """Checks whether the sample subset yaml file given
    to option -g, --sample-subsets is in the correct format.

    Parameters
    ----------
    yaml : dict
        Input parameters in loaded yaml format
    """
    opt = "[-g, --sample-subsets] YAML content must be a dict (keys=names)"
    if not isinstance(yaml, dict):
        sys.exit(opt)
    opt += ', containing nested dicts (keys=metadata variables)'
    if [1 for x in yaml.values() if not isinstance(x, dict)]:
        sys.exit(opt)
    opt += ', which values are lists (possible factors of the variables)'
    for name, val_vals in yaml.items():
        for val, vals in val_vals.items():
            if not isinstance(vals, list):
                sys.exit(opt)


def check_format(yaml, input_data: str):
    """Checks whether the yaml file used as input is in the correct format.

    Parameters
    ----------
    yaml : dict
        Input parameters in loaded yaml format
    input_data : str
        String telling what is the expected import format
    """
    if input_data == 'sample_subsets':
        check_sample_subset_format(yaml)


