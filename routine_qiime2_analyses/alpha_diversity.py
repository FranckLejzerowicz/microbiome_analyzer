# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import pandas as pd
from os.path import isfile, splitext
from skbio.tree import TreeNode
from routine_qiime2_analyses.analyses_prep import AnalysisPrep
from routine_qiime2_analyses._routine_q2_cmds import (
    run_import, run_export, write_rarefy, write_seqs_fasta,
    write_fragment_insertion
)
from routine_qiime2_analyses._routine_q2_filter import (
    get_thresholds, no_filtering, get_dat_filt, filtering_thresholds,
    harsh_filtering
)
from routine_qiime2_analyses._routine_q2_rarefy import (
    get_digit_depth,
)
from routine_qiime2_analyses._routine_q2_io_utils import (
    read_yaml_file
)
from routine_qiime2_analyses._routine_q2_taxonomy import (
    get_taxonomy_command, get_edit_taxonomy_command, get_gid_features,
    get_split_levels, fix_collapsed_data, write_collapse_taxo

)
from routine_qiime2_analyses._routine_q2_phylo import (
    get_wol_tree, get_sepp_tree
)
from routine_qiime2_analyses._routine_q2_xpbs import (
    print_message
)


class AlphaDiversity(object):
    """
    """

    def __init__(self, config, project) -> None:
        self.config = config
        self.project = project
        self.analysis = ''
        self.cmds = {}

    def register_command(self):
        AnalysisPrep.analyses_commands[self.analysis] = self.cmds[self.analysis]
