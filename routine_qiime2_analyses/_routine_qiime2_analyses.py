# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import subprocess
from os.path import abspath, isfile, exists

from routine_qiime2_analyses._routine_q2_io_utils import get_prjct_nm, get_datasets
from routine_qiime2_analyses._routine_q2_phylo import shear_tree, run_sepp
from routine_qiime2_analyses._routine_q2_filter import import_datasets, filter_rare_samples
from routine_qiime2_analyses._routine_q2_beta import run_beta, export_beta, run_pcoas, run_emperor
from routine_qiime2_analyses._routine_q2_alpha import (run_alpha, merge_meta_alpha, export_meta_alpha,
                                                       run_correlations, run_volatility,
                                                       run_alpha_group_significance)
from routine_qiime2_analyses._routine_q2_deicode import run_deicode
from routine_qiime2_analyses._routine_q2_permanova import run_permanova
from routine_qiime2_analyses._routine_q2_adonis import run_adonis


def routine_qiime2_analyses(i_datasets: tuple, i_datasets_folder: str, project_name: str,
                            p_longi_column: str, thresh: int, p_perm_tests: tuple,
                            p_perm_groups: str, p_formulas: str, force: bool,
                            i_wol_tree: str, i_sepp_tree: str, qiime_env: str) -> None:
    """
    Main qiime2 functions writer.

    :param i_datasets: Internal name identifying the datasets in the input folder.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param project_name: Nick name for your project.
    :param p_qiime2_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param p_perm_tests: Subsets for PERMANOVA.
    :param p_perm_groups: Groups to test between in each PERMANOVA subset (yml file path).
    :param p_formulas: Formula for Adonis tests for each PERMANOVA subset (yml file path).
    :param p_longi_column: If data is longitudinal; provide the time metadata column for volatility analysis.
    :param thresh: Minimum number of reads per sample to be kept.
    :param force: Force the re-writing of scripts for all commands.
    :param i_wol_tree: default to ./routine_qiime2_analyses/resources/wol_tree.nwk.
    :param i_sepp_tree: path to the SEPP database artefact. Default to None.
    """

    # check input
    if not exists(i_datasets_folder):
        print('%s is not an existing folder\nExiting...' % i_datasets_folder)
        sys.exit(1)

    i_datasets_folder = abspath(i_datasets_folder)
    if isfile(i_datasets_folder):
        print('%s is a file. Needs a folder as input\nExiting...' % i_datasets_folder)
        sys.exit(1)

    # check Xpbs
    ret_code, ret_path = subprocess.getstatusoutput('which Xpbs')
    if ret_code:
        print('Xpbs is not installed (and make sure to edit its config.txt)\nExiting...')
        sys.exit(1)
    else:
        with open(ret_path) as f:
            for line in f:
                break
        if line.startswith('$HOME'):
            print('Xpbs is installed but its config.txt need editing!\nExiting...')
            sys.exit(1)

    prjct_nm = get_prjct_nm(project_name)

    # INIT -------------------------------------------------------------------------------------
    datasets, datasets_read, datasets_features, datasets_phylo = get_datasets(i_datasets,
                                                                              i_datasets_folder)
    import_datasets(i_datasets_folder, datasets, datasets_phylo, force, prjct_nm, qiime_env)
    if thresh:
        filter_rare_samples(i_datasets_folder, datasets, datasets_read, datasets_features,
                            datasets_phylo, force, prjct_nm, qiime_env, thresh)
    trees = {}
    shear_tree(i_datasets_folder, datasets_phylo, datasets_features, prjct_nm,
               i_wol_tree, trees, force, qiime_env)
    run_sepp(i_datasets_folder, datasets, datasets_read, datasets_phylo, prjct_nm,
             i_sepp_tree, trees, force, qiime_env)
    # ------------------------------------------------------------------------------------------

    # ALPHA ------------------------------------------------------------
    diversities = run_alpha(i_datasets_folder, datasets, datasets_phylo,
                            trees, force, prjct_nm, qiime_env)
    to_export = merge_meta_alpha(i_datasets_folder, diversities,
                                 force, prjct_nm, qiime_env)
    export_meta_alpha(i_datasets_folder, to_export,
                      force, prjct_nm, qiime_env)
    run_correlations(i_datasets_folder, datasets, diversities,
                     force, prjct_nm, qiime_env)
    if p_longi_column:
        run_volatility(i_datasets_folder, datasets, p_longi_column,
                       force, prjct_nm, qiime_env)
    # ------------------------------------------------------------------

    # BETA ----------------------------------------------------
    betas = run_beta(i_datasets_folder, datasets, datasets_phylo,
                     trees, force, prjct_nm, qiime_env)
    export_beta(i_datasets_folder, betas,
                force, prjct_nm, qiime_env)
    pcoas = run_pcoas(i_datasets_folder, betas,
                      force, prjct_nm, qiime_env)
    run_emperor(i_datasets_folder, pcoas, prjct_nm, qiime_env)
    # ---------------------------------------------------------

    # STATS -----------------------------------------------------------------------
    run_deicode(i_datasets_folder, datasets, p_perm_groups,
                force, prjct_nm, qiime_env)
    run_alpha_group_significance(i_datasets_folder, diversities, p_perm_groups,
                                 force, prjct_nm, qiime_env)

    if p_perm_tests:
        run_permanova(i_datasets_folder, datasets, betas,
                      p_perm_tests, p_perm_groups,
                      force, prjct_nm, qiime_env)
    if p_formulas:
        run_adonis(p_formulas, i_datasets_folder, datasets, betas,
                   p_perm_groups, force, prjct_nm, qiime_env)
    # ------------------------------------------------------------------------------
