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
from routine_qiime2_analyses._routine_q2_phylo import shear_tree
from routine_qiime2_analyses._routine_q2_filter import filter_rare_samples
from routine_qiime2_analyses._routine_q2_beta import run_beta, export_beta, run_pcoas, run_emperor
from routine_qiime2_analyses._routine_q2_alpha import (run_alpha, merge_meta_alpha, export_meta_alpha,
                                                       run_correlations, run_volatility,
                                                       run_alpha_group_significance)
from routine_qiime2_analyses._routine_q2_permanova import run_permanova
from routine_qiime2_analyses._routine_q2_deicode import run_deicode


def routine_qiime2_analyses(i_datasets: tuple, i_folder: str, project_name: str,
                            gid: bool, p_longi_column: str, thresh: int,
                            p_perm_subsets: str, p_perm_groups: str,
                            force: bool, i_wol_tree: str, qiime_env: str, biom: bool):
    """
    Main qiime2 functions writer.

    :param i_datasets: Internal name identifying the datasets in the input folder.
    :param i_folder: Path to the folder containing the .tsv datasets.
    :param i_wol_tree: default on barnacle /projects/wol/profiling/dbs/wol/phylogeny/web_of_life_tree.nwk.
    :param project_name: Nick name for your project.
    :param p_qiime2_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param p_perm_subsets: Subsets for PERMANOVA.
    :param p_perm_groups: Groups to tests between in each PERMANOVA subset (yml file path).
    :param p_longi_column: If data is longitudinal; provide the time metadata column for volatility analysis.
    :param thresh: Minimum number of reads per sample to be kept.
    :param force: Force the re-writing of scripts for all commands.
    :param gid: If feature names have the genome ID (to use the Web of Life tree).
    :param biom: Use biom files in the input folder
    :return None
    """

    # check input
    if not exists(i_folder):
        print('%s is not an existing folder\nExiting...' % i_folder)
        sys.exit(1)

    i_folder = abspath(i_folder)
    if isfile(i_folder):
        print('%s is a file. Needs a folder as input\nExiting...' % i_folder)
        sys.exit(1)

    # check Xpbs
    ret_code, ret_path = subprocess.getstatusoutput('which Xpbs')
    if ret_code:
        print('Xpbs is not installed (and ake sure to edit its config.txt)\nExiting...')
        sys.exit(1)

    prjct_nm = get_prjct_nm(project_name)

    datasets, datasets_read, datasets_features = get_datasets(i_datasets, i_folder, gid, biom)

    if thresh:
        filter_rare_samples(i_folder, datasets, datasets_read, datasets_features,
                            force, prjct_nm, qiime_env, thresh, gid)

    distances = {'alpha': ['observed_otus', 'pielou_e', 'shannon'],
                 'beta': ['jaccard', 'braycurtis', 'aitchison']}
    wol_trees = {}
    if datasets_features:
        distances['alpha'].append('faith_pd')
        distances['beta'].extend(['unweighted_unifrac', 'weighted_unifrac'])
        wol_trees = shear_tree(datasets_features, i_folder, prjct_nm,
                               i_wol_tree, force, qiime_env)

    betas = run_beta(i_folder, datasets, distances['beta'],
                     wol_trees, force, prjct_nm, qiime_env)
    export_beta(i_folder, betas,
                force, prjct_nm, qiime_env)
    pcoas = run_pcoas(i_folder, betas,
                      force, prjct_nm, qiime_env)
    run_emperor(i_folder, pcoas, prjct_nm, qiime_env)

    diversities = run_alpha(i_folder, datasets, distances['alpha'],
                            wol_trees, force, prjct_nm, qiime_env)
    to_export = merge_meta_alpha(i_folder, diversities,
                                 force, prjct_nm, qiime_env)
    export_meta_alpha(i_folder, to_export,
                      prjct_nm, qiime_env)
    run_correlations(i_folder, datasets, diversities,
                     force, prjct_nm, qiime_env)
    if p_longi_column:
        run_volatility(i_folder, datasets, p_longi_column,
                       force, prjct_nm, qiime_env)

    if p_perm_groups:
        run_deicode(i_folder, datasets, p_perm_groups,
                    force, prjct_nm, qiime_env)
        run_alpha_group_significance(i_folder, diversities, p_perm_groups,
                                     distances['alpha'], force, prjct_nm, qiime_env)

        if p_perm_subsets:
            run_permanova(i_folder, datasets, betas,
                          distances['beta'], p_perm_subsets, p_perm_groups,
                          force, prjct_nm, qiime_env)
