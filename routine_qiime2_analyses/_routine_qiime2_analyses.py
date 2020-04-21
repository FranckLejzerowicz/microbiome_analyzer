# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import subprocess
from os.path import abspath, exists, isdir, isfile

from routine_qiime2_analyses._routine_q2_io_utils import get_prjct_nm, get_datasets
from routine_qiime2_analyses._routine_q2_filter import import_datasets, filter_rare_samples
from routine_qiime2_analyses._routine_q2_rarefy import run_rarefy
from routine_qiime2_analyses._routine_q2_phylo import shear_tree, run_sepp, get_precomputed_trees
from routine_qiime2_analyses._routine_q2_qemistree import run_qemistree
from routine_qiime2_analyses._routine_q2_taxonomy import run_taxonomy, run_barplot, get_precomputed_taxonomies
from routine_qiime2_analyses._routine_q2_alpha import (run_alpha, merge_meta_alpha, export_meta_alpha,
                                                       run_correlations, run_volatility,
                                                       run_alpha_group_significance)
from routine_qiime2_analyses._routine_q2_beta import run_beta, export_beta, run_pcoas, run_emperor
from routine_qiime2_analyses._routine_q2_deicode import run_deicode
from routine_qiime2_analyses._routine_q2_permanova import run_permanova
from routine_qiime2_analyses._routine_q2_adonis import run_adonis
from routine_qiime2_analyses._routine_q2_songbird import run_songbird
from routine_qiime2_analyses._routine_q2_mmvec import run_mmvec
from routine_qiime2_analyses._routine_q2_mmbird import run_mmbird


def routine_qiime2_analyses(
        i_datasets: tuple,
        i_datasets_folder: str,
        project_name: str,
        p_longi_column: str,
        thresh: int,
        p_alpha_subsets: str,
        p_perm_tests: tuple,
        p_perm_groups: str,
        p_formulas: str,
        force: bool,
        i_classifier: str,
        i_wol_tree: str,
        i_sepp_tree: str,
        i_qemistree: str,
        p_diff_models: str,
        p_mmvec_pairs: str,
        qiime_env: str,
        chmod: str,
        p_skip: tuple,
        gpu: bool,
        standalone: bool,
        raref: bool) -> None:
    """
    Main qiime2 functions writer.

    :param i_datasets: Internal name identifying the datasets in the input folder.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param project_name: Nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param p_alpha_subsets: Subsets for alpha diversity.
    :param p_perm_tests: Subsets for PERMANOVA.
    :param p_perm_groups: Groups to test between in each PERMANOVA subset (yml file path).
    :param p_formulas: Formula for Adonis tests for each PERMANOVA subset (yml file path).
    :param p_longi_column: If data is longitudinal; provide the time metadata column for volatility analysis.
    :param thresh: Minimum number of reads per sample to be kept.
    :param force: Force the re-writing of scripts for all commands.
    :param i_classifier: Path to the taxonomic classifier.
    :param i_wol_tree: default to ./routine_qiime2_analyses/resources/wol_tree.nwk.
    :param i_sepp_tree: path to the SEPP database artefact. Default to None.
    :param i_qemistree: path to the tree generated using Qemistree (for metabolomics datasets).
    :param p_diff_models: Formulas for multinomial regression-based differential abundance ranking.
    :param p_mmvec_pairs: Pairs of datasets for which to compute co-occurrences probabilities.
    :param chmod: whether to change permission of output files (defalt: 775).
    :param p_skip: steps to skip.
    :param gpu: Use GPUs instead of CPUs for MMVEC.
    :param standalone:
    :param raref: Whether to only perform the routine analyses on the rarefied datasets.
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

    # --> datasets_read <--
    # path_pd : indexed with feature name
    # meta_pd : not indexed -> "sample_name" as first column
    import_datasets(i_datasets_folder, datasets, datasets_phylo, force, prjct_nm, qiime_env, chmod)
    if raref:
        run_rarefy(i_datasets_folder, datasets, datasets_read,
                   datasets_features, datasets_phylo,
                   force, prjct_nm, qiime_env, chmod)
    if thresh:
        filter_rare_samples(i_datasets_folder, datasets, datasets_read, datasets_features,
                            datasets_phylo, prjct_nm, qiime_env, thresh, chmod)
    trees = {}
    get_precomputed_trees(i_datasets_folder, datasets, datasets_phylo, trees)
    if 'wol' not in p_skip:
        shear_tree(i_datasets_folder, datasets_read, datasets_phylo, datasets_features, prjct_nm,
                   i_wol_tree, trees, force, qiime_env, chmod)
    if i_sepp_tree and 'sepp' not in p_skip:
        run_sepp(i_datasets_folder, datasets, datasets_read, datasets_phylo, prjct_nm,
                 i_sepp_tree, trees, force, qiime_env, chmod)

    taxonomies = {}
    get_precomputed_taxonomies(i_datasets_folder, datasets, taxonomies)
    if i_qemistree and 'qemistree' not in p_skip:
        if isdir(i_qemistree):
            run_qemistree(i_datasets_folder, datasets, prjct_nm,
                          i_qemistree, taxonomies, force, qiime_env, chmod)
        else:
            print('[Warning] The Qemistree path %s is not a folder.')

    if 'taxonomy' not in p_skip:
        run_taxonomy(i_datasets_folder, datasets, datasets_read, datasets_phylo, datasets_features,
                     i_classifier, taxonomies, force, prjct_nm, qiime_env, chmod)

        if 'barplot' not in p_skip:
            run_barplot(i_datasets_folder, datasets, taxonomies, force, prjct_nm, qiime_env, chmod)
    # ------------------------------------------------------------------------------------------

    # ALPHA ------------------------------------------------------------
    if 'alpha' not in p_skip:
        diversities = run_alpha(i_datasets_folder, datasets, datasets_read,
                                datasets_phylo, p_alpha_subsets, trees,
                                force, prjct_nm, qiime_env, chmod)
        if 'merge_alpha' not in p_skip:
            to_export = merge_meta_alpha(i_datasets_folder, datasets, diversities,
                                         force, prjct_nm, qiime_env, chmod)
            if 'export_alpha' not in p_skip:
                export_meta_alpha(datasets, to_export)
        if 'alpha_correlations' not in p_skip:
            run_correlations(i_datasets_folder, datasets, diversities,
                             force, prjct_nm, qiime_env, chmod)
        if p_longi_column:
            if 'volatility' not in p_skip:
                run_volatility(i_datasets_folder, datasets, p_longi_column,
                               force, prjct_nm, qiime_env, chmod)
    # ------------------------------------------------------------------

    # BETA ----------------------------------------------------
    if 'beta' not in p_skip:
        betas = run_beta(i_datasets_folder, datasets, datasets_phylo,
                         trees, force, prjct_nm, qiime_env, chmod)
        if 'export_beta' not in p_skip:
            export_beta(i_datasets_folder, betas,
                        force, prjct_nm, qiime_env, chmod)
        if 'emperor' not in p_skip:
            pcoas = run_pcoas(i_datasets_folder, betas,
                              force, prjct_nm, qiime_env, chmod)
            run_emperor(i_datasets_folder, pcoas, prjct_nm, qiime_env, chmod)
    # ---------------------------------------------------------

    # STATS -----------------------------------------------------------------------
    if 'beta' not in p_skip and 'deicode' not in p_skip:
        run_deicode(i_datasets_folder, datasets, p_perm_groups,
                    force, prjct_nm, qiime_env, chmod)

    if 'alpha' not in p_skip and 'alpha_kw' not in p_skip:
        run_alpha_group_significance(i_datasets_folder, datasets, diversities,
                                     p_perm_groups, force, prjct_nm, qiime_env, chmod)
    if p_perm_tests:
        if 'beta' not in p_skip and 'permanova' not in p_skip:
            run_permanova(i_datasets_folder, datasets, betas,
                          p_perm_tests, p_perm_groups,force, prjct_nm, qiime_env, chmod)
    if p_formulas:
        if 'beta' not in p_skip and 'adonis' not in p_skip:
            run_adonis(p_formulas, i_datasets_folder, datasets, betas,
                       p_perm_groups, force, prjct_nm, qiime_env, chmod)
    # ------------------------------------------------------------------------------

    # MMVEC AND SONGBIRD -----------------------------------------------------------
    songbird_outputs, mmvec_outputs = {}, {}
    if p_diff_models:
        if 'songbird' not in p_skip:
            songbird_outputs = run_songbird(p_diff_models, i_datasets_folder, datasets,
                                            force, prjct_nm, qiime_env, chmod)
    if p_mmvec_pairs:
        if 'mmvec' not in p_skip:
            mmvec_outputs = run_mmvec(p_mmvec_pairs, i_datasets_folder, datasets,
                                      datasets_read, force, gpu, standalone,
                                      prjct_nm, qiime_env, chmod)
    if p_diff_models and p_mmvec_pairs:
        run_mmbird(i_datasets_folder, songbird_outputs, mmvec_outputs)
    # ------------------------------------------------------------------------------