# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from os.path import isdir

from routine_qiime2_analyses.dataset_collection import Datasets
from routine_qiime2_analyses.create_scripts import CreateScripts
from routine_qiime2_analyses.analyses_config import AnalysesConfig
from routine_qiime2_analyses.analyses_prep import AnalysisPrep
from routine_qiime2_analyses.alpha_diversity import AlphaDiversity
from routine_qiime2_analyses.paired_data import PairedData
from routine_qiime2_analyses.differential_abundance import DiffModels
from routine_qiime2_analyses.post_analyses import PostAnalysis

from routine_qiime2_analyses._routine_q2_io_utils import (
    check_input, get_filt_raref_suffix, check_xpbs_install, get_conda_envs,
    get_prjct_anlss_nm, get_datasets, get_run_params, summarize_songbirds,
    get_analysis_folder, get_train_test_dict)
from routine_qiime2_analyses._routine_q2_filter import (
    import_datasets, filter_rare_samples, get_filt3d_params,
    explore_filtering,
    deleted_non_filt)
from routine_qiime2_analyses._routine_q2_rarefy import run_rarefy
from routine_qiime2_analyses._routine_q2_phylo import (
    shear_tree, run_sepp, get_precomputed_trees)
from routine_qiime2_analyses._routine_q2_qemistree import run_qemistree
from routine_qiime2_analyses._routine_q2_taxonomy import (
    run_taxonomy, run_barplot, run_collapse, make_pies, get_taxo_levels,
    edit_taxonomies, get_precomputed_taxonomies,
    create_songbird_feature_metadata)
from routine_qiime2_analyses._routine_q2_doc import run_doc
from routine_qiime2_analyses._routine_q2_sourcetracking import (
    run_sourcetracking)
from routine_qiime2_analyses._routine_q2_alpha import (
    run_alpha, merge_meta_alpha, export_meta_alpha, run_correlations,
    run_alpha_rarefaction, run_volatility, run_alpha_group_significance)
from routine_qiime2_analyses._routine_q2_beta import (
    run_beta, export_beta, run_pcoas, run_biplots, run_emperor,
    run_emperor_biplot, run_empress, run_empress_biplot)
from routine_qiime2_analyses._routine_q2_decay import (
    run_distance_decay, distance_decay_figure)
from routine_qiime2_analyses._routine_q2_procrustes_mantel import (
    run_procrustes, run_mantel)
from routine_qiime2_analyses._routine_q2_deicode import run_deicode
from routine_qiime2_analyses._routine_q2_permanova import (
    run_permanova, summarize_permanova)
from routine_qiime2_analyses._routine_q2_nestedness import (
    run_nestedness, nestedness_graphs, nestedness_nodfs)
from routine_qiime2_analyses._routine_q2_adonis import run_adonis
from routine_qiime2_analyses._routine_q2_phate import run_phate
from routine_qiime2_analyses._routine_q2_songbird import run_songbird
from routine_qiime2_analyses._routine_q2_mmvec import run_mmvec
from routine_qiime2_analyses._routine_q2_mmbird import run_mmbird


def routine_qiime2_analyses(
        i_datasets: tuple,
        i_datasets_folder: str,
        project_name: str,
        qiime_env: str,
        p_longi_column: str,
        p_filt_threshs: str,
        p_raref_depths: str,
        eval_rarefs: bool,
        p_alpha_subsets: str,
        p_beta_subsets: str,
        p_perm_tests: tuple,
        p_perm_tests_min: int,
        p_beta_groups: str,
        p_nestedness_groups: str,
        p_beta_type: tuple,
        p_procrustes: str,
        p_mantel: str,
        p_distance_decay: str,
        p_collapse_taxo: str,
        p_train_test: str,
        p_adonis_formulas: str,
        p_doc_config: str,
        p_sourcetracking_config: str,
        p_phate_config: str,
        do_biplots: bool,
        force: bool,
        i_classifier: str,
        i_wol_tree: str,
        i_sepp_tree: str,
        i_qemistree: str,
        p_diff_models: str,
        p_mmvec_pairs: str,
        p_mmvec_highlights: str,
        p_xmmvec: str,
        p_run_params: str,
        p_chmod: str,
        p_skip: tuple,
        gpu: bool,
        standalone: bool,
        raref: bool,
        noloc: bool,
        As: tuple,
        Bs: tuple,
        split: bool,
        dropout: bool,
        doc_phate: bool,
        filt3d: bool,
        p_filt3d_config: str,
        filt_only: bool,
        jobs: bool,
        chunkit: int) -> None:
# def routine_qiime2_analyses(*args, **kwargs):
    """Main qiime2 functions writer.

    Parameters
    ----------
        i_datasets : tuple
            Names identifying the datasets in the input folder
        i_datasets_folder : str
            Path to the folder containing the data/metadata sub-folders
        project_name : str
            Nick name for the project.
        qiime_env : str
            Name of a qiime2 conda environment where analysis
            tools to be run are installed
        p_longi_column : str
            Time metadata column if data is longitudinal
        p_filt_threshs : str
            Minimum sample read abundance to be kept in the sample
        p_raref_depths : str
        eval_rarefs : bool
        p_alpha_subsets : str
            Features subsets config for alpha diversity
        p_beta_subsets : str
        p_perm_tests : tuple
            Subsets for PERMANOVA
        p_perm_tests_min : int
        p_beta_groups : str
            Config defining the samples groups to test between in each
            PERMANOVA subset
        p_nestedness_groups : str
        p_beta_type : tuple
        p_procrustes : str
        p_mantel : str
        p_distance_decay : str
        p_collapse_taxo : str
        p_train_test : str
        p_formulas : str
            Config defining the formula for Adonis tests
            for each PERMANOVA subset
        p_doc_config : str
        p_sourcetracking_config : str
        p_phate_config : str
        do_biplots : bool
        force : bool
            Force the re-writing of scripts for all commands
        i_classifier : str
            Path to the taxonomic classifier
        i_wol_tree : str
            Path to the Web of Life (WOL) tree for Woltka data or any dataset
            relying on WOL genomes (default to
            ./routine_qiime2_analyses/resources/wol_tree.nwk)
        i_sepp_tree : str
            Path to the SEPP database QIIME2 artifact (no
            default: no read placement performed by default)
        i_qemistree : str
            Path to the tree generated using Qemistree (for metabolomics
            datasets)
        p_diff_models : str
            Config defining the formulas for multinomial regression-based
            differential abundance ranking (using songbird)
        p_mmvec_pairs : str
            Config defining the pairs of datasets for which to compute
            co-occurrences probabilities (using mmvec)
        p_mmvec_highlights : str
        p_xmmvec : str
        p_run_params : str
        p_chmod : str
            Whether to change permission of output files (default: 744)
        p_skip : tuple
            Analysis steps to skip
        gpu : bool
            Whether to use GPUs instead of CPUs for mmvec
        standalone : bool
            Whether to run songbird and mmvec as standalone
            (in order to access live, tensorboard visualizer)
        raref : bool
            Whether to repeat all analyses on rarefied datasets (must be
            set for a config passed to `-r` to be used)
        noloc : bool
        As : tuple
        Bs : tuple
        split : bool
        dropout : bool
        doc_phate : bool
        filt3d : bool
            Whether to build prevalence-abundance-features 3D plots for each
            dataset (based on prevalence-abundance filters passed in configs)
            instead of running any other analysis.
        p_filt3d_config : str
            Prevalence-abundance filters
        filt_only : bool
        jobs : bool
        chunkit : int
    """

    # INITIALIZATION ----------------------------------------------------------
    # check input
    # i_datasets, i_datasets_folder, project_name, qiime_env = args
    # config = AnalysesConfig(args, kwargs)
    config = AnalysesConfig(
        i_datasets, i_datasets_folder, project_name, qiime_env, raref,
        p_filt_threshs, p_longi_column, p_raref_depths, eval_rarefs,
        p_alpha_subsets, p_beta_subsets, p_perm_tests, p_perm_tests_min,
        p_beta_groups, p_nestedness_groups, p_beta_type, p_procrustes,
        p_mantel, p_distance_decay, p_collapse_taxo, p_train_test,
        p_adonis_formulas, p_doc_config, p_sourcetracking_config,
        p_phate_config, do_biplots, force, i_classifier, i_wol_tree,
        i_sepp_tree, i_qemistree, p_diff_models, p_mmvec_pairs,
        p_mmvec_highlights, p_xmmvec, p_run_params, p_chmod, p_skip, gpu,
        standalone, noloc, As, Bs, split, dropout, doc_phate, filt3d,
        p_filt3d_config, filt_only, jobs, chunkit)

    i_datasets_folder = check_input(i_datasets_folder)
    # Initialization checks
    if jobs:
        # Xpbs must be installed (to make .pbs scripts) if required
        check_xpbs_install()
    # get jobs' project name for display in torque scheduler
    # prjct_nm = get_prjct_nm(project_name)
    prjct_nm = config.get_prjct_anlss_nm()

    # get the names of the conda environments
    conda_envs = get_conda_envs(qiime_env)
    # get default, per-analysis run parameters (i.e. memory, #cpus, etc)
    # run_params = get_run_params(p_run_params, conda_envs)
    run_params = config.get_run_params()

    # READ --------------------------------------------------------------------
    print('(Get datasets)')
    datasets_objects = get_datasets(i_datasets, i_datasets_folder)
    datasets = datasets_objects[0]
    datasets_read = datasets_objects[1]
    datasets_features = datasets_objects[2]
    datasets_phylo = datasets_objects[3]
    datasets_rarefs = datasets_objects[4]

    project = Datasets(i_datasets, i_datasets_folder)

    # get the "", "_flt", "_rrf" or "_flt_rrf" suffix for the sh/pbs scripts
    # filt_raref = get_filt_raref_suffix(p_filt_threshs, raref)
    filt_raref = config.get_filt_raref_suffix()
    scripting = CreateScripts(config, prjct_nm, run_params, filt_raref)

    # this creates the "need" to do procrustes in "eval_rarefs" mode
    if eval_rarefs and not p_procrustes:
        p_procrustes = 1

    # workflow = config.get_workflow()

    # PREPROCESSING ------------------------------------------------------------
    print('(import_datasets)')
    import_datasets(i_datasets_folder, datasets, datasets_phylo,
                    force, prjct_nm, qiime_env, p_chmod, noloc,
                    run_params['import'], filt_raref, jobs, chunkit)
    AnalysisPrep('import').import_datasets(config, project)

    datasets_filt = {}
    datasets_filt_map = {}
    if p_filt_threshs:
        print('(filter_rare_samples)')
        filter_rare_samples(
            i_datasets_folder, datasets, datasets_read, datasets_features,
            datasets_rarefs, datasets_filt, datasets_filt_map, datasets_phylo,
            prjct_nm, qiime_env, p_filt_threshs, p_chmod, noloc,
            run_params['filter'], filt_raref, jobs, chunkit)
        AnalysisPrep('filter').filter_rare_samples(config, project)

    eval_depths = {}
    if raref:
        print('(run_rarefy)')
        eval_depths = run_rarefy(
            i_datasets_folder, datasets, datasets_read, datasets_phylo,
            datasets_filt_map, datasets_rarefs, p_raref_depths, eval_rarefs,
            force, prjct_nm, qiime_env, p_chmod, noloc, run_params['rarefy'],
            filt_raref, filt_only, jobs, chunkit)
        AnalysisPrep('rarefy').rarefy(config, project)

    # TAXONOMY ------------------------------------------------------------
    taxonomies = {}
    method = 'sklearn'
    print('(get_precomputed_taxonomies)')
    get_precomputed_taxonomies(
        i_datasets_folder, datasets, datasets_filt_map, taxonomies, method)
    project.get_precomputed_taxonomy(config)

    if i_qemistree and 'qemistree' not in p_skip:
        if isdir(i_qemistree):
            print('(run_qemistree)')
            run_qemistree(
                i_datasets_folder, datasets, prjct_nm, i_qemistree,
                taxonomies, force, qiime_env, p_chmod, noloc,
                run_params['qemistree'], filt_raref, jobs, chunkit)
        else:
            print('[Warning] The Qemistree path %s is not a folder.')

    if 'taxonomy' not in p_skip:
        print('(run_taxonomy)')
        run_taxonomy(
            method, i_datasets_folder, datasets, datasets_read,
            datasets_phylo, datasets_features, datasets_filt_map,
            i_classifier, taxonomies, force, prjct_nm, qiime_env, p_chmod,
            noloc, run_params['taxonomy'], filt_raref, jobs, chunkit)
        AnalysisPrep('taxonomy').taxonomy(config, project)

        print('(run_edit_taxonomies)')
        edit_taxonomies(
            i_datasets_folder, taxonomies, force, prjct_nm, qiime_env, p_chmod,
            noloc, run_params['taxonomy'], filt_raref, jobs, chunkit)

        if 'barplot' not in p_skip:
            print('(run_barplot)')
            run_barplot(
                i_datasets_folder, datasets, taxonomies, force, prjct_nm,
                qiime_env, p_chmod, noloc, run_params['barplot'], filt_raref,
                jobs, chunkit)

    # TREES ------------------------------------------------------------
    trees = {}
    print('(get_precomputed_trees)')
    get_precomputed_trees(
        i_datasets_folder, datasets, datasets_filt_map, datasets_phylo, trees)
    project.get_precomputed_trees(config)
    print()
    print()
    print('datasets_phylo')
    print(datasets_phylo)
    print()
    print()
    print('datasets_filt_map')
    print(datasets_filt_map)
    print()
    print()
    print('trees')
    print(trees)

    if 'wol' not in p_skip:
        print('(shear_tree)')
        print('fff')
        shear_tree(
            i_datasets_folder, datasets, datasets_read, datasets_phylo,
            datasets_features, prjct_nm, i_wol_tree, trees, datasets_rarefs,
            force, qiime_env, p_chmod, noloc, run_params['wol'], filt_raref,
            jobs)
        AnalysisPrep('wol').shear_tree(config, project)

    if i_sepp_tree and 'sepp' not in p_skip:
        print('(run_sepp)')
        run_sepp(
            i_datasets_folder, datasets, datasets_read, datasets_phylo,
            datasets_rarefs, prjct_nm, i_sepp_tree, trees, force, qiime_env,
            p_chmod, noloc, run_params['sepp'], filt_raref, jobs)
        AnalysisPrep('sepp').sepp(config, project)

    if filt_only and datasets_filt_map:
        deleted_non_filt(
            datasets, datasets_read, datasets_features, datasets_phylo,
            datasets_rarefs, taxonomies, datasets_filt, datasets_filt_map)
    if filt_only:
        project.delete_non_filtered()

    split_taxa_pds = get_taxo_levels(taxonomies)
    project.get_taxo_levels()

    if 'do_pies' in p_skip:
        print('(run_do_pies)')
        pies_data = make_pies(
            i_datasets_folder, split_taxa_pds, datasets_rarefs, datasets_read)

    collapsed = {}
    datasets_collapsed = {}
    datasets_collapsed_map = {}
    if p_collapse_taxo and 'collapse' not in p_skip:
        print('(run_collapse)')
        collapsed = run_collapse(
            i_datasets_folder, datasets, datasets_filt, datasets_read,
            datasets_features, datasets_phylo, split_taxa_pds, taxonomies,
            p_collapse_taxo, datasets_rarefs, datasets_collapsed,
            datasets_collapsed_map, force, prjct_nm, qiime_env, p_chmod, noloc,
            run_params["collapse"], filt_raref, jobs)
        AnalysisPrep('collapse').collapse(config, project)

    # ALPHA ------------------------------------------------------------
    alpha = AlphaDiversity(config, project)
    if 'alpha' not in p_skip:
        print('(alpha)')
        diversities = run_alpha(
            i_datasets_folder, datasets, datasets_read, datasets_phylo,
            datasets_rarefs, p_alpha_subsets, trees, force, prjct_nm,
            qiime_env, p_chmod, noloc, As, dropout, run_params['alpha'],
            filt_raref, eval_depths, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

        if 'merge_alpha' not in p_skip:
            print('(to_export)')
            to_export = merge_meta_alpha(
                i_datasets_folder, datasets, datasets_rarefs, diversities,
                force, prjct_nm, qiime_env, p_chmod, noloc,  dropout,
                run_params['merge_alpha'], filt_raref, eval_depths, jobs,
                chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
            if 'export_alpha' not in p_skip:
                print('(export_meta_alpha)')
                export_meta_alpha(
                    datasets, filt_raref, datasets_rarefs, to_export, dropout)
                # AnalysisPrep('Alpha diversity').alpha(config, project)
                # print(AnalysisPrep.analyses_commands)
                # print(AnalysisPrep.analyses_commands.keys())
                # print(AnalysisPrepdsa)
        if 'alpha_correlations' not in p_skip:
            print('(run_correlations)')
            run_correlations(
                i_datasets_folder, datasets, diversities, datasets_rarefs,
                force, prjct_nm, qiime_env, p_chmod, noloc,
                run_params['alpha_correlations'], filt_raref, jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
        if 'alpha_rarefaction' not in p_skip:
            print('(run_alpha_rarefaction)')
            run_alpha_rarefaction(
                i_datasets_folder, datasets, datasets_rarefs, datasets_phylo,
                trees, force, prjct_nm, qiime_env, p_chmod, noloc, As,
                run_params['alpha_rarefaction'], filt_raref, jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
        if p_longi_column and 'volatility' not in p_skip:
            print('(run_volatility)')
            run_volatility(
                i_datasets_folder, datasets, p_longi_column,
                datasets_rarefs, force, prjct_nm, qiime_env, p_chmod,
                noloc, run_params['volatility'], filt_raref, jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)

    # BETA ----------------------------------------------------------------------
    if 'beta' not in p_skip:
        print('(betas)')
        betas = run_beta(
            i_datasets_folder, datasets, datasets_phylo, datasets_read,
            datasets_rarefs, p_beta_subsets, p_beta_groups, trees, force,
            prjct_nm, qiime_env, p_chmod, noloc, Bs, dropout, run_params['beta'],
            filt_raref, eval_depths, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)
        if 'export_beta' not in p_skip:
            print('(export_beta)')
            export_beta(
                i_datasets_folder, betas, datasets_rarefs, force, prjct_nm,
                qiime_env, p_chmod, noloc, run_params['export_beta'],
                filt_raref, jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
        if 'pcoa' not in p_skip:
            print('(run_pcoas)')
            pcoas = run_pcoas(
                i_datasets_folder, betas, datasets_rarefs, force, prjct_nm,
                qiime_env, p_chmod, noloc, run_params['pcoa'], filt_raref,
                jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
            if 'emperor' not in p_skip:
                print('(run_emperor)')
                run_emperor(
                    i_datasets_folder, pcoas, datasets_rarefs, prjct_nm,
                    qiime_env, p_chmod, noloc, run_params['emperor'],
                    filt_raref, jobs, chunkit)
                # AnalysisPrep('Alpha diversity').alpha(config, project)
                # print(AnalysisPrep.analyses_commands)
                # print(AnalysisPrep.analyses_commands.keys())
                # print(AnalysisPrepdsa)
            if 'empress' not in p_skip:
                print('(run_empress)')
                run_empress(
                    i_datasets_folder, pcoas, trees, datasets_phylo,
                    datasets_rarefs, taxonomies, prjct_nm, qiime_env, p_chmod,
                    noloc, run_params['empress'], filt_raref, jobs, chunkit)
                # AnalysisPrep('Alpha diversity').alpha(config, project)
                # print(AnalysisPrep.analyses_commands)
                # print(AnalysisPrep.analyses_commands.keys())
                # print(AnalysisPrepdsa)
        if do_biplots and 'biplot' not in p_skip:
            print('(run_biplots)')
            biplots, biplots_raw = run_biplots(
                i_datasets_folder, betas, datasets_rarefs,  taxonomies,
                force, prjct_nm, qiime_env, p_chmod, noloc,
                run_params['biplot'], filt_raref, jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
            if 'emperor_biplot' not in p_skip:
                print('(run_emperor_biplot)')
                run_emperor_biplot(
                    i_datasets_folder, biplots, biplots_raw, taxonomies,
                    split_taxa_pds, datasets_rarefs, prjct_nm, qiime_env,
                    p_chmod, noloc, run_params['emperor_biplot'], filt_raref,
                    jobs, chunkit)
                # AnalysisPrep('Alpha diversity').alpha(config, project)
                # print(AnalysisPrep.analyses_commands)
                # print(AnalysisPrep.analyses_commands.keys())
                # print(AnalysisPrepdsa)
            if 'empress_biplot' not in p_skip:
                print('(run_empress_biplot)')
                run_empress_biplot(
                    i_datasets_folder, biplots, biplots_raw, trees,
                    datasets_phylo, taxonomies, datasets_rarefs, prjct_nm,
                    qiime_env, p_chmod, noloc, run_params['empress_biplot'],
                    filt_raref, jobs, chunkit)
                # AnalysisPrep('Alpha diversity').alpha(config, project)
                # print(AnalysisPrep.analyses_commands)
                # print(AnalysisPrep.analyses_commands.keys())
                # print(AnalysisPrepdsa)

    # STATS ------------------------------------------------------------------
    if 'alpha' not in p_skip and 'alpha_group_significance' not in p_skip and 'alpha_kw' not in p_skip:
        print('(run_alpha_group_significance)')
        run_alpha_group_significance(
            i_datasets_folder, datasets, diversities, datasets_rarefs,
            p_beta_groups, force, prjct_nm, qiime_env, p_chmod, noloc, As,
            split, run_params['alpha_kw'], filt_raref, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and 'deicode' not in p_skip:
        print('(run_deicode)')
        run_deicode(
            i_datasets_folder, datasets, datasets_rarefs, p_beta_groups,
            force, prjct_nm, qiime_env, p_chmod, noloc, run_params['deicode'],
            filt_raref, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and p_perm_tests and 'permanova' not in p_skip:
        print('(run_permanova)')
        permanovas = run_permanova(
            i_datasets_folder, betas, p_perm_tests, p_perm_tests_min,
            p_beta_type, datasets_rarefs, p_beta_groups, force, prjct_nm,
            qiime_env, p_chmod, noloc, split, run_params['permanova'],
            filt_raref, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

        summarize_permanova(
            i_datasets_folder, permanovas, prjct_nm, qiime_env, p_chmod, noloc,
            split, run_params['permanova'], filt_raref, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and p_adonis_formulas and 'adonis' not in p_skip:
        print('(run_adonis)')
        run_adonis(
            p_adonis_formulas, i_datasets_folder, betas, datasets_rarefs,
            p_beta_groups, force, prjct_nm, qiime_env, p_chmod, noloc, split,
            run_params['adonis'], filt_raref, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and p_procrustes and 'procrustes' not in p_skip:
        print('(run_procrustes)')
        run_procrustes(
            i_datasets_folder, datasets_filt, p_procrustes, betas, force,
            prjct_nm, qiime_env, p_chmod, noloc, split,
            run_params['procrustes'], filt_raref, filt_only, eval_depths,
            jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and p_mantel and 'mantel' not in p_skip:
        print('(run_mantel)')
        run_mantel(
            i_datasets_folder, datasets_filt, p_mantel, betas, force,
            prjct_nm, qiime_env, p_chmod, noloc, split, run_params['mantel'],
            filt_raref, filt_only, eval_depths, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and p_nestedness_groups and 'nestedness' not in p_skip:
        print('(run_nestedness)')
        nestedness_res, colors, nodfs_fps = run_nestedness(
            i_datasets_folder, betas, datasets_collapsed_map,
            p_nestedness_groups, datasets_rarefs, force, prjct_nm, qiime_env,
            p_chmod, noloc, split, run_params['nestedness'], filt_raref, jobs,
            chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)
        if nestedness_res:
            print('(making_nestedness_figures (graphs))')
            nestedness_graphs(
                i_datasets_folder, nestedness_res, datasets, split_taxa_pds,
                datasets_rarefs, colors, datasets_collapsed_map, collapsed,
                filt_raref, prjct_nm, qiime_env, p_chmod, noloc, split,
                run_params['nestedness'], jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
        if nodfs_fps:
            print('(making_nestedness_figures (nodfs))')
            nestedness_nodfs(
                i_datasets_folder, nodfs_fps, collapsed, filt_raref,
                prjct_nm, qiime_env, p_chmod, noloc, split,
                run_params['nestedness'], jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)

    if 'beta' not in p_skip and p_distance_decay and 'decay' not in p_skip:
        print('(run_distance_decay)')
        distance_decay_res = run_distance_decay(
            i_datasets_folder, betas, p_distance_decay, datasets_rarefs,
            force, prjct_nm, qiime_env, p_chmod, noloc, split,
            run_params['decay'], filt_raref, jobs, chunkit)
        # if distance_decay_res:
        #     print('(making_distance_decay_figures)')
        #     distance_decay_figure(i_datasets_folder, distance_decay_res,
        #                           datasets_rarefs, filt_raref)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    # PHATE ---------------------------------------------------------------------
    if p_phate_config and 'phate' not in p_skip:
            print('(run_phate)')
            phates = run_phate(
                p_phate_config, i_datasets_folder, datasets, datasets_rarefs,
                force, prjct_nm, qiime_env, p_chmod, noloc, split,
                run_params['phate'], filt_raref, jobs, chunkit)
            # AnalysisPrep('Alpha diversity').alpha(config, project)
            # print(AnalysisPrep.analyses_commands)
            # print(AnalysisPrep.analyses_commands.keys())
            # print(AnalysisPrepdsa)
    else:
        phates = {}

    # DISSIMILARITY OVERLAP --------------------------------------------
    if 'doc' not in p_skip and p_doc_config:
        print('(run_doc)')
        run_doc(
            i_datasets_folder, datasets, p_doc_config, datasets_rarefs,
            force, prjct_nm, qiime_env, p_chmod, noloc, run_params['doc'],
            filt_raref, phates, doc_phate, split, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    # SOURCETRACKING --------------------------------------------
    if p_sourcetracking_config and 'sourcetracking' not in p_skip:
        print('(run_sourcetracking)')
        run_sourcetracking(
            i_datasets_folder, datasets, p_sourcetracking_config,
            datasets_rarefs, force, prjct_nm, qiime_env, p_chmod, noloc,
            run_params['sourcetracking'], filt_raref, split, jobs, chunkit)
        # AnalysisPrep('Alpha diversity').alpha(config, project)
        # print(AnalysisPrep.analyses_commands)
        # print(AnalysisPrep.analyses_commands.keys())
        # print(AnalysisPrepdsa)

    # MMVEC AND SONGBIRD --------------------------------------------------------
    train_test_dict = {}
    if p_train_test:
        # read config to create train-test metadata columns
        train_test_dict = get_train_test_dict(p_train_test)
    config.get_train_test_dict()

    filts = {}
    input_to_filtered = {}
    mmvec_outputs = []
    paired_datasets = PairedData(config, project)
    if p_mmvec_pairs:
        if filt3d:
            filts.update(get_filt3d_params(p_mmvec_pairs, 'mmvec'))
        elif 'mmvec' not in p_skip:
            # print('(run_mmvec)')
            # mmvec_outputs = run_mmvec(
            #     p_mmvec_pairs, i_datasets_folder, datasets, datasets_filt,
            #     datasets_read, train_test_dict, force, gpu, standalone,
            #     prjct_nm, qiime_env, p_chmod, noloc, split, filt_raref,
            #     run_params['mmvec'], input_to_filtered, jobs, chunkit)
            paired_datasets.make_train_test(config)
            paired_datasets.mmvec(config)

    songbird_outputs = []
    differentials = DiffModels(config, project)
    if p_diff_models:
        if filt3d:
            filts.update(get_filt3d_params(p_diff_models, 'songbird'))
        elif 'songbird' not in p_skip:
            # print('(run_songbird)')
            # songbird_outputs = run_songbird(
            #     p_diff_models, i_datasets_folder, datasets, datasets_read,
            #     datasets_filt, train_test_dict, input_to_filtered,
            #     mmvec_outputs, force, prjct_nm, qiime_env, p_chmod, noloc,
            #     split, run_params['songbird'], filt_raref, jobs, chunkit)
            differentials.get_songbirds_matrix(paired_datasets.mmvec_pd)
            differentials.make_train_test(config)
            differentials.songbird(config, project)
            # q2s_pd = summarize_songbirds(config.i_datasets_folder)
            # out_folder = get_analysis_folder(i_datasets_folder, 'songbird')
            # q2s_fp = '%s/songbird_q2.tsv' % out_folder
            # q2s_pd.to_csv(q2s_fp, index=False, sep='\t')
            # print('\t\t==> Written:', q2s_fp)
            # create_songbird_feature_metadata(
            #     i_datasets_folder, taxonomies, q2s_pd)

    if filt3d:
        print('(run_filt3d)')
        explore_filtering(
            i_datasets_folder, datasets, datasets_read, datasets_filt,
            datasets_filt_map, filts, p_filt3d_config)
    elif p_mmvec_pairs and 'mmbird' not in p_skip:
        # print('(run_mmbird)')
        # run_mmbird(
        #     i_datasets_folder, songbird_outputs, p_mmvec_highlights,
        #     p_xmmvec, mmvec_outputs, force, prjct_nm, qiime_env, p_chmod,
        #     noloc, filt_raref, run_params['mmbird'],
        #     input_to_filtered, jobs, chunkit)
        post_analyses = PostAnalysis(config, project)
        post_analyses.mmbird(paired_datasets, differentials)

    scripting.write_scripts(AnalysisPrep.analyses_commands)
