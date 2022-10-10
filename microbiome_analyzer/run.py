# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from microbiome_analyzer.core.config import AnalysesConfig
from microbiome_analyzer.core.datasets import Datasets
from microbiome_analyzer.core.jobs import CreateScripts
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.analyses.nestedness import Nestedness
from microbiome_analyzer.analyses.dm_decay import DmDecay
from microbiome_analyzer.analyses.doc import DOC
from microbiome_analyzer.analyses.post_analyses import PostAnalysis
from microbiome_analyzer.analyses.differential_abundance import DiffModels
from microbiome_analyzer.analyses.paired_data import PairedData
# from microbiome_analyzer.analyses.sourcetracking import Sourcetracking
# from routine_qiime2_analyses._routine_q2_qemistree import run_qemistree


def runner(**kwargs):
    """
    Main qiime2 functions writer.
    """
    config = AnalysesConfig(**kwargs)
    config. init()

    project = Datasets(config)
    project.collect_datasets()

    analysis = AnalysisPrep(config, project)
    analysis.import_datasets()
    if config.filt3d:
        analysis.explore_filtering()
    analysis.filter()
    if config.rarefy:
        analysis.rarefy()

    # if config.qemistree and 'qemistree' not in config.skip:
    #     analysis.run_qemistree()

    project.get_precomputed_taxonomy()
    if 'taxonomy' not in config.skip:
        analysis.taxonomy()
    project.get_taxo_levels()
    project.get_precomputed_trees()
    analysis.import_trees()

    if 'wol' not in config.skip:
        analysis.shear_tree()
    if config.sepp_tree and 'sepp' not in config.skip:
        analysis.sepp()
    if config.filt_only:
        project.delete_non_filtered()
    # if 'pies' not in config.skip:
    #     analysis.make_pies()

    if config.feature_subsets and 'feature_subsets' not in config.skip:
        analysis.subset_features()
    if config.collapse and 'collapse' not in config.skip:
        analysis.collapse_taxa()
    analysis.edits()
    if 'alpha' not in config.skip:
        analysis.alpha()
        if 'alpha_correlations' not in config.skip:
            analysis.alpha_correlations()
        if 'alpha_merge' not in config.skip:
            analysis.merge_alpha()
            analysis.merge_metadata()
        if 'alpha_rarefaction' not in config.skip:
            analysis.alpha_rarefaction()
        if config.longi_column and 'volatility' not in config.skip:
            analysis.volatility()
        if 'alpha_group_significance' not in config.skip:
            # make a python script to fill with groups to test
            pass

    project.get_meta_subsets()

    if config.phate and 'phate' not in config.skip:
        analysis.phate()
    if 'barplot' not in config.skip:
        analysis.barplot()

    if 'beta' not in config.skip:
        analysis.beta()
        if 'deicode' not in config.skip:
            analysis.deicode()
        if 'pcoa' not in config.skip:
            analysis.pcoa()
        if 'tsne' not in config.skip:
            analysis.tsne()
        if 'umap' not in config.skip:
            analysis.umap()
        if 'emperor' not in config.skip:
            analysis.emperor()
        if config.biplot and 'biplot' not in config.skip:
            analysis.biplots()
            if 'emperor_biplot' not in config.skip:
                analysis.emperor_biplot()
        if config.tests and 'permanova' not in config.skip:
            analysis.permanova()
            # analysis.permanova_r()
        # summarize_permanova(
        #     datasets_folder, permanovas, prjct_nm, qiime_env, p_chmod, noloc,
        #     slurm, split, run_params['permanova'], filt_raref, jobs, chunkt)
        if config.adonis and 'adonis' not in config.skip:
            analysis.adonis()
        if config.procrustes and 'procrustes' not in config.skip:
            analysis.procrustes_mantel('procrustes')
        if config.mantel and 'mantel' not in config.skip:
            analysis.procrustes_mantel('mantel')
        if config.nestedness and 'nestedness' not in config.skip:
            Nestedness(config, project)
        if config.dm_decay and 'dm_decay' not in config.skip:
            DmDecay(config, project)
        if config.geo_decay and 'geo_decay' not in config.skip:
            pass

    if config.doc and 'doc' not in config.skip:
        DOC(config, project)
    # if config.sourcetracking and 'sourcetracking' not in config.skip:
    #     Sourcetracking(config, project)

    paired_datasets = PairedData(config, project)
    print()
    print("paired_datasets.mmvec_pd")
    print(paired_datasets.mmvec_pd)
    if config.mmvec_pairs and 'mmvec' not in config.skip:
        paired_datasets.make_train_test()
        paired_datasets.mmvec()

    differentials = DiffModels(config, project)
    if config.diff_models and 'songbird' not in config.skip:
        differentials.prep_songbirds(paired_datasets.mmvec_pd)
        differentials.make_train_test()
        differentials.songbird()
        if 'qurro' not in config.skip:
            differentials.make_qurros()

    if 'empress' not in config.skip:
        analysis.empress()
    if 'empress_biplot' not in config.skip:
        analysis.empress_biplot()

    elif config.mmvec_pairs and 'mmbird' not in config.skip:
        post_analyses = PostAnalysis(config, project)
        post_analyses.mmbird(paired_datasets, differentials)

    scripting = CreateScripts(config, project)
    scripting.make_dirs()
    scripting.writing(AnalysisPrep)
    if len(scripting.run):
        m = '\n< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY >'
        print(m)
        scripting.display()  # show the scripts to run
    # scripting.write_scripts(AnalysisPrep.analyses_commands)
