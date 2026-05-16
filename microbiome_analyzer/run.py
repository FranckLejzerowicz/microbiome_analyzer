# ----------------------------------------------------------------------------
# Copyright (c) 2026, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import logging
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
from microbiome_analyzer.analyses.sourcetracking import Sourcetracking
# from routine_qiime2_analyses._routine_q2_qemistree import run_qemistree


def set_log():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    # Using StreamHandler writing to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    # Add the two Handlers
    logger.addHandler(ch)


def runner(**kwargs):
    """
    Main qiime2 functions writer.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    set_log()
    kwargs['logging'] = logging
    kwargs['command'] = 'run'
    logging.info('\n>>> `microbiome_analyzer run` started >>>\n')

    logging.info('* Reading configuration for tools and analyses')
    config = AnalysesConfig(**kwargs)
    config.init(logging)

    logging.info('* Collect the data associated with each dataset')
    project = Datasets(config)
    project.collect_datasets(logging)

    logging.info('* Preparing analyses for:')
    analysis = AnalysisPrep(config, project)
    logging.info('  - Importing datasets to QIIME2')
    analysis.import_datasets()
    if config.filt3d:
        logging.info('  - Exploring feature counts (vs prevalence/abundance)')
        analysis.explore_filtering()
    logging.info('  - Filtering samples and features')
    analysis.filter()
    if config.rarefy:
        logging.info('  - Rarefying features in samples')
        analysis.rarefy()
    # if config.qemistree and 'qemistree' not in config.skip:
    #     analysis.run_qemistree()
    logging.info('  - Obtaining taxonomy (user-defined or analysis)')
    project.get_precomputed_taxonomy()
    if 'taxonomy' not in config.skip:
        analysis.taxonomy()
    logging.info('  - Get taxonomic levels')
    project.get_taxo_levels()
    logging.info('  - Obtaining phylogeny')
    project.get_precomputed_trees()
    analysis.import_trees()

    if 'wol' not in config.skip:
        logging.info('  - Shear the Web Of Life phylogeny (woltka input only)')
        analysis.shear_tree()
    if config.sepp_tree and 'sepp' not in config.skip:
        logging.info('  - Placing ASV into reference phylogeny')
        analysis.sepp()
    if config.filt_only:
        project.delete_non_filtered()
    # if 'pies' not in config.skip:
    #     analysis.make_pies()

    if config.feature_subsets and 'feature_subsets' not in config.skip:
        logging.info('  - Making feature subsets')
        analysis.subset_features()

    if config.collapse and 'collapse' not in config.skip:
        logging.info('  - Collapsing features based on taxonomy')
        analysis.collapse_taxa()

    analysis.show_datasets()
    analysis.make_edits()

    if 'alpha' not in config.skip:
        logging.info('  - Measuring alpha diversity')
        analysis.alpha()
        if 'alpha_correlations' not in config.skip:
            logging.info('    - Correlations vs numeric metadata variables')
            analysis.alpha_correlations()
        if 'alpha_merge' not in config.skip:
            logging.info('    - Merging across diversity indices')
            analysis.merge_alpha()
            analysis.merge_metadata()
        if 'alpha_rarefaction' not in config.skip:
            logging.info('    - Drawing rarefaction curves')
            analysis.alpha_rarefaction()
        if 'alpha_group_significance' not in config.skip:
            logging.info('    - Testing for categorical metadata variables')
            analysis.alpha_group_significance()

    logging.info('  - Making sample subsets')
    project.get_meta_subsets()

    if config.phate and 'phate' not in config.skip:
        logging.info('  - Running PHATE ordinations')
        analysis.phate()
    if 'barplot' not in config.skip:
        logging.info('  - Drawing barplots')
        analysis.barplot()
    if 'krona' not in config.skip:
        logging.info('  - Drawing Krona plots')
        analysis.krona()

    if config.predict and 'classify' not in config.skip:
        logging.info('  - Predicting metadata (classification and regression)')
        analysis.classify()

    if 'beta' not in config.skip:
        logging.info('  - Measuring beta diversity')
        analysis.beta()
        if 'pcoa' not in config.skip:
            logging.info('    - Making PCoA ordinations')
            analysis.pcoa()
        if 'tsne' not in config.skip:
            logging.info('    - Making t-SNE ordinations')
            analysis.tsne()
        if 'umap' not in config.skip:
            logging.info('    - Making UMAP ordinations')
            analysis.umap()
        if 'emperor' not in config.skip:
            logging.info('    - Plotting ordinations in EMPeror')
            analysis.emperor()
        if config.biplot and 'biplot' not in config.skip:
            logging.info('    - Making PCoA ordinations (incl. features)')
            analysis.biplots()
            if 'emperor_biplot' not in config.skip:
                logging.info('    - Plotting ordinations in EMPeror (biplots)')
                analysis.emperor_biplot()
        if 'rpca' not in config.skip:
            logging.info('    - Running Robust PCA')
            analysis.rpca()
        if config.permanova and 'permanova' not in config.skip:
            logging.info('    - Testing with PERMANOVA/PERMDISP (per variable)')
            analysis.permanova()
        # summarize_permanova(
        #     datasets_folder, permanovas, prjct_nm, qiime_env, p_chmod, noloc,
        #     slurm, split, run_params['permanova'], filt_raref, jobs, chunkt)
        if config.adonis and 'adonis' not in config.skip:
            logging.info('    - Testing PERMANOVA in adonis (>1 variable)')
            analysis.adonis()
        if config.procrustes and 'procrustes' not in config.skip:
            logging.info('    - Procrustes and PROTests')
            analysis.procrustes_mantel('procrustes')
        if config.mantel and 'mantel' not in config.skip:
            logging.info('    - Mantel correlations tests')
            analysis.procrustes_mantel('mantel')
        if config.nestedness and 'nestedness' not in config.skip:
            logging.info('    - Measuring nestedness')
            Nestedness(config, project)
        if config.dm_decay and 'dm_decay' not in config.skip:
            logging.info('    - Measuring Distance-based decay vs sampling')
            DmDecay(config, project)
        if config.geo_decay and 'geo_decay' not in config.skip:
            logging.info('    - Measuring Distance-based decay vs geography')
            pass
        if config.time_subject:
            if 'ctf' not in config.skip:
                logging.info('    - CTF analysis')
                analysis.ctf()
            if 'volatility' not in config.skip:
                logging.info('    - Temporal plots')
                analysis.volatility()

    if config.doc and 'doc' not in config.skip:
        logging.info('  - Dissimilarity-Overlap curves')
        DOC(config, project)
    if config.sourcetracking and 'sourcetracking' not in config.skip:
        logging.info('  - Sourcetracking (sourcetracker2 and FEAST)')
        Sourcetracking(config, project)

    paired_datasets = PairedData(config, project)
    if config.mmvec_pairs and 'mmvec' not in config.skip:
        logging.info('  - Co-occurrence probability measurement with MMVEC')
        paired_datasets.make_train_test()
        paired_datasets.mmvec()

    differentials = DiffModels(config, project)
    if config.diff_abund and 'songbird' not in config.skip:
        logging.info('  - Differential abundance measurement with SONGBIRD')
        differentials.prep_songbirds(paired_datasets.mmvec_pd)
        differentials.make_train_test()
        differentials.songbird()
        if 'qurro' not in config.skip:
            differentials.make_qurros()

    if 'empress' not in config.skip:
        logging.info('  - Plotting phylogeny incl. feature metadata')
        analysis.empress()
    # if 'empress_biplot' not in config.skip:
    #     analysis.empress_biplot()

    if config.mmvec_pairs and 'mmbird' not in config.skip:
        logging.info('  - Integrating differentials in co-occurrence biplots')
        post_analyses = PostAnalysis(config, project)
        post_analyses.mmbird(paired_datasets, differentials)

    scripting = CreateScripts(config, project)
    logging.info('* Creating output folders')
    scripting.make_dirs()
    logging.info('* Writing command lines')
    scripting.writing(AnalysisPrep)
    if len(scripting.run):
        m = '\n< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY >'
        logging.info(m)
        scripting.display()  # show the scripts to run
    # scripting.write_scripts(AnalysisPrep.analyses_commands)
    logging.info('\n<<< `microbiome_analyzer run` completed <<<\n')
