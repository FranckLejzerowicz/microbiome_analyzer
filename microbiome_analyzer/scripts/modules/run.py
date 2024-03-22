# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from microbiome_analyzer.run import runner
from microbiome_analyzer import __version__


@click.command()
@click.option(
    "-i", "--analysis-folder", required=True,
    help="Path to the folder containing the data and metadata sub-folders")
@click.option(
    "-d", "--datasets", multiple=True, required=True,
    help="Dataset(s) identifier(s). Multiple is possible: e.g. "
         "-d dataset_number_1 and -d dataset_number_2 for "
         "'table.tsv' and table.tsv'")
@click.option(
    "-n", "--project-name", required=True, show_default=True,
    help="Nick name for your project")
@click.option(
    "-f", "--filtering", show_default=True, default=False,
    help="Samples to remove, min. read abundance and feature prevalence "
         "(>1 = based on absolute reads, 0-1 = based on relative reads). "
         "(yaml file)")
@click.option(
    "-r", "--rarefactions", required=False, show_default=False,
    help="rarefaction depth per dataset (yaml file)")
@click.option(
    "-e", "--qiime2-env", required=True, show_default=True,
    help="name of your qiime2 conda environment (e.g. qiime2-2021.11)")
@click.option(
    "-v", "--sepp-tree", required=False, show_default=True, default=None,
    help="Qiime2 SEPP reference database to use for 16S reads placement: "
         "https://docs.qiime2.org/2019.10/data-resources/#sepp-reference-"
         "databases (auto detection of datasets' tables with sequences as "
         "features)")
@click.option(
    "-w", "--wol-tree", required=False, show_default=True,
    default='resources/wol_tree.nwk',
    help="path to the tree containing the genome IDs (will check if exist in "
         "features names) (On barnacle, it is there: "
         "/projects/wol/profiling/dbs/wol/phylogeny/tree.nwk)")
@click.option(
    "-q", "--qemistree", required=False, show_default=True, default=None,
    help="Path to a folder containing Qemistree's feature data "
         "(named 'feature-data_<dataset_identifier>.qza'), "
         "and tree for each metabolomics dataset "
         "(named 'qemistree_<dataset_identifier>.qza')")
@click.option(
    "-z", "--classifier", required=False, show_default=True, default=None,
    help="Qiime2 reference taxonomic classifier database to use for 16S"
         "reads assignment: https://docs.qiime2.org/2020.2/"
         "data-resources/#taxonomy-classifiers-for-use-with-q2-"
         "feature-classifier")
@click.option(
    "-u", "--run-params", required=False, show_default=True, default=None,
    help="server run parameters")
@click.option(
    "-k", "--feature-subsets", required=False, show_default=True,
    default=None, help="Regex to use for subsetting features (yml file)")
@click.option(
    "-g", "--sample-subsets", required=False, show_default=True,
    default=False, help="Subsets for DMs, PCoAs, PERMANOVAs, etc (yml file)")
@click.option(
    "-t", "--time_subject", required=False, default=False, show_default=True,
    help="Time and subject variables for CTF/Volatility analyses (yaml file)")
@click.option(
    "-p", "--permanova", required=False, default=False, show_default=True,
    help="Groups for PERMANOVA/Kruskal-Wallis for each subset (yaml file)")
@click.option(
    "-a", "--adonis", required=False, default=False, show_default=True,
    help="Formula for Adonis2 tests for each subset (yaml file)")
@click.option(
    "-nstd", "--nestedness", required=False, show_default=True,
    default=False, help="Nestedness analysis config  (yml file)")
@click.option(
    "-bt", "--beta-type", required=False, show_default=False, multiple=True,
    default=('permanova', 'permdisp',), help="Type of beta group significance",
    type=click.Choice(['permanova', 'anosim', 'permdisp']))
@click.option(
    "-prc", "--procrustes", required=False, show_default=True, default=False,
    help="Pairs and subsets for procrustes/protests (yaml file)")
@click.option(
    "-mtl", "--mantel", required=False, show_default=True, default=False,
    help="Pairs and subsets for mantel test (yaml file)")
@click.option(
    "-ddecay", "--dm-decay", required=False, show_default=True, default=False,
    help="Parameters for (not geographic) distance decay analysis (yaml file)")
@click.option(
    "-gdecay", "--geo-decay", required=False, show_default=True,
    default=False, help="Parameters for geographic distance decay "
                        "analysis (yml file)")
@click.option(
    "-c", "--collapse", required=False, show_default=True, default=False,
    help="Nominative or rank-based taxonmic collapse per dataset (yaml file)")
@click.option(
    "-tt", "--train-test", required=False, show_default=True, default=False,
    help="Train test split per dataset (yaml file)")
@click.option(
    "-pht", "--phate", required=False, default=False, show_default=True,
    help="Filters, subsets, parameters and stratifications for the PHATE latent"
         " space analysis (yaml file)")
@click.option(
    "-st", "--sourcetracking", required=False, default=False,
    show_default=True, help="Filters, subsets, parameters and sink/sources for "
                            "sourcetracking (yaml file)")
@click.option(
    "-doc", "--doc", required=False, default=False, show_default=True,
    help="Filters and subsets for the dissimilarity overlap curves analyses "
         " (yaml file)")
@click.option(
    "-s", "--diff-abund", required=False, default=False, show_default=True,
    help="Formulas for multinomial regression-based differential abundance "
         "ranking (songbird) (yaml file)")
@click.option(
    "-m", "--mmvec-pairs", required=False, default=False, show_default=True,
    help="Pairs of datasets for which to compute co-occurrences "
         "probabilities (mmvec) (yaml file)")
@click.option(
    "-hlg", "--mmvec-highlights", required=False, default=False,
    show_default=True, help="Features to highlights on mmvec biplot (per "
                            "dataset) (yaml file)")
@click.option(
    "-mm", "--xmmvec", required=False, default=False, show_default=True,
    help="Config for Xmmvec (yaml file)")
@click.option(
    "-chmod", "--chmod", default=None, show_default=True,
    help="Change output files permission (default = 664 [= -rw-rw-r--])")
@click.option(
    "-skip", "--skip", default=None, show_default=True, multiple=True,
    help="Steps to skip (e.g. if already done or not necessary)",
    type=click.Choice([
        'taxonomy',
        'barplot',
        'wol',
        'sepp',
        'pies',
        'collapse',
        'feature_subsets',
        'alpha',
        'alpha_merge',
        'alpha_rarefactions',
        'alpha_correlations',
        'alpha_group_significance',
        'volatility',
        'phate',
        'krona',
        'rpca',
        'beta',
        'deicode',
        'pcoa',
        'umap',
        'tsne',
        'emperor',
        'empress',
        'biplot',
        'emperor_biplot',
        'empress_biplot',
        'permanova',
        'adonis',
        'doc',
        'procrustes',
        'mantel',
        'nestedness',
        'dm_decay',
        'geo_decay',
        'sourcetracking',
        'doc',
        'mmvec',
        'songbird',
        'mmbird'
    ]))
@click.option(
    "-As", "--alphas", default=None, show_default=True, multiple=True,
    help="Alpha diversity indices to use")
@click.option(
    "-Bs", "--betas", default=None, show_default=True, multiple=True,
    help="Beta diversity metrics to use")
@click.option(
    "-acc", "--account", default=None, show_default=True,
    help="Account name")
@click.option(
    "--biplot/--no-biplot", default=False, show_default=True,
    help="Whether to do the PCoA biplots or not")
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists)")
@click.option(
    "--gpu/--no-gpu", default=False, show_default=True,
    help="Use GPUs instead of CPUs for MMVEC")
@click.option(
    "--rarefy/--no-rarefy", default=False, show_default=True,
    help="Whether to rarefy and only perform the routine "
         "analyses on the rarefied dataset(s)")
@click.option(
    "--filt-only/--no-filt-only", default=False, show_default=True,
    help="Only process the filtered version (and not also the raw) "
         "version of each dataset")
@click.option(
    "-filt3d", "--filt3d", required=False, show_default=False,
    help="Levels for the exploration of filtering. Must be a yaml file")
@click.option(
    "--jobs/--no-jobs", default=True, show_default=True,
    help="Whether to prepare Torque jobs from scripts")
@click.option(
    "--torque/--no-torque", default=False, show_default=True,
    help="Whether to prepare Torque and not Slurm jobs")
@click.option(
    "-l", "--localscratch", type=int, show_default=False, default=None,
    help="Use localscratch with the provided memory amount (in GB)")
@click.option(
    "--scratch/--no-scratch", default=False, show_default=True,
    help="Use the scratch folder to move files and compute")
@click.option(
    "--userscratch/--no-userscratch", default=False, show_default=True,
    help="Use the userscratch folder to move files and compute")
@click.option(
    "--move-back/--no-move-back", default=True, show_default=True,
    help="Do not move back from scratch (makes sense only for --userscratch)")
@click.option(
    "--cleanup/--no-cleanup", default=False, show_default=True,
    help="Whether to cleanup the TMPDIR and SCRATCH_FOLDER (specific to NRIS)")
@click.option(
    "-x", "--chunks", required=False, show_default=False,
    type=int, default=None,
    help="Maximum number of jobs at which extra jobs will be added in chunks")
@click.version_option(__version__, prog_name="microbiome_analyzer")

def run(
        datasets,
        analysis_folder,
        project_name,
        qiime2_env,
        filtering,
        rarefactions,
        feature_subsets,
        time_subject,
        permanova,
        sample_subsets,
        nestedness,
        beta_type,
        procrustes,
        mantel,
        dm_decay,
        geo_decay,
        collapse,
        train_test,
        adonis,
        doc,
        sourcetracking,
        phate,
        biplot,
        force,
        classifier,
        wol_tree,
        sepp_tree,
        qemistree,
        diff_abund,
        mmvec_pairs,
        mmvec_highlights,
        xmmvec,
        run_params,
        skip,
        chmod,
        gpu,
        rarefy,
        alphas,
        betas,
        account,
        filt3d,
        filt_only,
        jobs,
        torque,
        localscratch,
        scratch,
        userscratch,
        move_back,
        cleanup,
        chunks
):
    """Write jobs for your pipeline configuration."""
    runner(
        datasets=datasets,
        dir=analysis_folder,
        project_name=project_name,
        qiime_env=qiime2_env,
        rarefy=rarefy,
        filter_fp=filtering,
        rarefs_fp=rarefactions,
        feature_subsets_fp=feature_subsets,
        time_subject_fp=time_subject,
        permanova_fp=permanova,
        sample_subsets_fp=sample_subsets,
        nestedness_fp=nestedness,
        beta_type=beta_type,
        procrustes_fp=procrustes,
        mantel_fp=mantel,
        dm_decay_fp=dm_decay,
        geo_decay_fp=geo_decay,
        collapse_fp=collapse,
        train_test_fp=train_test,
        adonis_fp=adonis,
        doc_fp=doc,
        sourcetracking_fp=sourcetracking,
        phate_fp=phate,
        biplot=biplot,
        force=force,
        classifier=classifier,
        wol_tree=wol_tree,
        sepp_tree=sepp_tree,
        qemistree=qemistree,
        diff_abund_fp=diff_abund,
        mmvec_pairs_fp=mmvec_pairs,
        mmvec_highlights_fp=mmvec_highlights,
        xmmvec_fp=xmmvec,
        run_params_fp=run_params,
        skip=skip,
        chmod=chmod,
        gpu=gpu,
        alphas=alphas,
        betas=betas,
        account=account,
        filt3d_fp=filt3d,
        filt_only=filt_only,
        jobs=jobs,
        torque=torque,
        localscratch=localscratch,
        scratch=scratch,
        userscratch=userscratch,
        move_back=move_back,
        cleanup=cleanup,
        chunks=chunks
    )
