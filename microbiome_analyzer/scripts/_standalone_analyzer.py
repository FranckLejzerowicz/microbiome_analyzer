# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from microbiome_analyzer.microbiome_analyzer import microbiome_analyzer
from microbiome_analyzer import __version__


@click.command()
@click.option(
    "-i", "--data-folder", required=True,
    help="Path to the folder containing the data tables")
@click.option(
    "-j", "--metadata-folder", required=True,
    help="Path to the folder containing the metadata tables")
@click.option(
    "-o", "--output-folder", required=True,
    help="Path to a folder containing the outputs")
@click.option(
    "-d", "--datasets", multiple=True, required=True,
    help="Dataset(s) identifier(s). Multiple is possible: e.g. "
         "-d dataset_number_1 and -d dataset_number_2 for "
         "'tab_dataset_number_1.tsv' and tab_dataset_number_2.tsv'")
@click.option(
    "-n", "--project-name", required=True, show_default=True,
    help="Nick name for your project")
@click.option(
    "-e", "--qiime2-env", required=True, show_default=True,
    help="name of your qiime2 conda environment (e.g. qiime2-2021.11)")
@click.option(
    "-w", "--wol-tree", required=False, show_default=True,
    default='resources/wol_tree.nwk',
    help="path to the tree containing the genome IDs (will check if exist in "
         "features names) (On barnacle, it is there: "
         "/projects/wol/profiling/dbs/wol/phylogeny/tree.nwk)")
@click.option(
    "-x", "--sepp-tree", required=False, show_default=True, default=None,
    help="Qiime2 SEPP reference database to use for 16S reads placement: "
         "https://docs.qiime2.org/2019.10/data-resources/#sepp-reference-"
         "databases (auto detection of datasets' tables with sequences as "
         "features)")
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
    "-t", "--test", multiple=True, required=False, show_default=True,
    default=False, help="Groups to tests between in each PERMANOVA subset "
                        "(multiple values possible, e.g. '-d sex -d age_cat')")
@click.option(
    "-g", "--sample-subsets", required=False, show_default=True,
    default=False, help="Subsets for DMs, PCoAs, PERMANOVAs, etc (yml file)")
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
    "-a", "--adonis", required=False, default=False, show_default=True,
    help="Formula for Adonis tests for each PERMANOVA subset (yaml file)")
@click.option(
    "-phate", "--phate", required=False, default=False, show_default=True,
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
    "-s", "--diff-models", required=False, default=False, show_default=True,
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
    "-l", "--longi-column", required=False, default=False, show_default=True,
    help="If data is longitudinal; provide the time metadata column"
         "for volatility analysis")
@click.option(
    "-f", "--filter", show_default=True, default=False,
    help="Samples to remove, min. read abundance and feature prevalence "
         "(>1 = based on absolute reads, 0-1 = based on relative reads). "
         "(yaml file)")
@click.option(
    "-r", "--rarefs", required=False, show_default=False,
    help="rarefaction depth per dataset (yaml file)")
@click.option(
    "-chmod", "--chmod", default='664', show_default=True,
    help="Change output files permission (default = 664 [= -rw-rw-r--])")
@click.option(
    "-skip", "--skip", default=None, show_default=True, multiple=True,
    help="Steps to skip (e.g. if already done or not necessary)",
    type=click.Choice(['alpha', 'merge_alpha', 'export_alpha',
                       'alpha_correlations', 'alpha_group_significance',
                       'wol', 'taxonomy', 'barplot', 'volatility', 'beta',
                       'export_beta', 'pcoa', 'biplot', 'emperor',
                       'emperor_biplot', 'empress', 'empress_biplot',
                       'phate', 'doc', 'deicode', 'sepp', 'do_pies',
                       'alpha_kw', 'permanova', 'procrustes', 'mantel', 'decay',
                       'nestedness', 'adonis', 'songbird', 'mmvec', 'mmbird']))
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
    "--raref/--no-raref", default=False, show_default=True,
    help="Whether to rarefy and only perform the routine "
         "analyses on the rarefied dataset(s)")
@click.option(
    "--loc/--no-loc", default=True, show_default=True,
    help="whether to do compute on scratch")
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
    "-chunks", "--chunks", required=False, show_default=False,
    type=int, default=None,
    help="Maximum number of jobs at which extra jobs will be added in chunks")
@click.version_option(__version__, prog_name="microbiome_analyzer")


def standalone_analyzer(
        datasets,
        data_folder,
        metadata_folder,
        output_folder,
        project_name,
        qiime2_env,
        longi_column,
        filter,
        rarefs,
        feature_subsets,
        test,
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
        diff_models,
        mmvec_pairs,
        mmvec_highlights,
        xmmvec,
        run_params,
        chmod,
        skip,
        gpu,
        raref,
        loc,
        alphas,
        betas,
        account,
        filt3d,
        filt_only,
        jobs,
        torque,
        chunks
):

    microbiome_analyzer(
        datasets=datasets,
        data_folder=data_folder,
        metadata_folder=metadata_folder,
        output_folder=output_folder,
        project_name=project_name,
        qiime_env=qiime2_env,
        longi_column=longi_column,
        raref=raref,
        filter_fp=filter,
        rarefs_fp=rarefs,
        feature_subsets_fp=feature_subsets,
        tests=test,
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
        diff_models_fp=diff_models,
        mmvec_pairs_fp=mmvec_pairs,
        mmvec_highlights_fp=mmvec_highlights,
        xmmvec_fp=xmmvec,
        run_params_fp=run_params,
        chmod=chmod,
        skip=skip,
        gpu=gpu,
        loc=loc,
        alphas=alphas,
        betas=betas,
        account=account,
        filt3d_fp=filt3d,
        filt_only=filt_only,
        jobs=jobs,
        torque=torque,
        chunkit=chunks
    )


if __name__ == "__main__":
    standalone_analyzer()
