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
    "-i", "--i-datasets-folder", required=True,
    help="Path to the folder containing the sub-folders 'data' and 'metadata'."
)
@click.option(
    "-d", "--i-datasets", multiple=True, required=True,
    help="Dataset(s) identifier(s). Multiple is possible: e.g. "
         "-d dataset_number_1 and -d dataset_number_2 for "
         "'tab_dataset_number_1.tsv' and tab_dataset_number_2.tsv'."
)
@click.option(
    "-n", "--p-project-name", required=True, show_default=True,
    help="Nick name for your project."
)
@click.option(
    "-e", "--p-qiime2-env", required=True, show_default=True,
    help="name of your qiime2 conda environment (e.g. qiime2-2019.10) "
)
@click.option(
    "-w", "--i-wol-tree", required=False, show_default=True,
    default='resources/wol_tree.nwk',
    help="path to the tree containing the genome IDs (will check if exist in "
         "features names) (On barnacle, it is there: "
         "/projects/wol/profiling/dbs/wol/phylogeny/tree.nwk)."
)
@click.option(
    "-x", "--i-sepp-tree", required=False, show_default=True, default = None,
    help="Qiime2 SEPP reference database to use for 16S reads placement: "
         "https://docs.qiime2.org/2019.10/data-resources/#sepp-reference-"
         "databases (auto detection of datasets' tables with sequences as "
         "features)."
)
@click.option(
    "-q", "--i-qemistree", required=False, show_default=True, default = None,
    help="Path to a folder containing Qemistree's feature data "
         "(named 'feature-data_<dataset_identifier>.qza'), "
         "and tree for each metabolomics dataset "
         "(named 'qemistree_<dataset_identifier>.qza')."
)
@click.option(
    "-z", "--i-classifier", required=False, show_default=True, default=None,
    help="Qiime2 reference taxonomic classifier database to use for 16S"
         "reads assignment: https://docs.qiime2.org/2020.2/"
         "data-resources/#taxonomy-classifiers-for-use-with-q2-"
         "feature-classifier"
)
@click.option(
    "-u", "--p-run-params", required=False, show_default=True,
    help="server run paramters"
)
@click.option(
    "-k", "--p-feature-subsets", required=False, show_default=True, default=None,
    help="Regex to use for subsetting features (yml file)."
)
@click.option(
    "-t", "--p-test", multiple=True, required=False, show_default=True,
    default=False, help="Groups to tests between in each PERMANOVA subset "
         "(multiple values are possible, e.g. '-d sex -d age_cat')."
)
@click.option(
    "-g", "--p-sample-subsets", required=False, show_default=True, default=False,
    help="Subsets for DMs, PCoAs, PERMANOVAs, etc. Must be a yaml file, e.g.\n"
         "(see example in 'examples/permanova_subsets.yml' and README)."
)
@click.option(
    "-nstd", "--p-nestedness", required=False, show_default=True,
    default=False, help="Nestedness analysis config. Must be a yaml file, "
                        "e.g.\n(see example in 'examples/nestedness.yml' "
                        "and README)."
)
@click.option(
    "-bt", "--p-beta-type", required=False, show_default=False,
    multiple=True, default=('permanova', 'permdisp',),
    type=click.Choice(['permanova', 'anosim', 'permdisp']),
    help="Type of beta group significance, one of"
         " 'permanova', 'anosim', 'permdisp'"
)
@click.option(
    "-prc", "--p-procrustes", required=False, show_default=True, default=False,
    help="Pairs and subsets for procrustes/protests. Must be a yaml file, "
         "e.g.\n(see example in 'examples/procrustes_subsets.yml' and README)."
)
@click.option(
    "-mtl", "--p-mantel", required=False, show_default=True, default=False,
    help="Pairs and subsets for mantel test. Must be a yaml file, e.g.\n"
         "(see example in 'examples/procrustes_subsets.yml' and README)."
)
@click.option(
    "-ddecay", "--p-dm-decay", required=False, show_default=True,
    default=False, help="Parameters for distance decay analysis (not "
                        "geographic). Must be a yaml file, e.g.\n(see "
                        "example in 'examples/dm_decay.yml' and README)."
)
@click.option(
    "-gdecay", "--p-geo-decay", required=False, show_default=True,
    default=False, help="Parameters for geographic distance decay analysis. "
                        "Must be a yaml file, e.g.\n(see example in "
                        "'examples/geographic_decay.yml' and README)."
)
@click.option(
    "-coll", "--p-collapse", required=False, show_default=True,
    default=False, help="Nominative or rank-based taxonmic collapse per "
                        "dataset. Must be a yaml file, e.g.\n(see example in "
                        "'examples/collapse.yml' and README)."
)
@click.option(
    "-tt", "--p-train-test", required=False, show_default=True, default=False,
    help="Train test split per dataset. Must be a yaml file, e.g.\n"
         "(see example in 'examples/train_test.yml' and README)."
)
@click.option(
    "-a", "--p-adonis", required=False, default=False,
    show_default=True, help="Formula for Adonis tests for each PERMANOVA "
                            "subset. Must be a yaml file, e.g.\n(see example "
                            "in examples/adonis.yml' and README)."
)
@click.option(
    "-phate", "--p-phate", required=False, default=False,
    show_default=True, help="Filters, subsets, parameters and stratifications "
                            "for the PHATE latent space analysis\nMust be a "
                            "yaml file, (see example in "
                            "'examples/phate.yml' and README)."
)
@click.option(
    "-st", "--p-sourcetracking", required=False, default=False,
    show_default=True, help="Filters, subsets, parameters and isn/sources for "
                            "sourcetracking\nMust be a yaml file, (see example "
                            "in 'examples/sourcetracking.yml' and "
                            "README)."
)
@click.option(
    "-doc", "--p-doc", required=False, default=False, show_default=True,
    help="Filters and subsets for the dissimilarity overlap curves analyses"
         "Must be a yaml file, (see example in 'examples/doc.yml'"
         "and README)."
)
@click.option(
    "-s", "--p-diff-models", required=False, default=False, show_default=True,
    help="Formulas for multinomial regression-based differential "
         "abundance ranking (songbird).\nMUST BE YAML FILE, "
         "see 'examples/songbird_models.yml' and README."
)
@click.option(
    "-m", "--p-mmvec-pairs", required=False, default=False, show_default=True,
    help="Pairs of datasets for which to compute co-occurrences "
         "probabilities (mmvec).\nMUST BE YAML FILE, see "
         "'examples/mmvec_pairs.yml' and README."
)
@click.option(
    "-hlg", "--p-mmvec-highlights", required=False, default=False,
    show_default=True, help="Features to highlights on mmvec biplot (per "
                            "dataset).\nMUST BE YAML FILE, "
                            "see 'examples/mmvec_highlight.yml' and README."
)
@click.option(
    "-mm", "--p-xmmvec", required=False, default=False, show_default=True,
    help="Config for Xmmvec.\n"
         "MUST BE YAML FILE, see 'examples/xmmvec.yml' and README."
)
@click.option(
    "-l", "--p-longi-column", required=False, default=False, show_default=True,
    help="If data is longitudinal; provide the time metadata column"
         "for volatility analysis."
)
@click.option(
    "-f", "--p-filter", show_default=True, default=False,
    help="Minimum sample read abundance to be kept in the sample (per dataset)"
         "(>1 = based on absolute reads, 0-1 = based on relative reads)."
         "Must be a yaml file, e.g. (see example in "
         "'examples/permanova_subsets.yml' and README)."
)
@click.option(
    "-r", "--p-rarefs", required=False, show_default=False,
    help="rarefaction depth per dataset. Must be a yaml file, e.g. "
         "(see example in 'examples/permanova_subsets.yml' and README)."
)
@click.option(
    "-c", "--p-chmod", default='664', show_default=True,
    help="Change output files permission (default = 664 [= -rw-rw-r--])."
)
@click.option(
    "-skip", "--p-skip", default=None, show_default=True, multiple=True,
    type=click.Choice(['alpha', 'merge_alpha', 'export_alpha',
                       'alpha_correlations', 'alpha_group_significance',
                       'wol', 'taxonomy', 'barplot', 'volatility', 'beta',
                       'export_beta', 'pcoa', 'biplot', 'emperor',
                       'emperor_biplot', 'empress', 'empress_biplot',
                       'phate', 'doc', 'deicode', 'sepp', 'do_pies',
                       'alpha_kw', 'permanova', 'procrustes', 'mantel', 'decay',
                       'nestedness', 'adonis', 'songbird', 'mmvec', 'mmbird']),
    help="Steps to skip (e.g. if already done or not necessary)."
         "\nSkipping 'alpha' will also skip 'merge_alpha', 'export_alpha',"
         "'alpha_correlations', 'alpha_kw' and 'volatility'."
         "\nSkipping 'beta' will also skip 'export_beta', 'emperor',"
         "'doc', 'emperor_biplot','deicode', 'permanova', 'adonis',"
         "'procrustes'."
)
@click.option(
    "-As", "--p-alphas", default=None, show_default=True, multiple=True,
    help="Alpha diversity indices to use."
)
@click.option(
    "-Bs", "--p-betas", default=None, show_default=True, multiple=True,
    help="Beta diversity metrics to use."
)
@click.option(
    "--biplot/--no-biplot", default=False, show_default=True,
    help="Whether to do the PCoA biplots or not"
)
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists)."
)
@click.option(
    "--gpu/--no-gpu", default=False, show_default=True,
    help="Use GPUs instead of CPUs for MMVEC."
)
@click.option(
    "--standalone/--no-standalone", default=False, show_default=True,
    help="Whether to run MMVEC using the standalone "
         "version (to check tensorboard)."
)
@click.option(
    "--raref/--no-raref", default=False, show_default=True,
    help="Whether to rarefy and only perform the routine "
         "analyses on the rarefied dataset(s)."
)
@click.option(
    "--loc/--no-loc", default=True, show_default=True,
    help="whether to do compute on scratch."
)
@click.option(
    "--split/--no-split", default=False, show_default=True,
    help="whether to split the commands to a "
         "lower looping level (i.e. more jobs)."
)
@click.option(
    "--filt-only/--no-filt-only", default=False, show_default=True,
    help="Only process the filtered version (and not also the raw) "
         "version of each dataset."
)
@click.option(
    "-filt3d", "--p-filt3d", required=False, show_default=False,
    help="Levels for the exploration of filtering. Must be a yaml file."
)
@click.option(
    "--jobs/--no-jobs", default=True, show_default=True,
    help="Whether to prepare Torque jobs from scripts."
)
@click.option(
    "--slurm/--no-slurm", default=False, show_default=True,
    help="Whether to prepare Slurm and not Torque jobs."
)
@click.option(
    "-chunkit", "--p-chunkit", required=False, show_default=False,
    type=int, default=None,
    help="Maximum number of jobs at which extra jobs will be added in chunks"
)
@click.version_option(__version__, prog_name="routine_qiime2_analyses")


def standalone_analyzer(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        p_qiime2_env,
        p_longi_column,
        p_filter,
        p_rarefs,
        p_feature_subsets,
        p_test,
        p_sample_subsets,
        p_nestedness,
        p_beta_type,
        p_procrustes,
        p_mantel,
        p_dm_decay,
        p_geo_decay,
        p_collapse,
        p_train_test,
        p_adonis,
        p_doc,
        p_sourcetracking,
        p_phate,
        biplot,
        force,
        i_classifier,
        i_wol_tree,
        i_sepp_tree,
        i_qemistree,
        p_diff_models,
        p_mmvec_pairs,
        p_mmvec_highlights,
        p_xmmvec,
        p_run_params,
        p_chmod,
        p_skip,
        gpu,
        standalone,
        raref,
        loc,
        p_alphas,
        p_betas,
        split,
        p_filt3d,
        filt_only,
        jobs,
        slurm,
        p_chunkit
):

    microbiome_analyzer(
        datasets=i_datasets,
        datasets_folder=i_datasets_folder,
        project_name=p_project_name,
        qiime_env=p_qiime2_env,
        longi_column=p_longi_column,
        raref=raref,
        filter_fp=p_filter,
        rarefs_fp=p_rarefs,
        feature_subsets_fp=p_feature_subsets,
        tests=p_test,
        sample_subsets_fp=p_sample_subsets,
        nestedness_fp=p_nestedness,
        beta_type=p_beta_type,
        procrustes_fp=p_procrustes,
        mantel_fp=p_mantel,
        dm_decay_fp=p_dm_decay,
        geo_decay_fp=p_geo_decay,
        collapse_fp=p_collapse,
        train_test_fp=p_train_test,
        adonis_fp=p_adonis,
        doc_fp=p_doc,
        sourcetracking_fp=p_sourcetracking,
        phate_fp=p_phate,
        biplot=biplot,
        force=force,
        classifier=i_classifier,
        wol_tree=i_wol_tree,
        sepp_tree=i_sepp_tree,
        qemistree=i_qemistree,
        diff_models_fp=p_diff_models,
        mmvec_pairs_fp=p_mmvec_pairs,
        mmvec_highlights_fp=p_mmvec_highlights,
        xmmvec_fp=p_xmmvec,
        run_params_fp=p_run_params,
        chmod=p_chmod,
        skip=p_skip,
        gpu=gpu,
        standalone=standalone,
        loc=loc,
        alphas=p_alphas,
        betas=p_betas,
        split=split,
        filt3d_fp=p_filt3d,
        filt_only=filt_only,
        jobs=jobs,
        slurm=slurm,
        chunkit=p_chunkit
    )


if __name__ == "__main__":
    standalone_analyzer()
