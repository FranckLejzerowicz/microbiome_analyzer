# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from routine_qiime2_analyses._routine_qiime2_analyses import routine_qiime2_analyses
from routine_qiime2_analyses import __version__


@click.command()
@click.option(
    "-i", "--i-datasets-folder", required=True,
    help="Path to the folder containing the sub-folders 'data' and 'metadata'."
)
@click.option(
    "-d", "--i-datasets", multiple=True, required=True,
    help="Dataset(s) identifier(s). Multiple is possible: e.g. -d dataset_number_1 and "
         "-d dataset_number_2 for 'tab_dataset_number_1.tsv' and tab_dataset_number_2.tsv'."
)
@click.option(
    "-w", "--i-wol-tree", required=False, show_default=True, default = 'resources/wol_tree.nwk',
    help="path to the tree containing the genome IDs (will check if exist in features names)"
         "(On barnacle, it is there: /projects/wol/profiling/dbs/wol/phylogeny/tree.nwk)."
)
@click.option(
    "-x", "--i-sepp-tree", required=False, show_default=True, default = None,
    help="Qiime2 SEPP reference database to use for 16S reads placement: "
         "https://docs.qiime2.org/2019.10/data-resources/#sepp-reference-databases "
         "(auto detection of datasets' tables with sequences as features)."
)
@click.option(
    "-q", "--i-qemistree", required=False, show_default=True, default = None,
    help="Path to a folder containing Qemistree's feature data "
         "(named 'feature-data_<dataset_identifier>.qza'), "
         "and tree for each metabolomics dataset "
         "(named 'qemistree_<dataset_identifier>.qza')."
)
@click.option(
    "-z", "--i-classifier", required=False, show_default=True, default = None,
    help="Qiime2 reference taxonomic classifier database to use for 16S reads assignment: "
         "https://docs.qiime2.org/2020.2/data-resources/#taxonomy-classifiers-for-use-with-q2-feature-classifier"
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
    "-u", "--p-run-params", required=False, show_default=True,
    help="server run paramters"
)
@click.option(
    "-b", "--p-alpha-subsets", required=False, show_default=True, default=None,
    help="Regex to use for subsetting features for alpha diversity subsets (yml file)."
)
@click.option(
    "-k", "--p-beta-subsets", required=False, show_default=True, default=None,
    help="Regex to use for subsetting features for beta diversity subsets (yml file)."
)
@click.option(
    "-t", "--p-perm-tests", multiple=True, required=False, show_default=True, default=False,
    help="Groups to tests between in each PERMANOVA subset "
         "(multiple values are possible, e.g. '-d sex -d age_cat')."
)
@click.option(
    "-g", "--p-perm-groups", required=False, show_default=True, default=False,
    help="Subsets for PERMANOVA. Must be a yaml file, e.g.\n"
         "(see example in 'examples/permanova_subsets.yml' and README)."
)
@click.option(
    "-bt", "--p-beta-type", required=False, show_default=False, multiple=True, default=('permanova',),
    type=click.Choice(['permanova', 'anosim', 'permdisp']),
    help="Type of beta group significance, one of"
         " 'permanova', 'anosim', 'permdisp'"
)
@click.option(
    "-prc", "--p-procrustes", required=False, show_default=True, default=False,
    help="Pairs and subsets for procrustes/protests. Must be a yaml file, e.g.\n"
         "(see example in 'examples/procrustes_subsets.yml' and README)."
)
@click.option(
    "-a", "--p-adonis-formulas", required=False, default=False, show_default=True,
    help="Formula for Adonis tests for each PERMANOVA subset. Must be a yaml file, e.g.\n"
         "(see example in 'examples/adonis_formulas.yml' and README)."
)
@click.option(
    "-s", "--p-diff-models", required=False, default=False, show_default=True,
    help="Formulas for multinomial regression-based differential abundance ranking (songbird).\n"
         "MUST BE YAML FILE, see 'examples/songbird_models.yml' and README."
)
@click.option(
    "-m", "--p-mmvec-pairs", required=False, default=False, show_default=True,
    help="Pairs of datasets for which to compute co-occurrences probabilities (mmvec).\n"
         "MUST BE YAML FILE, see 'examples/mmvec_pairs.yml' and README."
)
@click.option(
    "-h", "--p-mmvec-highlights", required=False, default=False, show_default=True,
    help="Features to highlights on mmvec biplot (per dataset).\n"
         "MUST BE YAML FILE, see 'examples/mmvec_highlights.yml' and README."
)
@click.option(
    "-l", "--p-longi-column", required=False, default=False, show_default=True,
    help="If data is longitudinal; provide the time metadata column"
         "for volatility analysis."
)
@click.option(
    "-f", "--p-filt-threshs", show_default=True, default=False,
    help="Minimum reads abundance per sample to be kept in the sample / per dataset"
         "(>1 = based on absolute reads, 0-1 = based on relative reads)."
         "Must be a yaml file, e.g. (see example in 'examples/permanova_subsets.yml' and README)."
)
@click.option(
    "-r", "--p-raref-depths", required=False, show_default=False,
    help="rarefaction depth per dataset. Must be a yaml file, e.g. "
         "(see example in 'examples/permanova_subsets.yml' and README)."
)
@click.option(
    "-c", "--p-chmod", default='664', show_default=True,
    help="Change output files permission (default = 664 [= -rw-rw-r--])."
)
@click.option(
    "--eval-rarefs/--no-eval-rarefs", required=False, show_default=False,
    help="Perform an evaluation of the various rarefaction depths "
         "and the passed one."
)
@click.option(
    "-skip", "--p-skip", default=None, show_default=True, multiple=True,
    type=click.Choice(['alpha', 'merge_alpha', 'export_alpha', 'alpha_correlations',
                       'wol', 'taxonomy', 'barplot', 'volatility', 'beta',
                       'export_beta', 'emperor', 'emperor_biplot', 'deicode', 'sepp', 'alpha_kw',
                       'permanova', 'adonis', 'songbird', 'mmvec', 'mmbird']),
    help="Steps to skip (e.g. if already done or not necessary)."
         "\nSkipping 'alpha' will also skip 'merge_alpha', 'export_alpha',"
         "'alpha_correlations', 'alpha_kw' and 'volatility'."
         "\nSkipping 'beta' will also skip 'export_beta', 'emperor',"
         "'emperor_biplot','deicode', 'permanova', 'adonis'."
)
@click.option(
    "-As", "--p-alphas", default=None, show_default=True, multiple=True,
    type=click.Choice(['observed_otus','pielou_e','shannon','faith_pd']),
    help="Alpha to use: 'observed_otus','pielou_e','shannon','faith_pd'."
)
@click.option(
    "-Bs", "--p-betas", default=None, show_default=True, multiple=True,
    type=click.Choice(['jaccard','braycurtis','aitchison','unweighted_unifrac','weighted_unifrac']),
    help="Beta to use: 'jaccard','braycurtis','aitchison','unweighted_unifrac','weighted_unifrac'."
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
    help="Whether to run MMVEC using the standalone version (to check tensorboard)."
)
@click.option(
    "--raref/--no-raref", default=False, show_default=True,
    help="Rarefy and only perform the routine analyses on the rarefied dataset(s)."
)
@click.option(
    "--loc/--no-loc", default=True, show_default=True,
    help="whether to do compute on scratch."
)
@click.option(
    "--split/--no-split", default=False, show_default=True,
    help="whether to split the commands to a lower looping level (i.e. more jobs)."
)
@click.option(
    "--dropout/--no-dropout", default=True, show_default=True,
    help="Allow samples dropouts for the features subsetting (alpha and beta)."
)
@click.option(
    "--filt3d/--no-filt3d", default=False, show_default=True,
    help="Only check the effect of filtering thresholds passed to songbird "
         "and mmvec instead of running all of them."
)
@click.option(
    "-filt3d_config", "--p-filt3d-config", required=False, show_default=False,
    help="Levels for the exploration of filtering. Must be a yaml file."
)
@click.version_option(__version__, prog_name="routine_qiime2_analyses")


def standalone_routine(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        p_longi_column,
        p_filt_threshs,
        p_raref_depths,
        eval_rarefs,
        p_alpha_subsets,
        p_beta_subsets,
        p_perm_tests,
        p_perm_groups,
        p_beta_type,
        p_procrustes,
        p_adonis_formulas,
        force,
        i_classifier,
        i_wol_tree,
        i_sepp_tree,
        i_qemistree,
        p_diff_models,
        p_mmvec_pairs,
        p_mmvec_highlights,
        p_qiime2_env,
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
        dropout,
        filt3d,
        p_filt3d_config
):

    routine_qiime2_analyses(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        p_longi_column,
        p_filt_threshs,
        p_raref_depths,
        eval_rarefs,
        p_alpha_subsets,
        p_beta_subsets,
        p_perm_tests,
        p_perm_groups,
        p_beta_type,
        p_procrustes,
        p_adonis_formulas,
        force,
        i_classifier,
        i_wol_tree,
        i_sepp_tree,
        i_qemistree,
        p_diff_models,
        p_mmvec_pairs,
        p_mmvec_highlights,
        p_qiime2_env,
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
        dropout,
        filt3d,
        p_filt3d_config
    )


if __name__ == "__main__":
    standalone_routine()
