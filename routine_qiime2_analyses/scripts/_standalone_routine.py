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
    "-n", "--p-project-name", required=True, show_default=True,
    help="Nick name for your project."
)
@click.option(
    "-e", "--p-qiime2-env", required=True, show_default=True,
    help="name of your qiime2 conda environment (e.g. qiime2-2019.10) "
)
@click.option(
    "-t", "--p-perm-tests", multiple=True, required=False, show_default=True, default=False,
    help="Groups to tests between in each PERMANOVA subset "
         "(multiple values are possible, e.g. '-d sex -d age_cat')."
)
@click.option(
    "-g", "--p-perm-groups", required=False, show_default=True, default=False,
    help="Subsets for PERMANOVA. Must be a yaml file, e.g.\n"
         "(see example in 'examples/example_PERMANOVA_subsets.yml' and the README)."
)
@click.option(
    "-a", "--p-adonis-formulas", required=False, default=False, show_default=True,
    help="Formula for Adonis tests for each PERMANOVA subset. Must be a yaml file, e.g.\n"
         "(see example in 'examples/example_ADONIS_formulas.yml' and the README)."
)
@click.option(
    "-l", "--p-longi-column", required=False, default=False, show_default=True,
    help="If data is longitudinal; provide the time metadata column"
         "for volatility analysis."
)
@click.option(
    "-f", "--p-reads-filter", default=0, show_default=True, type=int,
    help="Minimum number of reads per sample to be kept."
)
@click.option(
    "-skip", "--p-skip", default=None, show_default=True, multiple=True,
    type=click.Choice(['alpha', 'merge_alpha', 'export_alpha', 'alpha_correlations',
                       'volatility', 'beta', 'export_beta', 'emperor', 'deicode',
                       'alpha_kw', 'permanova', 'adonis']),
    help="Steps to skip (e.g. if already done or not necessary)."
         "\nSkipping 'alpha' will also skip 'merge_alpha', 'export_alpha',"
         "'alpha_correlations', 'alpha_kw' and 'volatility'."
         "\nSkipping 'beta' will also skip 'export_beta', 'emperor',"
         "'deicode', 'permanova', 'adonis''."
)
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists)."
)
@click.version_option(__version__, prog_name="routine_qiime2_analyses")


def standalone_routine(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        p_longi_column,
        p_reads_filter,
        p_perm_tests,
        p_perm_groups,
        p_adonis_formulas,
        force, i_wol_tree,
        i_sepp_tree,
        p_qiime2_env,
        p_skip
):

    routine_qiime2_analyses(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        p_longi_column,
        p_reads_filter,
        p_perm_tests,
        p_perm_groups,
        p_adonis_formulas,
        force, i_wol_tree,
        i_sepp_tree,
        p_qiime2_env,
        p_skip
    )


if __name__ == "__main__":
    standalone_routine()
