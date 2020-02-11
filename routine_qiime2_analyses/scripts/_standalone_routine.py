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
    help="Identifier(s) of the dataset(s) (e.g. '-d test1 -d test2' for inputs"
         "'data/tab_test1.tsv + metadata/meta_test1.tsv' as"
         " well as 'data/tab_test2.tsv + metadata/meta_test2.tsv')"
)
@click.option(
    "-t", "--i-wol-tree", required=False, show_default=True, default = None,
    help="default on barnacle /projects/wol/profiling/dbs/wol/phylogeny/tree.nwk"
)
@click.option(
    "-n", "--p-project-name", required=True, show_default=True,
    help="Nick name for your project."
)
@click.option(
    "-e", "--p-qiime2-env", required=True, show_default=True,
    help="name of your qiime2 conda environment (e.g. qiime2-2019.10)."
)
@click.option(
    "-s", "--p-perm-subsets", multiple=True, required=False, show_default=True, default=False,
    help="Groups to tests between in each PERMANOVA subset."
)
@click.option(
    "-g", "--p-perm-groups", required=False, show_default=True,
    default=False, help="Subsets for PERMANOVA. Must be a yaml file, e.g.\n"
         "(see example in 'example_PERMANOVA.yml' and the README)."
)
@click.option(
    "-a", "--p-adonis-formulas", default=False, show_default=True,
    help="Formula for Adonis tests for each PERMANOVA subset. Must be a yaml file, e.g.\n"
         "(see example in 'example_ADONIS.yml' and the README)."
)
@click.option(
    "-l", "--p-longi-column", default=False, show_default=True,
    help="If data is longitudinal; provide the time metadata column"
         "for volatility analysis."
)
@click.option(
    "-f", "--p-reads-filter", default=0, show_default=True, type=int,
    help="Minimum number of reads per sample to be kept."
)
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists)."
)
@click.option(
    "--gid/--no-gid", default=False, show_default=True,
    help="If feature names have the genome ID (to use the Web of Life tree)."
)
@click.option(
    "--biom/--no-biom", default=False, show_default=True,
    help="Use biom files in the input folder"
         "(automatic if there's no .tsv but only .biom file(s))."
)
@click.version_option(__version__, prog_name="routine_qiime2_analyses")


def standalone_routine(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        gid, p_longi_column,
        p_reads_filter,
        p_perm_subsets,
        p_perm_groups,
        p_adonis_formulas,
        force, i_wol_tree,
        p_qiime2_env, biom
):

    routine_qiime2_analyses(
        i_datasets,
        i_datasets_folder,
        p_project_name,
        gid, p_longi_column,
        p_reads_filter,
        p_perm_subsets,
        p_perm_groups,
        p_adonis_formulas,
        force, i_wol_tree,
        p_qiime2_env, biom
    )


if __name__ == "__main__":
    standalone_routine()
