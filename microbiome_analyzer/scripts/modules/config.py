# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from microbiome_analyzer.config import configure
from microbiome_analyzer import __version__


@click.command()
@click.option(
    "-i", "--analysis-folder", required=True,
    help="Path to the folder containing the data and metadata sub-folders")
@click.option(
    "-o", "--configs-folder", required=True,
    help="Path to the folder to contain the config files")
@click.option(
    "-d", "--datasets", multiple=True, required=True,
    help="Dataset(s) identifier(s). Multiple is possible: e.g. "
         "-d dataset_number_1 and -d dataset_number_2 for "
         "'tab_dataset_number_1.tsv' and tab_dataset_number_2.tsv'")
@click.option(
    "--filter/--no-filter", default=False,
    help="Prepare a template configuration file for `--filter`")
@click.option(
    "--rarefy/--no-rarefy", default=False,
    help="Prepare a template configuration file for `--rarefy`")
@click.option(
    "--feature-subsets/--no-feature-subsets", default=False,
    help="Prepare a template configuration file for `--feature-subsets`")
@click.option(
    "--sample-subsets/--no-sample-subsets", default=False,
    help="Prepare a template configuration file for `--sample-subsets`")
@click.option(
    "--adonis/--no-adonis", default=False,
    help="Prepare a template configuration file for `--adonis`")
@click.option(
    "--nestedness/--no-nestedness", default=False,
    help="Prepare a template configuration file for `--nestedness`")
@click.option(
    "--procrustes/--no-procrustes", default=False,
    help="Prepare a template configuration file for `--procrustes`")
@click.option(
    "--mantel/--no-mantel", default=False,
    help="Prepare a template configuration file for `--mantel`")
@click.option(
    "--dm-decay/--no-dm-decay", default=False,
    help="Prepare a template configuration file for `--dm-decay`")
@click.option(
    "--geo-decay/--no-geo-decay", default=False,
    help="Prepare a template configuration file for `--geo-decay`")
@click.option(
    "--collapse/--no-collapse", default=False,
    help="Prepare a template configuration file for `--collapse`")
@click.option(
    "--train-test/--notrain-test-", default=False,
    help="Prepare a template configuration file for `--train-test`")
@click.option(
    "--phate/--no-phate", default=False,
    help="Prepare a template configuration file for `--phate`")
@click.option(
    "--sourcetracking/--no-sourcetracking", default=False,
    help="Prepare a template configuration file for `--sourcetracking`")
@click.option(
    "--doc/--no-doc", default=False,
    help="Prepare a template configuration file for `--doc`")
@click.option(
    "--xmmvec/--no-xmmvec", default=False,
    help="Prepare a template configuration file for `--xmmvec`")
@click.option(
    "--mmvec-highlights/--no-mmvec-highlights", default=False,
    help="Prepare a template configuration file for `--mmvec-highlights`")
@click.option(
    "--mmvec-pairs/--no-mmvec-pairs", default=False,
    help="Prepare a template configuration file for `--mmvec-pairs`")
@click.option(
    "--diff-models/--no-diff-models", default=False,
    help="Prepare a template configuration file for `--diff-models`")
@click.option(
    "--filt3d/--no-filt3d", default=False,
    help="Prepare a template configuration file for `--filt3d`")
@click.option(
    "--all-datasets/--no-all-datasets", default=False,
    help="Setup all template configuration files once for all datasets")
@click.version_option(__version__, prog_name="microbiome_analyzer")

def config(
        analysis_folder,
        configs_folder,
        datasets,
        filter,
        rarefy,
        feature_subsets,
        sample_subsets,
        adonis,
        nestedness,
        procrustes,
        mantel,
        dm_decay,
        geo_decay,
        collapse,
        train_test,
        phate,
        sourcetracking,
        doc,
        xmmvec,
        mmvec_highlights,
        mmvec_pairs,
        diff_models,
        filt3d,
        all_datasets
):
    """Write template config files for different analyses."""
    configure(
        analysis_folder=analysis_folder,
        configs_folder=configs_folder,
        datasets=datasets,
        filtering=filter,
        rarefactions=rarefy,
        feature_subsets=feature_subsets,
        sample_subsets=sample_subsets,
        adonis=adonis,
        nestedness=nestedness,
        procrustes=procrustes,
        mantel=mantel,
        dm_decay=dm_decay,
        geo_decay=geo_decay,
        collapse=collapse,
        train_test=train_test,
        phate=phate,
        sourcetracking=sourcetracking,
        doc=doc,
        xmmvec=xmmvec,
        mmvec_highlights=mmvec_highlights,
        mmvec_pairs=mmvec_pairs,
        diff_models=diff_models,
        filt3d=filt3d,
        all_datasets=all_datasets
    )
