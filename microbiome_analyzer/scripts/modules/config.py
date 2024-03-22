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
    "-d", "--datasets", multiple=True, required=True,
    help="Dataset(s) identifier(s). Multiple is possible: e.g. "
         "-d dataset_number_1 and -d dataset_number_2 for "
         "'table.tsv' and table.tsv'")
@click.option(
    "--all-datasets/--no-all-datasets", default=False,
    help="Set configurations once for all datasets")
@click.option(
    "--all-configs/--no-all-configs", default=False,
    help="Go through the preparation of all configurations one by one")
@click.option(
    "--filtering/--no-filtering", default=False,
    help="Specify the sample/feature filtering parameters (for `--filtering`)")
@click.option(
    "--rarefy/--no-rarefy", default=False,
    help="Specify the rarefaction levels (for `--rarefy`)")
@click.option(
    "--feature-subsets/--no-feature-subsets", default=False,
    help="Specify the terms/regex to subset features (for `--feature-subsets`)")
@click.option(
    "--sample-subsets/--no-sample-subsets", default=False,
    help="Specify variables/factors to subset samples (for `--sample-subsets`)")
@click.option(
    "--time-subject/--no-time-subject", default=False,
    help="Specify variables for time and subject (for CTF/longitudinal)")
@click.option(
    "--permanova/--no-permanova", default=False,
    help="Specify variables for statistical comparisons (for `--permanova`)")
@click.option(
    "--adonis/--no-adonis", default=False,
    help="Specify models to make variance partitioning tests (for `--adonis`)")
@click.option(
    "--nestedness/--no-nestedness", default=False,
    help="Specify variables/params for nestedness testing (for `--nestedness`)")
@click.option(
    "--pairs/--no-pairs", default=False,
    help="Specify datasets pairs (for `--procrustes` and `--mantel`)")
@click.option(
    "--dm-decay/--no-dm-decay", default=False,
    help="Specify params for beta distance decay analysis (for `--dm-decay`)")
@click.option(
    "--geo-decay/--no-geo-decay", default=False,
    help="Specify params for beta vs geo distances decay (for `--geo-decay`)")
@click.option(
    "--collapse/--no-collapse", default=False,
    help="Specify taxonomic levels to collapse feature to (for `--collapse`)")
@click.option(
    "--train-test/--no-train-test", default=False,
    help="Specify variables for train/test stratification (for `--train-test`)")
@click.option(
    "--phate/--no-phate", default=False,
    help="Specify params for PHATE dimensionality reduction (for `--phate`)")
@click.option(
    "--sourcetracking/--no-sourcetracking", default=False,
    help="Specify source/sink in sourcetracker2/FEAST (for `--sourcetracking`)")
@click.option(
    "--doc/--no-doc", default=False,
    help="Specify params to fit the Dissimilarity-Overlab Curves (for `--doc`)")
@click.option(
    "--xmmvec/--no-xmmvec", default=False,
    help="Specify variables to plot conditional proba heatmap (for `--xmmvec`)")
@click.option(
    "--mmvec-highlights/--no-mmvec-highlights", default=False,
    help="Specify variables to highlight on biplots (for `--mmvec-highlights`)")
@click.option(
    "--mmvec-pairs/--no-mmvec-pairs", default=False,
    help="Specify MMVEC dataset pairs, subsets, filters (for `--mmvec-pairs`)")
@click.option(
    "--diff-abund/--no-diff-abund", default=False,
    help="Specify formulas, baselines, subsets (for `--diff-abund`)")
@click.option(
    "--filt3d/--no-filt3d", default=False,
    help="Prepare a configuration file (for `--filt3d`)")
@click.version_option(__version__, prog_name="microbiome_analyzer")

def config(
        analysis_folder,
        datasets,
        filtering,
        rarefy,
        feature_subsets,
        sample_subsets,
        time_subject,
        permanova,
        adonis,
        nestedness,
        pairs,
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
        diff_abund,
        filt3d,
        all_datasets,
        all_configs
):
    """Write configuration files for different analyses:"""
    configure(
        analysis_folder=analysis_folder,
        datasets=datasets,
        filtering=filtering,
        rarefactions=rarefy,
        feature_subsets=feature_subsets,
        sample_subsets=sample_subsets,
        time_subject=time_subject,
        permanova=permanova,
        adonis=adonis,
        nestedness=nestedness,
        pairs=pairs,
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
        diff_abund=diff_abund,
        filt3d=filt3d,
        all_datasets=all_datasets,
        all_configs=all_configs
    )
