# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_taxonomy_classifier,
    get_job_folder,
    get_analysis_folder,
    parse_g2lineage
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_barplots,
    write_seqs_fasta,
    write_taxonomy_sklearn,
    run_export, run_import
)


def run_taxonomy_others(force: bool, tsv_pd: pd.DataFrame,
                        out_qza: str, out_tsv: str) -> str:
    """
    :param force: Force the re-writing of scripts for all commands.
    :param tsv_pd: Current features table for the current dataset.
    :param out_qza: Taxonomy classification output to generate.
    :param out_tsv: Taxonomy classification output exported.
    """
    cmd = ''
    if force or not isfile(out_qza):
        with open(out_tsv, 'w') as o:
            o.write('Feature ID\tTaxon\n')
            for feat in tsv_pd.index:
                o.write('%s\t%s\n' % (feat, feat))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_wol(force: bool, tsv_pd: pd.DataFrame, out_qza: str,
                     out_tsv: str, cur_datasets_features: dict) -> str:
    """
    :param force: Force the re-writing of scripts for all commands.
    :param tsv_pd: Current features table for the current dataset.
    :param out_qza: Taxonomy classification output to generate.
    :param out_tsv: Taxonomy classification output exported.
    :param cur_datasets_features: gotu -> features name containing gotu
    """
    cmd = ''
    if force or not isfile(out_qza):
        g2lineage = parse_g2lineage()
        rev_cur_datasets_features = dict((y, x) for x, y in cur_datasets_features.items())
        with open(out_tsv, 'w') as o:
            o.write('Feature ID\tTaxon\n')
            for feat in tsv_pd.index:
                if rev_cur_datasets_features[feat] in g2lineage:
                    o.write('%s\t%s\n' % (feat, g2lineage[rev_cur_datasets_features[feat]]))
                else:
                    o.write('%s\t%s\n' % (feat, feat.replace('|', '; ')))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_amplicon(dat: str, i_datasets_folder: str, force: bool, tsv_pd: pd.DataFrame,
                          out_qza: str, out_tsv: str, i_classifier: str) -> str:
    """
    :param dat: Current dataset.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param force: Force the re-writing of scripts for all commands.
    :param tsv_pd: Current features table for the current dataset.
    :param out_qza: Taxonomy classification output to generate.
    :param out_tsv: Taxonomy classification output exported.
    :param i_classifier: Path to the taxonomic classifier.
    """
    cmd = ''
    ref_classifier_qza = get_taxonomy_classifier(i_classifier)
    odir_seqs = get_analysis_folder(i_datasets_folder, 'seqs/%s' % dat)
    out_fp_seqs_rad = '%s/seq_%s' % (odir_seqs, dat)
    out_fp_seqs_fasta = '%s.fasta' % out_fp_seqs_rad
    out_fp_seqs_qza = '%s.qza' % out_fp_seqs_rad
    if force or not isfile(out_fp_seqs_qza):
        cmd = write_seqs_fasta(out_fp_seqs_fasta, out_fp_seqs_qza, tsv_pd)
    if force or not isfile(out_qza):
        write_taxonomy_sklearn(out_qza, out_fp_seqs_qza, ref_classifier_qza)
        cmd = run_export(out_qza, out_tsv, '')
    return cmd


def run_taxonomy(i_datasets_folder: str, datasets_read: dict, datasets_phylo: dict,
                 datasets_features: dict, i_classifier: str, taxonomies: dict, force: bool,
                 prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    classify-sklearn: Pre-fitted sklearn-based taxonomy classifier

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param i_classifier: Path to the taxonomic classifier.
    :param taxonomies: dataset -> [method, assignment qza]
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'taxonomy')
    job_folder2 = get_job_folder(i_datasets_folder, 'taxonomy/chunks')

    amplicon_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'amplicon']
    wol_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'wol']

    method = 'sklearn'
    # method = 'hybrid-vsearch-sklearn'
    # method = 'consensus-blast'
    # method = 'consensus-vsearch'
    written = 0
    run_pbs = '%s/1_run_taxonomy.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets_read.items():
            if dat in taxonomies:
                continue
            tsv_pd, meta_pd = tsv_meta_pds
            out_sh = '%s/run_taxonomy_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                odir = get_analysis_folder(i_datasets_folder, 'taxonomy/%s' % dat)
                out_rad = '%s/tax_%s' % (odir, dat)
                if dat in amplicon_datasets:
                    out_qza = '%s/tax_%s.qza' % (out_rad, method)
                    out_tsv = '%s.tsv' % splitext(out_qza)[0]
                    taxonomies[dat] = [method, out_qza]
                    cmd = run_taxonomy_amplicon(dat, i_datasets_folder, force, tsv_pd,
                                                out_qza, out_tsv, i_classifier)
                else:
                    out_qza = '%s.qza' % out_rad
                    out_tsv = '%s.tsv' % out_rad
                    if dat in wol_datasets:
                        cur_datasets_features = datasets_features[dat]
                        taxonomies[dat] = ['wol', out_qza]
                        cmd = run_taxonomy_wol(force, tsv_pd, out_qza, out_tsv,
                                               cur_datasets_features)
                    else:
                        taxonomies[dat] = ['feat', out_qza]
                        cmd = run_taxonomy_others(force, tsv_pd, out_qza, out_tsv)
                if cmd:
                    cur_sh.write('echo "%s"\n' % cmd)
                    cur_sh.write('%s\n\n' % cmd)
                    written += 1
            run_xpbs(out_sh, out_pbs, '%s.tx.sklrn.%s' % (prjct_nm, dat),
                     qiime_env, '4', '1', '4', '200', 'mb',
                     chmod, written, 'single', o)
    if written:
        print('# Classify features using classify-sklearn')
        print('[TO RUN] sh', run_pbs)


def run_barplot(i_datasets_folder: str, datasets: dict, taxonomy: dict,
                force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> None:
    """
    barplot: Visualize taxonomy with an interactive bar plot

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param taxonomy: dataset -> [classification_method, tax_qza]
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'barplot')
    job_folder2 = get_job_folder(i_datasets_folder, 'barplot/chunks')

    written = 0
    run_pbs = '%s/1_run_barplot.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            method, tax_qza = taxonomy[dat]
            qza = '%s.qza' % splitext(tsv)[0]
            out_sh = '%s/run_barplot_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                odir = get_analysis_folder(i_datasets_folder, 'barplot/%s' % dat)
                out_qzv = '%s/bar_%s_%s.qzv' % (odir, dat, method)
                if force or not isfile(out_qzv):
                    write_barplots(out_qzv, qza, meta, tax_qza, cur_sh)
                    written += 1
            run_xpbs(out_sh, out_pbs, '%s.brplt.%s' % (prjct_nm, dat),
                     qiime_env, '4', '1', '1', '200', 'mb',
                     chmod, written, 'single', o)
    if written:
        print('# Make sample compositions barplots')
        print('[TO RUN] sh', run_pbs)
