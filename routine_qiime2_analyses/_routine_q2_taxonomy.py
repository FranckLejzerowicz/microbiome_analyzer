# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import basename, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_taxonomy_classifier,
    get_job_folder,
    get_analysis_folder
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_barplots,
    write_seqs_fasta,
    write_taxonomy_sklearn,
    run_export
)


def run_taxonomy(i_datasets_folder: str, datasets: dict, datasets_read: dict, i_classifier: str,
                 force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> dict:
    """
    classify-sklearn: Pre-fitted sklearn-based taxonomy classifier

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv, meta]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param i_classifier: Path to the taxonomic classifier.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Short nick name for your project.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    job_folder = get_job_folder(i_datasets_folder, 'taxonomy')
    job_folder2 = get_job_folder(i_datasets_folder, 'taxonomy/chunks')

    ref_classifier_qza = get_taxonomy_classifier(i_classifier)

    written = 0
    taxonomies = {}
    run_pbs = '%s/1_run_taxonomy.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            tsv_pd, meta_pd = datasets_read[dat]
            out_sh = '%s/run_taxonomy_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                qza = '%s.qza' % splitext(tsv)[0]
                odir = get_analysis_folder(i_datasets_folder, 'taxonomy/%s' % dat)
                out_qza = '%s/%s.qza' % (odir, basename(splitext(tsv)[0]).replace('tab_', 'tax_'))
                out_tsv = '%s.tsv' % splitext(out_qza)[0]

                odir_seqs = get_analysis_folder(i_datasets_folder, 'seqs/%s' % dat)
                out_fp_seqs_rad = '%s/seq_%s' % (odir_seqs, dat)
                out_fp_seqs_fasta = '%s.fasta' % out_fp_seqs_rad
                out_fp_seqs_qza = '%s.qza' % out_fp_seqs_rad
                if force or not isfile(out_fp_seqs_qza):
                    write_seqs_fasta(out_fp_seqs_fasta, out_fp_seqs_qza, tsv_pd, cur_sh)
                    written += 1

                taxonomies.setdefault(dat, []).append(('sklearn', out_qza))
                if force or not isfile(out_qza):
                    write_taxonomy_sklearn(out_qza, out_fp_seqs_qza, ref_classifier_qza, cur_sh)
                    cmd = run_export(out_qza, out_tsv, '')
                    cur_sh.write('echo "%s"\n' % cmd)
                    cur_sh.write('%s\n\n' % cmd)
                    written += 1

            run_xpbs(out_sh, out_pbs, '%s.tx.sklrn.%s' % (prjct_nm, dat),
                     qiime_env, '4', '1', '4', '200', 'mb',
                     chmod, written, 'single', o)
    if written:
        print('# Classify features using classify-sklearn')
        print('[TO RUN] sh', run_pbs)
    return taxonomies


def run_barplot(i_datasets_folder: str, datasets: dict, taxonomy: dict,
                force: bool, prjct_nm: str, qiime_env: str, chmod: str) -> dict:
    """
    classify-sklearn: Pre-fitted sklearn-based taxonomy classifier

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
                out_qza = '%s/bar_%s_%s.qza' % (odir, dat, method)
                if force or not isfile(out_qza):
                    write_barplots(out_qza, qza, meta, tax_qza, cur_sh)
                    written += 1
            run_xpbs(out_sh, out_pbs, '%s.brplt.%s' % (prjct_nm, dat),
                     qiime_env, '4', '1', '1', '200', 'mb',
                     chmod, written, 'single', o)
    if written:
        print('# Make sample compositions barplots')
        print('[TO RUN] sh', run_pbs)
