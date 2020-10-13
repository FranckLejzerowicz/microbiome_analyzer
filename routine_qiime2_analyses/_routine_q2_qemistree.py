# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_cmds import write_qemistree, run_export
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder
)


def run_qemistree(i_datasets_folder: str, datasets: dict, prjct_nm: str,
                  i_qemistree: str, taxonomies: dict, force: bool,  qiime_env: str,
                  chmod: str, noloc: bool, run_params: dict, filt_raref: str, jobs: bool) -> None:
    """
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets_read: dataset -> [tsv table, meta table]
    :param prjct_nm: Short nick name for your project.
    :param i_qemistree: path to qemistree folder (feature-data and tree).
    :param taxonomies: dataset -> [method, assignment qza]
    :param force: Force the re-writing of scripts for all commands.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """

    job_folder = get_job_folder(i_datasets_folder, 'qemistree')
    job_folder2 = get_job_folder(i_datasets_folder, 'qemistree/chunks')

    written = 0
    run_pbs = '%s/1_run_qemistree%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for dat, tsv_meta_pds in datasets.items():
            feature_data = '%s/feature-data_%s.qza' % (i_qemistree, dat)
            qemistree = '%s/qemistree_%s.qza' % (i_qemistree, dat)
            if not isfile(feature_data) or not isfile(qemistree):
                continue
            out_sh = '%s/run_qemistree_%s%s.sh' % (job_folder2, dat, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            odir = get_analysis_folder(i_datasets_folder, 'qemistree/%s' % dat)
            classyfire_qza = '%s/%s-classyfire.qza' % (odir, dat)
            classyfire_tsv = '%s.tsv' % splitext(classyfire_qza)[0]
            with open(out_sh, 'w') as cur_sh:
                if force or not isfile(classyfire_tsv):
                    write_qemistree(feature_data, classyfire_qza,
                                    classyfire_tsv, qemistree,
                                    cur_sh)
                    written += 1

            if isfile(classyfire_tsv):
                odir = get_analysis_folder(i_datasets_folder, 'taxonomy/%s' % dat)
                out_rad = '%s/tax_%s' % (odir, dat)
                tax_qza = '%s.qza' % out_rad
                tax_tsv = '%s.tsv' % out_rad
                classyfire_pd = pd.read_csv(classyfire_tsv, header=0, sep='\t')
                with open(tax_tsv, 'w') as o:
                    cols = ['id', 'kingdom', 'superclass', 'class', 'subclass', 'direct_parent']
                    o.write('Feature ID\tTaxon\n')
                    for row in classyfire_pd[cols].values:
                        o.write('%s\t%s\n' % (row[0], '; '.join(row[1:])))
                run_export(tax_tsv, tax_qza, 'FeatureData[Taxonomy]')
                taxonomies[dat] = ['direct_parent', tax_qza]
                written += 1
            else:
                print('[Warning] Maybe run qemistree first and then re-run pipeline to '
                      'have the classyfire taxonomy include in the barplots!')

            run_xpbs(out_sh, out_pbs, '%s.qmstr.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc, jobs)
    if written:
        print_message('# Make qemistree classyfire classifications', 'sh', run_pbs, jobs)
