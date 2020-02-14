# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
from skbio.tree import TreeNode
from os.path import isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder, get_analysis_folder, run_import, get_wol_tree, get_sepp_tree
)


def run_sepp(i_datasets_folder: str, datasets: dict, datasets_read: dict, datasets_phylo: dict,
             prjct_nm: str, i_sepp_tree: str, trees: dict, force: bool, qiime_env: str) -> None:
    """
    Run SEPP on the datasets composed or 16S deblur sequences (e.g. from redbiom/Qiita).

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param prjct_nm: Short nick name for your project.
    :param i_sepp_tree: database to use for sepp phylogeny reads placement.
    :param trees: to be update with tree to use for a dataset phylogenetic analyses.
    :param force: Force the re-writing of scripts for all commands.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :return:
    """
    # check whether there's dataset(s) that may use the reference tree (i.e. features are DNA sequences)
    sepp_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'amplicon']
    if len(sepp_datasets):

        ref_tree_qza = get_sepp_tree(i_sepp_tree)

        job_folder = get_job_folder(i_datasets_folder, 'sepp')
        job_folder2 = get_job_folder(i_datasets_folder, 'sepp/chunks')

        written = 0
        main_sh = '%s/1_run_sepp.sh' % job_folder
        with open(main_sh, 'w') as main_o:
            for dat in sepp_datasets:
                tsv, meta = datasets[dat]
                qza = '%s.qza' % splitext(tsv)[0]
                if not isfile(qza):
                    print('Need to first import %s to .qza to do reads placement '
                          '(see "# Import tables to qiime2")\nExiting...' % tsv)
                    sys.exit(1)
                qza_in = '%s_inTree.qza' % splitext(tsv)[0]
                qza_out = '%s_notInTree.qza' % splitext(tsv)[0]
                tsv_pd, meta_pd = datasets_read[dat]

                odir_seqs = get_analysis_folder(i_datasets_folder, 'seqs/%s' % dat)
                odir_sepp = get_analysis_folder(i_datasets_folder, 'sepp/%s' % dat)

                out_fp_seqs_rad = '%s/seq_%s' % (odir_seqs, dat)
                out_fp_seqs_fasta = '%s.fasta' % out_fp_seqs_rad
                out_fp_seqs_qza = '%s.qza' % out_fp_seqs_rad

                out_fp_sepp_tree = '%s/tree_%s.qza' % (odir_sepp, dat)
                trees[dat] = (qza_in, out_fp_sepp_tree)
                out_fp_sepp_plac = '%s/plac_%s.qza' % (odir_sepp, dat)

                out_sh = '%s/run_sepp_%s.sh' % (job_folder2, dat)
                out_pbs = '%s.pbs' % splitext(out_sh)[0]
                with open(out_sh, 'w') as sh:

                    if force or not isfile(out_fp_seqs_qza):
                        with open(out_fp_seqs_fasta, 'w') as fas_o:
                            for seq in tsv_pd.index:
                                fas_o.write('>%s\n%s\n' % (seq.strip(), seq.strip()))
                        cmd = run_import(out_fp_seqs_fasta, out_fp_seqs_qza, 'FeatureData[Sequence]')
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1

                    if force or not isfile(out_fp_sepp_tree):
                        cmd = 'qiime fragment-insertion sepp \\ \n'
                        cmd += '--i-representative-sequences %s \\ \n' % out_fp_seqs_qza
                        cmd += '--i-reference-database %s \\ \n' % ref_tree_qza
                        cmd += '--o-tree %s \\ \n' % out_fp_sepp_tree
                        cmd += '--o-placements %s \\ \n' % out_fp_sepp_plac
                        cmd += '--p-threads 24\n'
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)

                        cmd = 'qiime fragment-insertion filter-features \\ \n'
                        cmd += '--i-table %s \\ \n' % qza
                        cmd += '--i-tree %s \\ \n' % out_fp_sepp_tree
                        cmd += '--o-filtered-table %s \\ \n' % qza_in
                        cmd += '--o-removed-table %s\n' % qza_out
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n\n' % cmd)
                        written += 1
                if written:
                    xpbs_call(out_sh, out_pbs, '%s.spp.%s' % (prjct_nm, dat), qiime_env,
                              '100', '2', '12', '100', 'gb')
                    main_o.write('qsub %s\n' % out_pbs)
                else:
                    os.remove(out_sh)
        if written:
            print("# Fragment insertion using SEPP (%s)" % ', '.join(sepp_datasets))
            print('[TO RUN] sh', main_sh)


def shear_tree(i_datasets_folder: str, datasets_phylo: dict, datasets_features: dict, prjct_nm: str,
               i_wol_tree: str, trees: dict, force: bool, qiime_env: str) -> None:
    """
    Get the sub-tree from the Web of Life tree that corresponds to the gOTUs-labeled features.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param prjct_nm: Short nick name for your project.
    :param i_wol_tree: default on barnacle /projects/wol/profiling/dbs/wol/phylogeny/wol_tree.nwk.
    :param trees: to be update with tree to use for a dataset phylogenetic analyses.
    :param force: Force the re-writing of scripts for all commands.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :return:
    """
    # check whether there's dataset(s) that may use the Web of Life tree (i.e. features contain gID)
    wol_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'wol']
    if len(wol_datasets):

        job_folder = get_job_folder(i_datasets_folder, 'import_trees')
        job_folder2 = get_job_folder(i_datasets_folder, 'import_trees/chunks')

        i_wol_tree = get_wol_tree(i_wol_tree)
        wol = TreeNode.read(i_wol_tree)

        written = 0
        main_sh = '%s/0_run_import_trees.sh' % job_folder
        with open(main_sh, 'w') as main_o:
            for dat in wol_datasets:
                cur_datasets_features = datasets_features[dat]
                wol_features_ = wol.shear(cur_datasets_features)
                wol_features = wol_features_.copy()
                # rename the tip per the features names associated with each gID
                for tip in wol_features_.tips():
                    wol_features.find(tip).name = cur_datasets_features[str(tip.name)]

                analysis_folder = get_analysis_folder(i_datasets_folder, 'fasta_phylo/%s' % dat)
                wol_features_fpo = '%s/tree_%s.nwk' % (analysis_folder, dat)
                wol_features_qza = wol_features_fpo.replace('.nwk', '.qza')
                trees[dat] = ('', wol_features_qza)

                run_sh = '%s/run_import_tree_%s.sh' % (job_folder2, dat)
                run_pbs = run_sh.replace('.sh', '.pbs')
                wol_features.write(wol_features_fpo)

                with open(run_sh, 'w') as o:
                    if force or not isfile(wol_features_qza):
                        cmd = run_import(wol_features_fpo, wol_features_qza, "Phylogeny[Rooted]")
                        o.write("echo '%s'\n" % cmd)
                        o.write('%s\n\n' % cmd)
                        written += 1
                if written:
                    xpbs_call(run_sh, run_pbs, prjct_nm, qiime_env,
                              '1', '1', '1', '20', 'mb')
                    main_o.write('qsub %s\n' % run_pbs)
                else:
                    os.remove(run_sh)
        if written:
            print("# Shear Web of Life tree to features' genome IDs (%s)" % ', '.join(wol_datasets))
            print('[TO RUN] sh', main_sh)