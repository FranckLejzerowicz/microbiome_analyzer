# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from skbio.tree import TreeNode
from os.path import isfile, splitext
import pandas as pd

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_wol_tree,
    get_sepp_tree,
    get_raref_tab_meta_pds
)
from routine_qiime2_analyses._routine_q2_cmds import (
    write_fragment_insertion, write_seqs_fasta)
from routine_qiime2_analyses._routine_q2_cmds import run_import


def run_sepp(i_datasets_folder: str, datasets: dict, datasets_read: dict,
             datasets_phylo: dict, datasets_rarefs: dict, prjct_nm: str,
             i_sepp_tree: str, trees: dict, force: bool, qiime_env: str,
             chmod: str, noloc: bool, run_params: dict, filt_raref: str,
             jobs: bool) -> None:
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
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    # check whether there's dataset(s) that may use the reference tree (i.e. features are DNA sequences)
    sepp_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'amplicon']
    if len(sepp_datasets):
        ref_tree_qza = get_sepp_tree(i_sepp_tree)

        job_folder = get_job_folder(i_datasets_folder, 'phylo')
        job_folder2 = get_job_folder(i_datasets_folder, 'phylo/chunks')
        main_written = 0
        main_sh = '%s/1_run_sepp_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
        with open(main_sh, 'w') as main_o:
            for dat, tsv_metas_fps_ in datasets.items():
                written = 0
                if dat not in sepp_datasets:
                    continue
                out_sh = '%s/run_sepp_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
                out_pbs = '%s.pbs' % splitext(out_sh)[0]
                with open(out_sh, 'w') as cur_sh:
                    for idx, tsv_metas_fps in enumerate(tsv_metas_fps_):
                        tsv, meta = tsv_metas_fps
                        if not isinstance(datasets_read[dat][idx][0], pd.DataFrame) and datasets_read[dat][idx][0] == 'raref':
                            qza_raw_in = '%s/data/tab_%s_inTree.qza' % (i_datasets_folder, dat)
                            if isfile(qza_raw_in) and not force:
                                odir_sepp = get_analysis_folder(i_datasets_folder, 'phylo/%s' % dat)
                                out_fp_sepp_tree = '%s/tree_%s.qza' % (odir_sepp, dat)
                                # if idx:
                                #     trees[dat].append((qza_raw_in, out_fp_sepp_tree))
                                # else:
                                #     trees[dat] = [(qza_raw_in, out_fp_sepp_tree)]
                                if not idx:
                                    trees[dat] = (qza_raw_in, out_fp_sepp_tree)
                                print('Using the non rarefied tree (no need to recompute)...\nExiting')
                                continue
                            elif not isfile(tsv):
                                print('Must have run rarefaction to use it further...\nExiting')
                                sys.exit(0)
                            tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                            datasets_read[dat][idx] = [tsv_pd, meta_pd]
                        else:
                            tsv_pd, meta_pd = datasets_read[dat][idx]

                        qza = '%s.qza' % splitext(tsv)[0]
                        if not isfile(qza):
                            print('Need to first import %s to .qza to do reads placement '
                                  '(see "# Import tables to qiime2")\nExiting...' % tsv)
                            sys.exit(0)

                        cur_raref = datasets_rarefs[dat][idx]
                        qza_in = '%s_inTree%s.qza' % (splitext(tsv)[0], cur_raref)
                        qza_out = '%s_notInTree%s.qza' % (splitext(tsv)[0], cur_raref)

                        odir_seqs = get_analysis_folder(i_datasets_folder, 'seqs/%s' % dat)
                        odir_sepp = get_analysis_folder(i_datasets_folder, 'phylo/%s' % dat)

                        out_fp_seqs_rad = '%s/seq_%s%s' % (odir_seqs, dat, cur_raref)
                        out_fp_seqs_fasta = '%s.fasta' % out_fp_seqs_rad
                        out_fp_seqs_qza = '%s.qza' % out_fp_seqs_rad

                        out_fp_sepp_tree = '%s/tree_%s%s.qza' % (odir_sepp, dat, cur_raref)
                        # if idx:
                        #     trees[dat].append((qza_in, out_fp_sepp_tree))
                        # else:
                        #     trees[dat] = [(qza_in, out_fp_sepp_tree)]
                        if not idx:
                            trees[dat] = (qza_in, out_fp_sepp_tree)

                        written = 0
                        if force or not isfile(out_fp_seqs_qza):
                            cmd = write_seqs_fasta(out_fp_seqs_fasta,
                                                   out_fp_seqs_qza, tsv_pd)
                            cur_sh.write('echo "%s"\n' % cmd)
                            cur_sh.write('%s\n\n' % cmd)
                            written += 1
                        if force or not isfile(out_fp_sepp_tree):
                            cmd = write_fragment_insertion(
                                out_fp_seqs_qza, ref_tree_qza,
                                out_fp_sepp_tree, qza,
                                qza_in)
                            cur_sh.write('echo "%s"\n' % cmd)
                            cur_sh.write('%s\n\n' % cmd)
                            written += 1
                            main_written += 1
                run_xpbs(out_sh, out_pbs, '%s.spp.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                         run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                         run_params["mem_num"], run_params["mem_dim"],
                         chmod, written, 'single', main_o, noloc, jobs)
        if main_written:
            print_message("# Fragment insertion using SEPP (%s)" % ', '.join(sepp_datasets), 'sh', main_sh, jobs)


def shear_tree(
        i_datasets_folder: str, datasets: dict, datasets_read: dict,
        datasets_phylo: dict, datasets_features: dict, prjct_nm: str,
        i_wol_tree: str, trees: dict, datasets_rarefs: dict, force: bool,
        qiime_env: str, chmod: str, noloc: bool, run_params: dict,
        filt_raref: str, jobs: bool) -> None:
    """
    Get the sub-tree from the Web of Life tree that corresponds to the gOTUs-labeled features.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets_read: dataset -> [tsv table, meta table]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param prjct_nm: Short nick name for your project.
    :param i_wol_tree: default on barnacle /projects/wol/profiling/dbs/wol/phylogeny/wol_tree.nwk.
    :param trees: to be update with tree to use for a dataset phylogenetic analyses.
    :param force: Force the re-writing of scripts for all commands.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    # check whether there's dataset(s) that may use the Web of Life tree (i.e. features contain gID)
    wol_datasets = [dat for dat, (tree, correction) in datasets_phylo.items() if tree == 'wol']
    if len(wol_datasets):
        job_folder = get_job_folder(i_datasets_folder, 'phylo')
        job_folder2 = get_job_folder(i_datasets_folder, 'phylo/chunks')

        i_wol_tree = get_wol_tree(i_wol_tree)
        wol = TreeNode.read(i_wol_tree)

        main_written = 0
        main_sh = '%s/0_run_import_trees_%s%s.sh' % (job_folder, prjct_nm, filt_raref)
        with open(main_sh, 'w') as main_o:
            for dat, tsv_metas_fps_ in datasets.items():
                written = 0
                if dat not in wol_datasets:
                    continue
                out_sh = '%s/run_import_tree_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
                out_pbs = out_sh.replace('.sh', '.pbs')
                with open(out_sh, 'w') as o:
                    for idx, tsv_metas_fps in enumerate(tsv_metas_fps_):
                        tsv, meta = tsv_metas_fps
                        if not isinstance(datasets_read[dat][idx][0], pd.DataFrame) and datasets_read[dat][idx][0] == 'raref':
                            if not isfile(tsv):
                                print('Must have run rarefaction to use it further...\nExiting')
                                sys.exit(0)
                            tsv_pd, meta_pd = get_raref_tab_meta_pds(meta, tsv)
                            datasets_read[dat][idx] = [tsv_pd, meta_pd]
                        else:
                            tsv_pd, meta_pd = datasets_read[dat][idx]
                        cur_raref = datasets_rarefs[dat][idx]
                        cur_datasets_features = dict(
                            gid for gid in datasets_features[dat].items() if gid[1] in tsv_pd.index)

                        analysis_folder = get_analysis_folder(i_datasets_folder, 'phylo/%s' % dat)
                        wol_features_fpo = '%s/tree_%s%s.nwk' % (analysis_folder, dat, cur_raref)
                        wol_features_qza = wol_features_fpo.replace('.nwk', '.qza')

                        # if idx:
                        #     trees[dat].append(('', wol_features_qza))
                        # else:
                        #     trees[dat] = [('', wol_features_qza)]
                        if not idx:
                            trees[dat] = ('', wol_features_qza)

                        if force or not isfile(wol_features_qza):
                            wol_features = wol.shear(list(cur_datasets_features.keys()))
                            # rename the tip per the features names associated with each gID
                            for tip in wol_features.tips():
                                tip.name = cur_datasets_features[tip.name]
                            wol_features.write(wol_features_fpo)
                            cmd = run_import(wol_features_fpo, wol_features_qza, "Phylogeny[Rooted]")
                            o.write("echo '%s'\n" % cmd)
                            o.write('%s\n\n' % cmd)
                        written += 1
                        main_written += 1
                run_xpbs(out_sh, out_pbs, '%s.shr.%s%s' % (prjct_nm, dat, filt_raref), qiime_env,
                         run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                         run_params["mem_num"], run_params["mem_dim"],
                         chmod, written, 'single', main_o, noloc, jobs)
        if main_written:
            print_message("# Shear Web of Life tree to features' genome IDs (%s)" % ', '.join(wol_datasets), 'sh', main_sh, jobs)


def get_precomputed_trees(
        i_datasets_folder: str, datasets: dict, datasets_filt_map: dict,
        datasets_phylo: dict, trees: dict) -> None:
    """
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: dataset -> [tsv/biom path, meta path]
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    :param trees: to be update with tree to use for a dataset phylogenetic analyses.
    """
    for dat in datasets:
        if dat in datasets_filt_map:
            dat_tree = datasets_filt_map[dat]
        else:
            dat_tree = dat
        analysis_folder = get_analysis_folder(i_datasets_folder, 'phylo/%s' % dat_tree)
        tree_qza = '%s/tree_%s.qza' % (analysis_folder, dat_tree)
        if isfile(tree_qza):
            trees[dat] = ('', tree_qza)
            datasets_phylo[dat] = ('precpu', 0)
