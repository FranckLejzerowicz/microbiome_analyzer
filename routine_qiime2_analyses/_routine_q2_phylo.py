# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from skbio.tree import TreeNode
from os.path import isfile

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder, run_import

ROOT = pkg_resources.resource_filename('routine_qiime2_analyses', 'resources')


def shear_tree(datasets_features: dict, i_folder: str, prjct_nm: str,
               i_wol_tree: str, force: bool, qiime_env: str) -> dict:
    """
    Get the sub-tree from the Web of Life tree that corresponds to the gOTUs-labeled features.

    :param datasets_features: dataset -> list of features names in the dataset tsv / biom file.
    :param i_folder: Path to the folder containing the .tsv datasets.
    :param prjct_nm: Short nick name for your project.
    :param i_wol_tree: default on barnacle /projects/wol/profiling/dbs/wol/phylogeny/web_of_life_tree.nwk.
    :param force: Force the re-writing of scripts for all commands.
    :param qiime_env: name of your qiime2 conda environment (e.g. qiime2-2019.10).
    :return:
    """
    print("# Shear Web of Life tree to features' genome IDs")
    wol_trees = {}

    if not i_wol_tree:
        i_wol_tree = '%s/web_of_life_tree.nwk' % ROOT

    wol = TreeNode.read(i_wol_tree)
    for dat in datasets_features:
        wol_features_ = wol.shear(datasets_features[dat])
        wol_features = wol_features_
        for tip in wol_features_.tips():
            wol_features.find(tip).name = datasets_features[dat][str(tip.name)].replace(';', '|').replace(' ', '')

        analysis_folder = get_analysis_folder(i_folder, 'fasta_phylo/%s' % dat)
        wol_features_fpo = '%s/%s.nwk' % (analysis_folder, dat)
        wol_features_qza = wol_features_fpo.replace('.nwk', '.qza')
        wol_trees[dat] = wol_features_qza

        job_folder = get_job_folder(i_folder, 'import_tree_%s' % dat)
        run_sh = '%s/0_import_tree.sh' % job_folder
        run_pbs = run_sh.replace('.sh', '.pbs')
        wol_features.write(wol_features_fpo)

        written = 0
        with open(run_sh, 'w') as o:
            if force or not isfile(wol_features_qza):
                cmd = run_import(wol_features_fpo, wol_features_qza, "Phylogeny[Rooted]")
                o.write("echo '%s'\n" % cmd)
                o.write('%s\n' % cmd)
                written += 1
        if written:
            xpbs_call(run_sh, run_pbs, prjct_nm, qiime_env,
                      '1', '1', '1', '20', 'mb')
            print('[TO RUN] qsub', run_pbs)
    return wol_trees
