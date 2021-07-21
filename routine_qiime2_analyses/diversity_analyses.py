# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import pandas as pd
from os.path import isfile, splitext
from skbio.tree import TreeNode
from routine_qiime2_analyses.analyses_prep import AnalysisPrep

from routine_qiime2_analyses._routine_q2_io_utils import (
    get_metrics, read_yaml_file, get_job_folder, get_analysis_folder,
    write_main_sh, get_main_cases_dict, read_meta_pd,
    get_raref_tab_meta_pds, simple_chunks
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict
from routine_qiime2_analyses._routine_q2_cmds import (
    get_case, write_alpha_group_significance_cmd,
    get_new_meta_pd, get_new_alpha_div, write_metadata_tabulate,
    write_diversity_alpha, write_diversity_alpha_correlation,
    write_diversity_alpha_rarefaction,
    write_longitudinal_volatility, get_subset,
    # write_longitudinal_volatility, get_metric, get_subset,
    write_filter_features, run_export
)


class Diversity(object):
    """
    """

    def __init__(self, config, project) -> None:
        self.config = config
        self.project = project
        self.cmds = {}
        self.alpha_metrics = get_metrics('alpha_metrics', config.alphas)
        self.beta_metrics = get_metrics('beta_metrics', config.betas)
        self.alpha_subsets = read_yaml_file(self.config.alpha_subsets)
        self.alphas = {}
        self.betas = {}

    def alpha(self):
        for d, dat in project.datasets.items():
            for idx, tsv in enumerate(dat.tsv):
                qza = dat.qza[idx]
                meta = dat.meta[idx]
                tsv_pd = dat.data[idx]
                meta_pd = dat.metadata[idx]
                cur_raref = dat.rarefs[idx]
            #
            # divs = {}
            # for metric in alpha_metrics:
            #     odir = get_analysis_folder(i_datasets_folder,
            #                                'alpha/%s%s' % (dat, cur_raref))
            #     out_fp = '%s/%s_%s.qza' % (
            #     odir, basename(splitext(qza)[0]), metric)
            #     out_tsv = '%s.tsv' % splitext(out_fp)[0]
            #     if force or not isfile(out_fp):
            #         ret_continue = write_diversity_alpha(out_fp, datasets_phylo,
            #                                              trees,
            #                                              dat, qza, metric,
            #                                              cur_sh, qiime_env)
            #         if ret_continue:
            #             continue
            #         cmd = run_export(out_fp, out_tsv, '')
            #         cur_sh.write('echo "%s"\n' % cmd)
            #         cur_sh.write('%s\n\n' % cmd)
            #         written += 1
            #         main_written += 1
            #     divs.setdefault('', []).append((out_fp, metric))
            #
            # if alpha_subsets and dat in alpha_subsets:
            #     for subset, subset_regex in alpha_subsets[dat].items():
            #         odir = get_analysis_folder(i_datasets_folder,
            #                                    'alpha/%s%s/%s' % (
            #                                    dat, cur_raref, subset))
            #         if dropout:
            #             qza_subset_ = '%s/%s_%s.qza' % (
            #             odir, basename(splitext(qza)[0]), subset)
            #         else:
            #             qza_subset_ = '%s/%s_%s_noDropout.qza' % (
            #             odir, basename(splitext(qza)[0]), subset)
            #         feats_subset = '%s.meta' % splitext(qza_subset_)[0]
            #         feats = get_subset(tsv_pd, subset_regex)
            #         if not len(feats):
            #             continue
            #         subset_pd = pd.DataFrame(
            #             {'Feature ID': feats, 'Subset': [subset] * len(feats)})
            #         subset_pd.to_csv(feats_subset, index=False, sep='\t')
            #         write_filter_features(tsv_pd, feats, qza, qza_subset_,
            #                               feats_subset, cur_sh, dropout)
            #         for metric in alpha_metrics:
            #
            #             if metric in ['faith_pd'] and datasets_phylo[dat][
            #                 1] and dat in trees:
            #                 tree_in_qza = trees[dat][0]
            #                 tree_in_tsv = '%s.tsv' % splitext(tree_in_qza)[0]
            #                 if dropout:
            #                     qza_subset = '%s/%s_%s.qza' % (
            #                     odir, basename(splitext(tree_in_qza)[0]),
            #                     subset)
            #                 else:
            #                     qza_subset = '%s/%s_%s_noDropout.qza' % (
            #                     odir, basename(splitext(tree_in_qza)[0]),
            #                     subset)
            #                 write_filter_features(
            #                     pd.read_csv(tree_in_tsv, header=0, index_col=0,
            #                                 sep='\t'),
            #                     feats, tree_in_qza, qza_subset,
            #                     feats_subset, cur_sh, dropout)
            #             else:
            #                 qza_subset = qza_subset_
            #
            #             out_fp = '%s/%s__%s.qza' % (
            #             odir, basename(splitext(qza_subset)[0]), metric)
            #             out_tsv = '%s.tsv' % splitext(out_fp)[0]
            #
            #             if force or not isfile(out_fp):
            #                 ret_continue = write_diversity_alpha(out_fp,
            #                                                      {dat: [1, 0]},
            #                                                      trees,
            #                                                      dat,
            #                                                      qza_subset,
            #                                                      metric,
            #                                                      cur_sh,
            #                                                      qiime_env)
            #                 if ret_continue:
            #                     continue
            #                 cmd = run_export(out_fp, out_tsv, '')
            #                 cur_sh.write('echo "%s"\n' % cmd)
            #                 cur_sh.write('%s\n\n' % cmd)
            #                 written += 1
            #                 main_written += 1
            #             divs.setdefault(subset, []).append((out_fp, metric))
            # diversities[dat].append(divs)

    def register_command(self, analysis):
        AnalysisPrep.analyses_commands[analysis] = self.cmds[analysis]
