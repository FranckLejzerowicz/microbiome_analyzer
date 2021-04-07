# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import random
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from os.path import isfile, splitext
from routine_qiime2_analyses.analyses_prep import AnalysisPrep
from routine_qiime2_analyses._routine_q2_cmds import (
    get_case, get_new_meta_pd, run_import, run_export, filter_feature_table,
)
from routine_qiime2_analyses._routine_q2_metadata import (
    get_train_perc_from_numeric, get_meta_subset, get_cat_vars_and_vc,
    make_train_test_from_cat
)
from routine_qiime2_analyses._routine_q2_io_utils import (
    read_meta_pd, get_analysis_folder, filter_mb_table, filter_non_mb_table,
    write_filtered_tsv

)
from routine_qiime2_analyses._routine_q2_mmvec import (
    get_mmvec_dicts, rename_duplicate_columns,
)


class PairedData(object):

    analyses_commands = {}

    def __init__(self, config, project) -> None:
        self.cmds = {}
        paired_config = get_mmvec_dicts(config.p_mmvec_pairs)
        self.pairs = paired_config[0]
        self.filtering = paired_config[1]
        self.params = paired_config[2]
        self.subsets = paired_config[3]
        self.mmvecs = pd.DataFrame()

        self.get_mmvec_matrix()
        self.datasets_paths = self.make_datasets_paths(config, project)
        self.unstack_mmvecs()
        self.get_common_paths(config)

    def unstack_mmvecs(self):
        mmvecs_us = self.mmvecs.set_index(
            ['pair', 'filter', 'subset', 'omic']
        ).unstack().reset_index()
        mmvecs_us = mmvecs_us.drop(columns=[('variable', 1), ('factors', 1)])
        mmvecs_us.columns = ['%s%s' % (x, y) if x not in ['variable', 'factors']
                             else x for x, y in mmvecs_us.columns]
        self.mmvecs = mmvecs_us

    def get_dataset_path(self, dat, filter, subset):
        qza, meta = self.datasets_paths.loc[
            (self.datasets_paths['dataset'] == dat) &
            (self.datasets_paths['filter'] == filter) &
            (self.datasets_paths['subset'] == subset),
            ['qza', 'meta']].values[0]
        return qza, meta

    def get_traintests(self, meta_fp, new_meta_pd, vars, train, train_col):
        if train.isdigit() or train.replace('.', '').isdigit():
            train_perc = get_train_perc_from_numeric(train, new_meta_pd)
            vars_pd = new_meta_pd[vars].copy()
            cat_vars, cat_pd, vc, rep_d = get_cat_vars_and_vc(vars, vars_pd)
            if cat_vars and vc.size < cat_pd.shape[0] * 0.5:
                train_samples = make_train_test_from_cat(
                    cat_pd, vc, train_perc, meta_fp, cat_vars, train_col, rep_d)
            else:
                train_samples = random.sample(
                    new_meta_pd.index.tolist(),
                    k=int(train_perc * new_meta_pd.shape[0]))
            return train_samples
        return None

    def make_train_test_column(self, meta_fp, train_test_d,
                               meta_pd, dat1, dat2) -> dict:
        train_tests = {}
        train = train_test_d['train']
        meta_tt_pd = meta_pd.set_index('sample_name').copy()
        for dat in [dat1, dat2]:
            if 'datasets' in train_test_d and dat in train_test_d['datasets']:
                for tt, vars in train_test_d['datasets'][dat].items():
                    vars_pd = meta_tt_pd[vars].copy()
                    vars_pd = vars_pd.loc[~vars_pd.isna().any(1)]
                    vars_pd = rename_duplicate_columns(vars_pd)
                    trains = self.get_traintests(meta_fp, vars_pd, vars,
                                                 str(train), tt)
                    if trains:
                        train_tests[tt] = trains
                break
        return train_tests

    def get_new_fps(self, meta_dir, data_dir, qza1, qza2, dat1, prev1, abun1,
                    dat2, prev2, abun2, pair, nsams, config):
        meta_fp = '%s/meta_%s_%s_%s__%s_%s_%s__%s_%ss.tsv' % (
            meta_dir, dat1, prev1, abun1, dat2, prev2, abun2, pair, nsams)
        new_tsv1 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
            data_dir, dat1, prev1, abun1, pair, nsams)
        new_qza1 = '%s.qza' % splitext(new_tsv1)[0]
        new_tsv2 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
            data_dir, dat2, prev2, abun2, pair, nsams)
        new_qza2 = '%s.qza' % splitext(new_tsv2)[0]
        if config.force or not isfile(new_qza1):
            cmd = filter_feature_table(qza1, new_qza1, meta_fp)
            self.cmds.setdefault('import', []).append(cmd)
        if config.force or not isfile(new_tsv1):
            cmd = run_export(new_qza1, new_tsv1, 'FeatureTable')
            self.cmds.setdefault('import', []).append(cmd)
        if config.force or not isfile(new_qza2):
            cmd = filter_feature_table(qza2, new_qza2, meta_fp)
            self.cmds.setdefault('import', []).append(cmd)
        if config.force or not isfile(new_tsv2):
            cmd = run_export(new_qza2, new_tsv2, 'FeatureTable')
            self.cmds.setdefault('import', []).append(cmd)
        self.register_command('mmvec_paired_imports')
        return meta_fp, new_tsv1, new_qza1, new_tsv2, new_qza2

    def get_common_paths(self, config):
        paths = []
        pfs = ['pair', 'filter', 'subset']
        for (pair, filter, subset), mmvec in self.mmvecs.groupby(pfs):
            data_dir = get_analysis_folder(config.i_datasets_folder,
                'mmvec/common/data/%s/%s' % (pair, subset))
            meta_dir = get_analysis_folder(config.i_datasets_folder,
                'mmvec/common/metadata/%s/%s' % (pair, subset))
            mmvec_d = mmvec.iloc[0, :].to_dict()
            dat1, dat2 = mmvec_d['dataset1'], mmvec_d['dataset2']
            prev1, prev2 = mmvec_d['prevalence1'], mmvec_d['prevalence2']
            abun1, abun2 = mmvec_d['abundance1'], mmvec_d['abundance2']
            qza1, meta1 = self.get_dataset_path(dat1, filter, subset)
            qza2, meta2 = self.get_dataset_path(dat2, filter, subset)
            if not isfile(meta1) or not isfile(meta2):
                continue
            meta1_pd, meta2_pd = read_meta_pd(meta1), read_meta_pd(meta2)
            sams = set(meta1_pd.sample_name) & set(meta2_pd.sample_name)
            if len(sams) < 10:
                print('Not enough samples in pair %s: %s (%s) vs %s (%s)' % (
                    pair, mmvec_d['dataset1'], meta1_pd.shape[0],
                          mmvec_d['dataset2'], meta2_pd.shape[0]))
                continue
            meta_fp, new_tsv1, new_qza1, new_tsv2, new_qza2 = self.get_new_fps(
                meta_dir, data_dir, qza1, qza2, dat1, prev1, abun1,
                dat2, prev2, abun2, pair, len(sams), config)
            meta_subset = get_meta_subset(meta1_pd, meta2_pd, sams)
            meta_subset.to_csv(meta_fp, index=False, sep='\t')
            paths.append([pair, filter, subset, sams, meta_fp,
                          new_tsv1, new_tsv2, new_qza1, new_qza2])
            print('\t\t\t* [TODO]', pair, filter, subset, ':',
                  dat1, 'vs', dat2, '(%s samples)' % meta_subset.shape[0])
        if paths:
            common_paths_pd = pd.DataFrame(paths, columns=(pfs + ['common_sams',
                'meta_fp', 'new_tsv1', 'new_tsv2', 'new_qza1', 'new_qza2']))
            self.mmvecs = self.mmvecs.merge(common_paths_pd,
                                            on=['pair', 'filter', 'subset'])

    def make_train_test(self, config):
        for _, mmvec in self.mmvecs.groupby(['pair', 'filter', 'subset']):
            d = mmvec.iloc[0, :].to_dict()
            fps = ['dataset1', 'dataset2', 'meta_fp', 'new_tsv1',
                   'new_tsv2', 'new_qza1', 'new_qza2']
            dat1, dat2, meta_fp, tsv1, tsv2, qza1, qza2 = [d[x] for x in fps]
            meta_subset = read_meta_pd(meta_fp)
            train_tests = self.make_train_test_column(
                meta_fp, config.train_test_dict, meta_subset, dat1, dat2)
            rewrite = False
            meta_subset_cols = set(meta_subset.columns)
            for train_col, train_samples in train_tests.items():
                if train_col not in meta_subset_cols:
                    rewrite = True
                    meta_subset[train_col] = [
                        'Train' if x in set(train_samples) else
                        'Test' for x in meta_subset.sample_name.tolist()]
            if config.force or rewrite:
                meta_subset.to_csv(meta_fp, index=False, sep='\t')

    def make_datasets_paths(self, config, project):
        datasets_path = self.get_datasets_paths(config)
        for (dataset, filter, subset), row in datasets_path.groupby(
                ['dataset', 'filter', 'subset']):
            row_d = row.iloc[0, :].to_dict()
            tsv, qza, meta = row_d['tsv'], row_d['qza'], row_d['meta']
            if isfile(tsv) and isfile(qza) and isfile(meta):
                continue
            data = project.datasets[dataset]
            variable, factors = row_d['variable'], row_d['factors']
            meta_pd = get_new_meta_pd(data.metadata[0], subset,
                                      variable, factors)
            tsv_pd = data.data[0][meta_pd.sample_name.tolist()]
            preval, abund = row_d['prevalence'], row_d['abundance']
            if row_d['is_mb']:
                tsv_pd, res = filter_mb_table(preval, abund, tsv_pd)
            else:
                tsv_pd, res = filter_non_mb_table(preval, abund, tsv_pd)
            meta_pd.to_csv(meta, index=False, sep='\t')
            if config.force or not isfile(tsv):
                write_filtered_tsv(tsv, tsv_pd)
            if config.force or not isfile(qza):
                cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                self.cmds.setdefault(dataset, []).append(cmd)
        self.register_command('mmvec_single_imports')
        return datasets_path

    def get_datasets_paths(self, config):
        datasets_paths = self.mmvecs.copy()
        datasets_paths = datasets_paths.drop(columns=['pair', 'omic'])
        datasets_paths = datasets_paths.loc[
            ~datasets_paths.astype(str).duplicated()]
        paths = []
        for r, row in datasets_paths.iterrows():
            dataset = row['dataset']
            filter = row['filter']
            subset = row['subset']
            odir = get_analysis_folder(
                config.i_datasets_folder, 'mmvec/datasets/%s/%s' % (
                    dataset, subset))
            rad = '%s_%s' % (dataset, filter)
            tsv = '%s/tab_%s.tsv' % (odir, rad)
            qza = '%s.qza' % splitext(tsv)[0]
            meta = '%s/meta_%s.tsv' % (odir, rad)
            paths.append([tsv, qza, meta])
        datasets_paths = pd.concat([
            datasets_paths, pd.DataFrame(paths, columns=['tsv', 'qza', 'meta'])
        ], axis=1)
        return datasets_paths

    def get_mmvec_matrix(self):
        self.make_mmvec_pairs_table()
        self.merge_filtering_per_pair()
        self.merge_subsets_apply()

    def make_mmvec_pairs_table(self):
        # make the mmvec pairs table
        mmvec_dat_cols = ['pair', 'omic', 'dataset', 'is_mb']
        self.mmvecs = pd.DataFrame([([pair, (omic+1)] + list(dataset))
                                    for pair, paired in self.pairs.items()
                                    for omic, dataset in enumerate(paired)],
                                   columns=mmvec_dat_cols)

    def merge_filtering_per_pair(self):
        # merge the filtering per pair
        filts_df = []
        for pair, filts_dats in self.filtering.items():
            for filt, prevs_abunds in filts_dats.items():
                for dat_mb, prev_abun in prevs_abunds.items():
                    filts_df.append(([pair, filt] + list(dat_mb) + prev_abun))
        if filts_df:
            filts = pd.DataFrame(filts_df, columns=[
                'pair', 'filter', 'dataset', 'is_mb',
                'prevalence', 'abundance'])
            self.mmvecs = self.mmvecs.merge(
                filts, on=['pair', 'dataset', 'is_mb'])
        else:
            self.mmvecs['filter'] = '0_0'
            self.mmvecs['prevalence'] = '0'
            self.mmvecs['abundance'] = '0'

    def merge_subsets_apply(self):
        # merge the subsets to apply
        subsets_fp = [
            [pair, var, subset, get_case(subset, var)]
            for var, subsets in self.subsets.items()
            for subset in subsets
            for pair in self.mmvecs.pair.unique()]
        if subsets_fp:
            subsets = pd.DataFrame(
                subsets_fp, columns=['pair', 'variable', 'factors', 'subset'])
            self.mmvecs = self.mmvecs.merge(subsets, on=['pair'], how='outer')

    def get_paired_unique_datasets(self, project):
        unique_datasets = set()
        for dats in self.pairs.values():
            if len(dats) != 2:
                print('Only two datasets must be provided by pair')
                continue
            for dat in dats:
                if dat not in project.datasets:
                    continue
                unique_datasets.add(dat)
        return list(unique_datasets)

    def get_params_combinations(self):
        params = []
        params_list = ['train_column', 'n_examples', 'batches', 'learns',
                       'epochs', 'priors', 'thresh_feats', 'latent_dims']
        for param in params_list:
            print(param, self.params[param])
        return params

    def mmvec(self, config, project):
        params = self.get_params_combinations()
        for r, row in self.mmvecs.iterrows():
            break

    def register_command(self, analysis):
        AnalysisPrep.analyses_commands[analysis] = self.cmds

# class CreateScripts(object):
#
#     def __init__(self):
#         pass
#
#     def get_job_folders(self, analysis):
#         job_folder = get_job_folder(
#             self.config.i_datasets_folder, analysis)
#         job_folder2 = get_job_folder(
#             self.config.i_datasets_folder, '%s/chunks' % analysis)
#         return job_folder, job_folder2
#
#     def get_analysis_folder(self, dat):
#         return get_analysis_folder(
#             self.config.i_datasets_folder, '%s/%s' % (self.analysis, dat))
#
#     def print_message(self):
#         if self.main_written:
#             print_message(
#                 '# %s' % self % self.analysis,
#                 'sh', self.run_pbs, self.jobs)
