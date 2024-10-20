# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import sys
import random
import itertools
import pandas as pd
from os.path import isdir, isfile, splitext
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import (
    run_import, run_export, write_sample_filter, write_mmvec)
from microbiome_analyzer.core.metadata import (
    get_subset_pd, get_train_perc_from_numeric, get_meta_subset,
    get_cat_vars_and_vc, make_train_test_from_cat, rename_duplicate_columns)
from microbiome_analyzer.analyses.filter import (
    filter_mb_table, filter_non_mb_table)
from microbiome_analyzer._inputs import read_meta_pd
from microbiome_analyzer._scratch import io_update, to_do, rep
from microbiome_analyzer._io_utils import write_filtered_tsv


class PairedData(object):

    def __init__(self, config, project) -> None:
        if config.mmvec_pairs_fp:
            self.config = config
            self.project = project
            self.dir = project.dir
            self.dirs = project.dirs
            self.out = ''
            (self.pairs, self.filtering, self.params,
             self.subsets) = self.get_mmvec_dicts()
            self.mmvecs = pd.DataFrame()
            self.ios = {}
            self.cmds = {}
            self.analysis = ''
            self.params_list = [
                'train_column', 'n_examples', 'batches', 'learns', 'epochs',
                'input_prior', 'output_prior', 'thresh_feats', 'latent_dims',
                'summary_interval']
            self.get_mmvec_matrix(project)
            if self.mmvecs.shape[0]:
                self.datasets_paths = self.make_datasets_paths(project)
                self.unstack_mmvecs()
                self.get_common_paths()
        self.mmvec_pd = pd.DataFrame()

    def get_pairs(self) -> dict:
        """
        Get the parameters for mmvec pairs to process.
        """
        if 'pairs' not in self.config.mmvec_pairs:
            print('[mmvec] No datasets pairs specified:\nExiting\n')
            sys.exit(0)
        mmvec_pairs = {}
        for pair, paired_datasets in self.config.mmvec_pairs['pairs'].items():
            n_dats = len(paired_datasets)
            if n_dats != 2:
                print('[mmvec] Must be two datasets per mmvec pair (found %s)\n'
                      'Exiting\n' % n_dats)
                sys.exit(0)
            paired = []
            for dat in paired_datasets:
                if dat[-1] == '*':
                    paired.append((dat[:-1], 1))
                else:
                    paired.append((dat, 0))
            mmvec_pairs[pair] = paired
        return mmvec_pairs

    @staticmethod
    def get_dat_mb_or_not(dat: str) -> tuple:
        if dat[-1] == '*':
            return dat[:-1], 1
        else:
            return dat, 0

    def get_filtering(self, mmvec_pairs: dict) -> dict:
        """
        Get the parameters for songbird passed by the user.
        """
        dats = []
        filtering = {}
        for pair, dats_pair in mmvec_pairs.items():
            if pair not in filtering:
                filtering[pair] = {'0_0': {}}
            for dat in dats_pair:
                dats.append(dat)
                filtering[pair]['0_0'][dat] = ['0', '0']
        if 'filtering' not in self.config.mmvec_pairs:
            print('No filtering thresholds set in '
                  '%s\n:' % self.config.mmvec_pairs_fp)
        if self.config.mmvec_pairs.get('filtering'):
            if 'global' in self.config.mmvec_pairs['filtering']:
                for filt_name, prev_abund in self.config.mmvec_pairs[
                        'filtering']['global'].items():
                    for pair, dats_pair in mmvec_pairs.items():
                        if pair not in filtering:
                            filtering[pair] = {}
                        if filt_name not in filtering[pair]:
                            filtering[pair][filt_name] = {}
                        for dat in dats_pair:
                            dats.append(dat)
                            filtering[pair][filt_name][dat] = prev_abund
            for pair, pair_d in self.config.mmvec_pairs['filtering'].items():
                if pair == 'global':
                    continue
                filtering[pair] = {}
                for filt_name, dats_d in pair_d.items():
                    filtering[pair][filt_name] = {}
                    for dat_, prev_abund in dats_d.items():
                        dat = self.get_dat_mb_or_not(dat_)
                        if dat in dats:
                            filtering[pair][filt_name][dat] = prev_abund
        return filtering

    def get_mmvec_params(self) -> dict:
        """
        Get the parameters for songbird passed by the user.
        """
        params = {
            'train_column': ['None'],
            'n_examples': ['10'],
            'batches': ['2'],
            'learns': ['1e-4'],
            'epochs': ['5000'],
            'input_prior': ['0.8'],
            'output_prior': ['0.8'],
            'thresh_feats': ['0'],
            'latent_dims': ['3'],
            'summary_interval': ['10']
        }
        if 'params' not in self.config.mmvec_pairs:
            print('No parameters set in %s:\nUsing defaults: %s' % (
                self.config.mmvec_pairs_fp, ', '.join([
                    '%s: %s' % (k, v) for k, v in params.items()])))
        else:
            for param, cur_param in self.config.mmvec_pairs['params'].items():
                if not isinstance(cur_param, list):
                    print('Parameter %s should be a list (correct in %s)\n' % (
                        param, self.config.mmvec_pairs_fp))
                    sys.exit(0)
                params[param] = cur_param
        return params

    def get_subsets(self) -> dict:
        subsets = {'ALL': {}}
        if 'subsets' in self.config.mmvec_pairs:
            subsets.update(self.config.mmvec_pairs['subsets'])
        return subsets

    def get_mmvec_dicts(self):
        """
        Collect pairs of datasets from the passed yaml file:
        """
        mmvec_pairs = self.get_pairs()
        mmvec_filtering = self.get_filtering(mmvec_pairs)
        mmvec_params = self.get_mmvec_params()
        mmvec_subsets = self.get_subsets()
        return mmvec_pairs, mmvec_filtering, mmvec_params, mmvec_subsets

    def unstack_mmvecs(self):
        # print(self.mmvecs)
        mmvecs_us = self.mmvecs.set_index(
            ['pair', 'filter', 'subset', 'omic']
        ).unstack().reset_index()
        # print(1)
        # print(mmvecs_us)
        mmvecs_us.columns = ['%s%s' % (x, y) for x, y in mmvecs_us.columns]
        # print(2)
        # print(mmvecs_us)
        self.mmvecs = mmvecs_us

    def get_dataset_path(self, dat, filter_, subset):
        qza, meta = self.datasets_paths.loc[
            (self.datasets_paths['dataset'] == dat) &
            (self.datasets_paths['filter'] == filter_) &
            (self.datasets_paths['subset'] == subset),
            ['qza', 'meta']].values[0]
        return qza, meta

    @staticmethod
    def get_traintests(meta_fp, new_meta_pd, vars_, train, train_col):
        if train.isdigit() or train.replace('.', '').isdigit():
            train_perc = get_train_perc_from_numeric(train, new_meta_pd)
            vars_pd = new_meta_pd[vars_].copy()
            cat_vars, cat_pd, vc, rep_d = get_cat_vars_and_vc(vars_, vars_pd)
            if cat_vars and vc.size < cat_pd.shape[0] * 0.5:
                train_samples = make_train_test_from_cat(
                    new_meta_pd, cat_pd, vc, train_perc, meta_fp, cat_vars,
                    train_col, rep_d)
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
                for tt, vars_ in train_test_d['datasets'][dat].items():
                    vars_pd = meta_tt_pd[vars_].copy()
                    vars_pd = vars_pd.loc[~vars_pd.isna().any(axis=1)]
                    vars_pd = rename_duplicate_columns(vars_pd)
                    trains = self.get_traintests(meta_fp, vars_pd, vars_,
                                                 str(train), tt)
                    if trains:
                        train_tests[tt] = trains
                break
        return train_tests

    def get_new_fps(self, meta_dir, data_dir, qza1, qza2, dat1, prev1, abun1,
                    dat2, prev2, abun2, pair, nsams):
        meta_fp = '%s/meta_%s_%s_%s__%s_%s_%s__%s_%ss.tsv' % (
            meta_dir, dat1, prev1, abun1, dat2, prev2, abun2, pair, nsams)
        new_tsv1 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
            data_dir, dat1, prev1, abun1, pair, nsams)
        new_qza1 = '%s.qza' % splitext(new_tsv1)[0]
        new_tsv2 = '%s/tab_%s_%s_%s__%s_%ss.tsv' % (
            data_dir, dat2, prev2, abun2, pair, nsams)
        new_qza2 = '%s.qza' % splitext(new_tsv2)[0]
        cmd = ''
        if self.config.force or to_do(new_qza1):
            cmd += write_sample_filter(qza1, new_qza1, meta_fp)
            io_update(self, i_f=[qza1, meta_fp], o_f=new_qza1, key='import')
        if self.config.force or to_do(new_tsv1):
            cmd += run_export(new_qza1, new_tsv1, 'FeatureTable')
            io_update(self, o_f=new_tsv1, key='import')
            if isfile(new_qza1):
                io_update(self, i_f=new_qza1, key='import')
        if self.config.force or to_do(new_qza2):
            cmd += write_sample_filter(qza2, new_qza2, meta_fp)
            io_update(self, i_f=[qza2, meta_fp], o_f=new_qza2, key='import')
        if self.config.force or to_do(new_tsv2):
            cmd += run_export(new_qza2, new_tsv2, 'FeatureTable')
            io_update(self, o_f=new_tsv2, key='import')
            if not to_do(new_qza2):
                io_update(self, i_f=new_qza2, key='import')
        if cmd:
            self.cmds.setdefault('import', []).append(cmd)
        return meta_fp, new_tsv1, new_qza1, new_tsv2, new_qza2

    def get_common_paths(self):
        paths = []
        self.analysis = 'mmvec_paired_imports'
        # print(self.mmvecs)
        pfs = ['pair', 'filter', 'subset']
        for (pair, filter_, subset), mmvec in self.mmvecs.groupby(pfs):
            data_dir = self.get_output('common/data/%s/%s' % (pair, subset))
            meta_dir = self.get_output('common/metadata/%s/%s' % (pair, subset))
            mmvec_d = mmvec.iloc[0, :].to_dict()
            dat1_, dat2_ = mmvec_d['dataset1'], mmvec_d['dataset2']
            dat1 = dat1_.replace('/', '__')
            dat2 = dat2_.replace('/', '__')
            prev1, prev2 = mmvec_d['prevalence1'], mmvec_d['prevalence2']
            abun1, abun2 = mmvec_d['abundance1'], mmvec_d['abundance2']
            qza1, meta1 = self.get_dataset_path(dat1_, filter_, subset)
            qza2, meta2 = self.get_dataset_path(dat2_, filter_, subset)
            if to_do(meta1) or to_do(meta2):
                continue
            meta1_pd = read_meta_pd(rep(meta1))
            meta2_pd = read_meta_pd(rep(meta2))
            sams = set(meta1_pd.sample_name) & set(meta2_pd.sample_name)
            if len(sams) < 9:
                print('- Too few samples (%s) for "%s": %s (%s) vs %s (%s)' % (
                    pair, len(sams),
                    mmvec_d['dataset1'], meta1_pd.shape[0],
                    mmvec_d['dataset2'], meta2_pd.shape[0]))
                continue
            meta_fp, new_tsv1, new_qza1, new_tsv2, new_qza2 = self.get_new_fps(
                meta_dir, data_dir, qza1, qza2, dat1, prev1, abun1,
                dat2, prev2, abun2, pair, len(sams))
            meta_subset = get_meta_subset(meta1_pd, meta2_pd, sams)
            meta_subset.to_csv(rep(meta_fp), index=False, sep='\t')
            paths.append([pair, filter_, subset, sams, meta_fp,
                          new_tsv1, new_tsv2, new_qza1, new_qza2])
        if paths:
            common_paths_pd = pd.DataFrame(paths, columns=(
                    pfs + ['common_sams', 'meta_fp', 'new_tsv1', 'new_tsv2',
                           'new_qza1', 'new_qza2']))
            self.mmvecs = self.mmvecs.merge(
                common_paths_pd, on=['pair', 'filter', 'subset'])

        self.register_io_command()

    def make_train_test(self):
        if self.mmvecs.shape[0]:
            for _, mmvec in self.mmvecs.groupby(['pair', 'filter', 'subset']):
                d = mmvec.iloc[0, :].to_dict()
                fps = ['dataset1', 'dataset2', 'meta_fp', 'new_tsv1',
                       'new_tsv2', 'new_qza1', 'new_qza2']
                dat1, dat2, meta_fp, tsv1, tsv2, qza1, qza2 = [
                    d[x] for x in fps]
                meta_subset = read_meta_pd(rep(meta_fp))
                train_tests = self.make_train_test_column(
                    rep(meta_fp), self.config.train_test_dict,
                    meta_subset, dat1, dat2)
                rewrite = False
                meta_subset_cols = set(meta_subset.columns)
                for train_col, train_samples in train_tests.items():
                    if train_col not in meta_subset_cols:
                        rewrite = True
                        meta_subset[train_col] = [
                            'Train' if x in set(train_samples) else
                            'Test' for x in meta_subset.sample_name.tolist()]
                if self.config.force or rewrite:
                    meta_subset.to_csv(rep(meta_fp), index=False, sep='\t')

    def make_datasets_paths(self, project):
        datasets_path = self.get_datasets_paths()
        self.analysis = 'mmvec_single_imports'
        for (dat, _, subset), row in datasets_path.groupby(
                ['dataset', 'filter', 'subset']):
            row_d = row.iloc[0, :].to_dict()
            tsv, qza, meta = row_d['tsv'], row_d['qza'], row_d['meta']
            if not to_do(tsv) and not to_do(qza) and not to_do(meta):
                continue
            data = project.datasets[dat]
            subsets = row_d['subsets']
            meta_pd = get_subset_pd(data.metadata, dat, subset, subsets)
            if meta_pd is None:
                continue
            tsv_pd = data.data[''].to_dataframe(dense=True)
            sams = list(set(meta_pd.sample_name) & set(tsv_pd.columns))
            tsv_pd = tsv_pd[sams]
            preval, abund = row_d['prevalence'], row_d['abundance']
            if row_d['is_mb']:
                tsv_pd, _ = filter_mb_table(preval, abund, tsv_pd)
            else:
                tsv_pd, _ = filter_non_mb_table(preval, abund, tsv_pd)
            meta_pd.to_csv(rep(meta), index=False, sep='\t')
            if self.config.force or to_do(tsv):
                write_filtered_tsv(rep(tsv), tsv_pd)
                io_update(self, o_f=tsv, key=dat)
            if self.config.force or to_do(qza):
                cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                io_update(self, o_f=qza, key=dat)
                if not to_do(tsv):
                    io_update(self, i_f=tsv, key=dat)
                self.cmds.setdefault(dat, []).append(cmd)
        self.register_io_command()
        return datasets_path

    def get_datasets_paths(self):
        self.analysis = 'mmvec'
        datasets_paths = self.mmvecs.copy()
        datasets_paths = datasets_paths.drop(columns=['pair', 'omic'])
        datasets_paths = datasets_paths.loc[
            ~datasets_paths.astype(str).duplicated()]
        ps = []
        for r, row in datasets_paths.iterrows():
            dat = row['dataset']
            filter_ = row['filter']
            subset = row['subset']
            self.get_output('datasets/%s/%s' % (dat.replace('/', '__'), subset))
            rad = '%s_%s' % (dat.replace('/', '__'), filter_)
            tsv = '%s/tab_%s.tsv' % (self.out, rad)
            qza = '%s.qza' % splitext(tsv)[0]
            meta = '%s/meta_%s.tsv' % (self.out, rad)
            ps.append([dat, filter_, subset, tsv, qza, meta])

        paths_pd = pd.DataFrame(
            ps, columns=['dataset', 'filter', 'subset', 'tsv', 'qza', 'meta'])
        datasets_paths = datasets_paths.merge(
            paths_pd, on=['dataset', 'filter', 'subset'], how='left')
        return datasets_paths

    def get_mmvec_matrix(self, project):
        self.make_mmvec_pairs_table(project)
        # print(self.mmvecs)
        if self.mmvecs.shape[0]:
            self.merge_filtering_per_pair()
            self.merge_subsets_apply()

    def make_mmvec_pairs_table(self, project):
        mmvec_dat_cols = ['pair', 'omic', 'dataset', 'is_mb']
        mmvecs = []
        for pair, paired in self.pairs.items():
            cur_pair = []
            for omic, dataset in enumerate(paired):
                if dataset[0] not in project.datasets:
                    # cur_pair.append(([pair, (omic + 1)] + list(dataset)))
                    break
                cur_pair.append(([pair, (omic + 1)] + list(dataset)))
            else:
                mmvecs.extend(cur_pair)
        if mmvecs:
            self.mmvecs = pd.DataFrame(mmvecs, columns=mmvec_dat_cols)

    def merge_filtering_per_pair(self):
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
        subsets_fp = [[pair, subset, subsets]
                      for subset, subsets in self.subsets.items()
                      for pair in self.mmvecs.pair.unique()]
        if subsets_fp:
            subsets = pd.DataFrame(
                subsets_fp, columns=['pair', 'subset', 'subsets'])
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
        """Make a pandas data frame from the combinations
        of mmvec run/hyper-parameters. It includes the
        handling of user-specified 'train_column', which
        always take precedence over the default 'n_examples'.

        Returns
        -------
        params_pd : pd.DataFrame
            Comobinations of parameters as rows, and
            individual parameters as columns.
        """
        params = []
        to_combine = [self.params[param] for param in self.params_list]
        for params_combination in itertools.product(*to_combine):
            params.append(params_combination)
        params_pd = pd.DataFrame(params, columns=self.params_list).astype(str)
        params_pd.loc[params_pd['train_column'] != 'None', 'n_examples'] = ''
        params_pd = params_pd.drop_duplicates()
        return params_pd

    @staticmethod
    def process_params_combinations(
            row: pd.Series,
            params_pd: pd.DataFrame,
            mess: set):
        """Filter the combinations of parameters too remove
        those involving unusable train/test splits, e.g. not
        having the specified or too few samples therein.

        Parameters
        ----------
        row : pd.Series
            Info associated with onoe pair of datasets.
        params_pd : pd.Series
            Combinations of parameters (rows)
        mess : set
            Messages to print
        """
        valid_params = []
        ncommon = len(row['common_sams'])
        meta = row['meta_fp']
        meta_cols = pd.read_table(rep(meta), nrows=1, index_col=0).columns
        for p, params in params_pd.iterrows():
            n_test, train_column = params['n_examples'], params['train_column']
            if train_column != 'None':
                if train_column not in set(meta_cols):
                    valid_params.append(p)
                    m = 'Training column "%s" not in metadata' % train_column
                    if m not in mess:
                        mess.add(m)
                        print(m)
                    continue
                meta_pd = pd.read_table(
                    rep(meta), usecols=['sample_name', train_column])
                ntrain = meta_pd[train_column].value_counts()['Train']
                if ncommon < (1.2 * ntrain):
                    valid_params.append(p)
                    print('\tskipping pair "%s" (%s samples for %s training'
                          ' samples]):' % (row['pair'], ncommon, ntrain))
                    continue
            else:
                if int(ncommon) < (1.2 * int(n_test)):
                    valid_params.append(p)
                    print('\tskipping pair "%s" (too few samples [%s samples '
                          'for %s examples]):' % (row['pair'], ncommon, n_test))
                    continue
        if valid_params:
            params_pd.drop(index=valid_params, inplace=True)

    def get_res_dir(self, params):
        """

        Parameters
        ----------
        params : pd.Seris
        config : Class

        Returns
        -------

        """
        res_dir = 'b%s_l%s_e%s_ip%s_op%s_f%s_d%s_t%s_n%s_s%s_gpu-%s' % (
            params['batches'],
            params['learns'],
            params['epochs'],
            params['input_prior'],
            params['output_prior'],
            params['thresh_feats'],
            params['latent_dims'],
            params['train_column'],
            params['n_examples'],
            params['summary_interval'],
            str(self.config.gpu)[0]
        )
        return res_dir

    def get_out(self, model_null: str) -> tuple:
        """
        Parameters
        ----------
        model_null : str
            "model" or null""

        Returns
        -------
        mod_nul_dir : str
        mod_nul_rnk : str
        mod_nul_ord : str
        mod_nul_stt : str
        """
        mod_nul_dir = '%s/%s' % (self.out, model_null)
        if not isdir(rep(mod_nul_dir)):
            os.makedirs(rep(mod_nul_dir))
        mod_nul_rnk = '%s/ranks.tsv' % mod_nul_dir
        mod_nul_ord = '%s/ordination.txt' % mod_nul_dir
        mod_nul_stt = '%s/stats.qza' % mod_nul_dir
        return mod_nul_dir, mod_nul_rnk, mod_nul_ord, mod_nul_stt

    def get_mmvec_pd(self, mmvec):
        self.mmvec_pd = pd.DataFrame(
            mmvec, columns=[
                'pair', 'filter', 'subset',  # 'case',
                'omic1', 'omic2', 'pr1', 'ab1', 'pr2', 'ab2', 'n_common',
                'meta_common_fp', 'omic1_common_fp', 'omic2_common_fp',
                'omic1_common_qza', 'omic2_common_qza', 'mmvec_parameters',
                'mmvec_out'
            ])

    def get_output(self, dat) -> str:
        out = '%s/%s/%s' % (self.dir, self.analysis, dat)
        if not isdir(rep(out)):
            os.makedirs(rep(out))
        self.out = out
        return out

    def mmvec(self) -> None:
        """Main script for the creation of mmvec jobs.
        It iterates over the rows of the table created
        upfront and over each combination of parameters
        and collect the output info for potential reuse
        in figure generation and post-analysis.

        Parameters
        ----------
        config : Class instance of AnalysesConfig
            Contains all the routine analyses config info.
        """
        mess = set()
        mmvec = []
        self.analysis = 'mmvec'
        params_pd = self.get_params_combinations()
        for r, row in self.mmvecs.iterrows():
            self.process_params_combinations(row, params_pd, mess)
            pair, filter_, subset = row['pair'], row['filter'], row['subset']
            d1, p1, a1 = row['dataset1'], row['prevalence1'], row['abundance1']
            d2, p2, a2 = row['dataset2'], row['prevalence2'], row['abundance2']
            for p, params in params_pd.iterrows():
                res_dir = self.get_res_dir(params)
                self.get_output('paired/%s/%s/%s_%s-%s__%s_%s-%s/%s' % (
                    pair, subset, d1, p1, a1, d2, p2, a2, res_dir))
                mod_dir, mod_rnk, mod_rdn, mod_stt = self.get_out('model')
                nul_dir, nul_rnk, nul_rdn, nul_stt = self.get_out('null')
                summary = '%s/paired-summary.qzv' % self.out
                mmvec.append([
                    pair, filter_, subset, d1, d2, p1, a1, p2, a2,
                    len(row['common_sams']), row['meta_fp'], row['new_tsv1'],
                    row['new_tsv2'], row['new_qza1'], row['new_qza2'],
                    'mmvec_out__%s' % res_dir, self.out])
                if self.config.force or to_do(summary):
                    cmd = write_mmvec(
                        self, pair, row['meta_fp'], row['new_qza1'],
                        row['new_qza2'], res_dir, mod_dir, nul_dir, mod_rnk,
                        mod_rdn, mod_stt, nul_rnk, nul_rdn, nul_stt, summary,
                        params['batches'], params['learns'], params['epochs'],
                        params['input_prior'], params['output_prior'],
                        params['thresh_feats'], params['latent_dims'],
                        params['train_column'], params['n_examples'],
                        params['summary_interval'], self.config.gpu,
                        self.config.qiime_env)
                    self.cmds.setdefault(pair, []).append(cmd)
        if mmvec:
            self.get_mmvec_pd(mmvec)
        self.register_io_command()

    def register_io_command(self):
        AnalysisPrep.analyses_ios[self.analysis] = dict(self.ios)
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.ios = {}
        self.cmds = {}
