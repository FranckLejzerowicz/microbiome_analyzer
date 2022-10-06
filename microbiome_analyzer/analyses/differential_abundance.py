# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import sys
import random
import itertools
import pandas as pd
from os.path import dirname, isdir, splitext
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import run_import, write_songbird
from microbiome_analyzer.core.metadata import (
    get_subset, get_subset_pd, rename_duplicate_columns,
    get_train_perc_from_numeric, get_cat_vars_and_vc, make_train_test_from_cat)
from microbiome_analyzer.analyses.filter import (
    filter_mb_table, filter_non_mb_table)
from microbiome_analyzer._io_utils import write_filtered_tsv
from microbiome_analyzer._inputs import read_meta_pd
from microbiome_analyzer._scratch import io_update, to_do, rep


class DiffModels(object):

    def __init__(self, config, project) -> None:
        self.ios = {}
        self.cmds, self.bcmds, self.fcmds = {}, {}, {}
        self.analysis = None
        self.config = config
        self.project = project
        self.dir = project.dir
        self.dirs = project.dirs
        self.out = ''
        if config.diff_models:
            (self.songbird_models, self.filtering, self.params, self.baselines,
             self.songbird_subsets) = self.get_songbird_dicts()
            self.models, self.models_issues = {}, {}
            self.songbirds = pd.DataFrame(dtype='object', columns=[
                'dataset', 'is_mb', 'filter', 'prevalence', 'abundance'])
            self.params_list = [
                'train', 'batches', 'learns', 'epochs', 'diff_priors',
                'thresh_feats', 'thresh_samples', 'summary_interval']
            self.q2s_pd = pd.DataFrame()
        self.songbird_pd = pd.DataFrame()

    def get_output(self, dat: str = '') -> str:
        out = '%s/%s' % (self.dir, self.analysis)
        if dat:
            out += '/%s' % dat
        if not isdir(rep(out)):
            os.makedirs(rep(out))
        self.out = out
        return out

    def get_songbird_models(self):
        if 'models' not in self.config.diff_models:
            print('No models in %s' % self.config.diff_models)
            sys.exit(0)
        is_mb = {}
        models = {}
        for dat in self.config.diff_models['models'].keys():
            if dat[-1] == '*':
                is_mb[dat[:-1]] = (dat[:-1], 1)
                models[dat[:-1]] = self.config.diff_models['models'][dat]
            else:
                is_mb[dat] = (dat, 0)
                models[dat] = self.config.diff_models['models'][dat]
        return models, is_mb

    def get_songbird_baselines(self):
        baselines = {}
        if 'baselines' in self.config.diff_models:
            return self.config.diff_models['baselines']
        return baselines

    def get_songbird_params(self):
        params = {
            'train': ['0.7'],
            'batches': ['2'],
            'learns': ['1e-4'],
            'epochs': ['5000'],
            'thresh_feats': ['0'],
            'thresh_samples': ['0'],
            'diff_priors': ['0.5'],
            'summary_interval': ['1']
        }
        if 'params' not in self.config.diff_models:
            print('No parameters set in %s:\nUsing defaults: %s' % (
                self.config.diff_models,
                ', '.join(['%s: %s' % (k, v) for k, v in params.items()])))
        else:
            for param in self.config.diff_models['params']:
                cur_param = self.config.diff_models['params'][param]
                if not isinstance(cur_param, list):
                    print('Parameter %s should be a list (correct in %s)\n' % (
                        param, self.config.diff_models))
                    sys.exit(0)
                params[param] = cur_param
        return params

    def get_filtering(self, models, is_mb):
        filtering = {
            '0_0': dict((is_mb[dat], ['0', '0']) for dat in models.keys())}
        if 'filtering' not in self.config.diff_models:
            pass
            # print('No filtering config in %s\n:' % self.config.diff_models)
        else:
            if 'global' in self.config.diff_models['filtering']:
                for fname, p_a in self.config.diff_models[
                        'filtering']['global'].items():
                    for dat in models.keys():
                        if fname not in filtering:
                            filtering[fname] = {}
                        filtering[fname][is_mb[dat]] = p_a
            for dat, filts in self.config.diff_models['filtering'].items():
                if dat == 'global':
                    continue
                for fname, p_a in filts.items():
                    if fname not in filtering:
                        filtering[fname] = {}
                    if dat in models.keys():
                        filtering[fname][is_mb[dat]] = p_a
        unique_filtering = {}
        for filter_name, dat_d in filtering.items():
            for dat, (preval, abund) in dat_d.items():
                unique_filtering.setdefault(dat, set()).add(
                    (filter_name, preval, abund))
        return unique_filtering

    def get_songbird_dicts(self):
        models, is_mb = self.get_songbird_models()
        filtering = self.get_filtering(models, is_mb)
        params = self.get_songbird_params()
        baselines = self.get_songbird_baselines()
        subsets = {'ALL': [[]]}
        if 'subsets' in self.config.diff_models:
            subsets.update(self.config.diff_models['subsets'])
        return models, filtering, params, baselines, subsets

    def merge_subsets_apply(self):
        subsets_fp = [
            [dataset, var, subset, get_subset(subset, var), '']
            for var, subsets in self.songbird_subsets.items()
            for subset in subsets
            for dataset in self.songbirds.dataset.unique()]
        if subsets_fp:
            subsets = pd.DataFrame(
                subsets_fp, columns=['dataset', 'variable', 'factors',
                                     'subset', 'pair'])
            self.songbirds = self.songbirds.merge(
                subsets, on=['dataset'], how='outer')

    def get_songbirds_filts(self):
        filts_df = []
        for (dat, is_mb), filts_dats in self.filtering.items():
            if dat not in self.project.datasets:
                continue
            for (filt, prev, abund) in filts_dats:
                filts_df.append([dat, is_mb, filt, prev, abund])
        if filts_df:
            self.songbirds = pd.DataFrame(filts_df, columns=[
                'dataset', 'is_mb', 'filter', 'prevalence', 'abundance'])

    def prep_songbirds(self, mmvec_pd):
        self.get_songbirds_filts()
        self.merge_subsets_apply()
        self.make_datasets_paths()
        self.merge_mmvecs(mmvec_pd)

    def merge_mmvecs(self, mmvec_pd):
        mmvecs = []
        for row in mmvec_pd.values:
            pair, filt, subset, dat1, dat2, prev1, abun1, prev2, abun2 = row[:9]
            meta_common_fp = row[10]
            omic1_common_qza = row[13]
            omic2_common_qza = row[14]
            mmvecs.append([dat1, filt, prev1, abun1, subset, pair,
                           omic1_common_qza, omic2_common_qza, meta_common_fp])
            mmvecs.append([dat2, filt, prev2, abun2, subset, pair,
                           omic1_common_qza, omic2_common_qza, meta_common_fp])
        if mmvecs and self.songbirds.shape[0]:
            self.songbirds.drop(
                columns=['is_mb', 'variable', 'factors'],
                inplace=True)
            self.songbirds = pd.concat([
                self.songbirds,
                pd.DataFrame(mmvecs, columns=self.songbirds.columns)])

    def make_datasets_paths(self):
        self.get_datasets_paths()
        self.analysis = 'songbird_imports'
        if self.songbirds.shape[0]:
            for (dat, _, subset), row in self.songbirds.groupby(
                    ['dataset', 'filter', 'subset']):
                row_d = row.iloc[0, :].to_dict()
                tsv, qza, meta = row_d['tsv'], row_d['qza'], row_d['meta']
                data = self.project.datasets[dat]
                variable, factors = row_d['variable'], row_d['factors']
                meta_pd = get_subset_pd(
                    data.metadata, subset, variable, factors)
                meta_pd.to_csv(rep(meta), index=False, sep='\t')
                if not self.config.force and not to_do(tsv) and not to_do(qza):
                    continue
                # tsv_pd = data.data[0].to_dataframe(dense=True)[
                tsv_pd = data.data[''].to_dataframe(dense=True)[
                    meta_pd.sample_name.tolist()]
                preval, abund = row_d['prevalence'], row_d['abundance']
                if row_d['is_mb']:
                    tsv_pd = filter_mb_table(preval, abund, tsv_pd)[0]
                else:
                    tsv_pd = filter_non_mb_table(preval, abund, tsv_pd)[0]
                if self.config.force or to_do(tsv):
                    write_filtered_tsv(rep(tsv), tsv_pd)
                    io_update(self, o_f=tsv, key=(dat, ''))
                if self.config.force or to_do(qza):
                    io_update(self, o_f=qza, key=(dat, ''))
                    if not to_do(tsv):
                        io_update(self, i_f=tsv, key=(dat, ''))
                    cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                    self.cmds.setdefault(dat, []).append(cmd)
        self.register_io_command()

    def get_datasets_paths(self):
        paths = []
        if self.songbirds.shape[0]:
            for r, row in self.songbirds.iterrows():
                dataset = row['dataset']
                filter_ = row['filter']
                subset = row['subset']
                for analysis in ['mmvec', 'songbird']:
                    self.analysis = analysis
                    self.get_output('datasets/%s/%s' % (dataset, subset))
                    rad = '%s_%s' % (dataset, filter_)
                    tsv = '%s/tab_%s.tsv' % (self.out, rad)
                    qza = '%s.qza' % splitext(tsv)[0]
                    meta = '%s/meta_%s.tsv' % (self.out, rad)
                    paths.append([tsv, qza, meta])
        if paths:
            self.songbirds = pd.concat([
                self.songbirds,
                pd.DataFrame(paths, columns=['tsv', 'qza', 'meta'])], axis=1)

    @staticmethod
    def get_traintests(meta_fp, new_meta_pd, vars_, train, train_col):
        if train.isdigit() or train.replace('.', '').isdigit():
            train_perc = get_train_perc_from_numeric(train, new_meta_pd)
            vars_pd = new_meta_pd[vars_].copy()
            cat_vars, cat_pd, vc, rep_d = get_cat_vars_and_vc(vars_, vars_pd)
            if cat_vars and vc.size < cat_pd.shape[0] * 0.5:
                train_samples = make_train_test_from_cat(
                    cat_pd, vc, train_perc, meta_fp, cat_vars, train_col, rep_d)
            else:
                train_samples = random.sample(
                    new_meta_pd.index.tolist(),
                    k=int(train_perc * new_meta_pd.shape[0]))
            return train_samples
        return None

    def make_train_test_column(self, meta_fp, meta_pd, dat) -> dict:
        train_test_d = self.config.train_test_dict
        train_tests = {}
        train = train_test_d['train']
        meta_tt_pd = meta_pd.set_index('sample_name').copy()
        if 'datasets' in train_test_d and dat in train_test_d['datasets']:
            for tt, vars_ in train_test_d['datasets'][dat].items():
                vars_pd = meta_tt_pd[vars_].copy()
                vars_pd = vars_pd.loc[~vars_pd.isna().any(1)]
                vars_pd = rename_duplicate_columns(vars_pd)
                trains = self.get_traintests(
                    meta_fp, vars_pd, vars_, str(train), tt)
                if trains:
                    train_tests[tt] = trains
        return train_tests

    def make_train_test(self):
        if self.songbirds.shape[0]:
            for _, sb in self.songbirds.groupby(
                    ['dataset', 'filter', 'subset']):
                d = sb.iloc[0, :].to_dict()
                fps = ['dataset', 'tsv', 'qza', 'meta']
                dat, tsv, qza, meta_fp = [d[x] for x in fps]
                meta_subset = read_meta_pd(rep(meta_fp))
                train_tests = self.make_train_test_column(
                    rep(meta_fp), meta_subset, dat)
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

    def get_params_combinations(self):
        """Make a pandas data frame from the combinations
        of songbird run/hyper-parameters. It includes the
        handling of user-specified 'train_column', which
        always take precedence over the default 'n_examples'.

        Returns
        -------
        params_pd : pd.DataFrame
            Combinations of params as rows, and individual params as columns.
        """
        params = []
        to_combine = [self.params[param] for param in self.params_list]
        for params_combination in itertools.product(*to_combine):
            params.append(params_combination)
        params_pd = pd.DataFrame(params, columns=self.params_list).astype(str)
        return params_pd

    @staticmethod
    def print_message_or_not(mess, m):
        if m not in mess:
            mess.add(m)

    def process_params_combinations(
            self,
            dataset: str,
            meta_pd: pd.DataFrame,
            params_pd: pd.DataFrame,
            mess: set):
        """Filter the combinations of parameters too remove
        those involving unusable train/test splits, e.g. not
        having the specified or too few samples therein.

        Parameters
        ----------
        dataset : str
            Dataset
        meta_pd : pd.DataFrame
            Dataset metadata table.
        params_pd : pd.DataFrame
            Combinations of parameters (rows)
        mess : set
            Messages to print
        """
        examples = []
        valid_params = []
        nsams = meta_pd.shape[0]
        meta_cols = meta_pd.columns
        for p, params in params_pd.iterrows():
            train = params['train']
            if train.replace('.', '').isdigit():
                if float(train) < 0.1:
                    valid_params.append(p)
                    m = '[songbird] "%s": %s %s for training is too low' % (
                        dataset, float(train) * 100, '%')
                    self.print_message_or_not(mess, m)
                elif float(train) > 0.95:
                    valid_params.append(p)
                    m = '[songbird] "%s": %s %s for training is too high' % (
                        dataset, float(train) * 100, '%')
                    self.print_message_or_not(mess, m)
                else:
                    examples.append(int(nsams * (1 - float(train))))
            else:
                if train not in set(meta_cols):
                    valid_params.append(p)
                    m = '[songbird] Training column "%s" not in metadata' % (
                        train)
                    self.print_message_or_not(mess, m)
                else:
                    train_vc = meta_pd[train].value_counts()
                    if {'Train', 'Test'}.issubset(set(train_vc.index)):
                        ntrain = train_vc['Train']
                        if nsams < (1.2 * ntrain):
                            valid_params.append(p)
                            m = '[songbird] "%s": %s samples for %s training ' \
                                'samples:' % (dataset, nsams, ntrain)
                            self.print_message_or_not(mess, m)
                    else:
                        valid_params.append(p)
                        m = '[songbird] "%s": no TrainTest in column "%s"' % (
                            dataset, train)
                        self.print_message_or_not(mess, m)
        if valid_params:
            params_pd.drop(index=valid_params, inplace=True)
        if examples:
            params_pd['examples'] = examples

    @staticmethod
    def get_filt_params(params):
        """
        Parameters
        ----------
        params : pd.Series

        Returns
        -------
        filt_list : list
        params_list : list
        """
        filt_list = [
            ('--min-feature-count', str(params['thresh_feats'])),
            ('--min-sample-count', str(params['thresh_samples']))]
        params_list = [
            ('--p-batch-size', str(params['batches'])),
            ('--p-learning-rate', str(params['learns'])),
            ('--p-epochs', str(params['epochs'])),
            ('--differential-prior: %s' % str(params['diff_priors']),
                str(params['diff_priors']).replace('.', '')),
            ('--p-training-column: %s' % str(params['train']),
                str(params['train']).replace('.', '')),
            ('--p-summary-interval: %s' % str(params['summary_interval']),
                str(params['summary_interval']).replace('.', ''))]

        return filt_list, params_list

    @staticmethod
    def get_params_dir(params):
        """

        Parameters
        ----------
        params : pd.Series

        Returns
        -------
        params_dir : str
        """
        params_dir = 'filt_f%s_s%s/%s_%s_%s_%s_%s_%s' % (
            str(params['thresh_feats']),
            str(params['thresh_samples']),
            str(params['batches']),
            str(params['learns']),
            str(params['epochs']),
            str(params['diff_priors'].replace('.', '')),
            str(params['train'].replace('.', '')),
            str(params['summary_interval'].replace('.', ''))
        )
        return params_dir

    @staticmethod
    def get_out(
            odir: str,
            model_null: str) -> tuple:
        """
        Parameters
        ----------
        odir : str
            Output dierctory for a mmvec model/null pair.
        model_null : str
            "model" or null""

        Returns
        -------
        mod_nul_dir : str
        mod_nul_rnk : str
        mod_nul_ord : str
        mod_nul_stt : str
        """
        mod_nul_dir = '%s/%s' % (odir, model_null)
        if not isdir(mod_nul_dir):
            os.makedirs(mod_nul_dir)
        mod_nul_rnk = '%s/ranks.tsv' % mod_nul_dir
        mod_nul_ord = '%s/ordination.txt' % mod_nul_dir
        mod_nul_stt = '%s/stats.qza' % mod_nul_dir
        return mod_nul_dir, mod_nul_rnk, mod_nul_ord, mod_nul_stt

    @staticmethod
    def get_dat_pair_dir(dat, pair):
        if pair:
            return '%s/%s' % (dat, pair)
        else:
            return '%s/unpaired' % dat

    def get_dirs(
            self,
            pair_dir,
            filt,
            subset,
            filt_list,
            params_list,
            model
    ) -> tuple:

        dat_dir = pair_dir
        p_dir = ''
        for (name, level) in [
            ('filter', filt),
            ('subset', subset),
            ('songbird_filt', filt_list),
            ('params', params_list),
            ('model', model)
        ]:
            text = ''
            add_dir = ''
            if name == 'filter':
                add_dir = level
                text += 'Preliminary feature filtering:\n\n'
                text += '%s: <prevalence>_<abundance> thresholds\n' % filt
                text += '\nValue between 0 and <1 indicates a fraction:\n'
                text += '- prevalence: min. fraction of sample presences\n'
                text += '- abundance: min. fraction of samples reads\n'
                text += '  * e.g. "0.1" corresponds to 10 percent\n'
                text += '\nValue >=1 indicates an absolute number:\n'
                text += '- prevalence: min number of sample presences\n'
                text += '- abundance: min number of reads (per sample)\n'
                text += 'In both cases, filtering on prevalence occurs on '
                text += 'per-sample, abundance-filtered features, i.e.:\n'
                text += '\n`tab.loc[(tab_perc > abund).sum(1) > preval]`\n'
            elif name == 'subset':
                add_dir = level
                text += 'Sample subset:\n\n'
                if subset == 'ALL':
                    text += '%s: No sample subset\n' % subset
                else:
                    text += '%s: <variable>_<factor(s)>:\n' % subset
                text += '\n(see your config for formula of this model)\n'
            elif name == 'songbird_filt':
                text += 'Feature and sample filtering in songbird:\n\n'
                for f in filt_list:
                    text += '  %s\n' % ' '.join(f)
                text += '\n(see songbird command line usage)\n'
                add_dir = 'filt_f%s_s%s' % tuple([x[1] for x in filt_list])
                p_dir += add_dir
            elif name == 'params':
                add_dir = '%s_%s_%s_%s_%s_%s' % tuple(
                    [x[1] for x in params_list])
                p_dir += '/' + add_dir
                text += 'Songbird run parameters:\n\n'
                for param_list in params_list:
                    text += '  %s\n' % ' = '.join(list(param_list))
                text += '\n(see songbird command line usage)\n'
            elif name == 'model':
                add_dir = level
                text += 'Model: %s\n' % model
                text += '\n(see your config for formula of this model)\n'

            dat_dir = '%s/%s' % (dat_dir, add_dir)
            o_dir = self.get_output(dat_dir)
            readme = '%s/readme.txt' % o_dir
            with open(rep(readme), 'w') as o:
                o.write(text)

        new_qza = '%s/tab.qza' % o_dir
        new_meta = '%s/metadata.tsv' % o_dir
        return dat_dir, p_dir, o_dir, new_qza, new_meta

    # @staticmethod
    # def get_main_dirs(
    #         pair_dir, filt, subset, params_dir, model, config) -> tuple:
    #     datdir = '%s/%s/%s/%s/%s' % (pair_dir, filt, subset,
    #                                  params_dir, model)
    #     odir = get_output(config.folder,
    #                                'songbird/%s' % datdir)
    #     new_qza = '%s/tab.qza' % odir
    #     new_meta = '%s/metadata.tsv' % odir
    #     return datdir, odir, new_qza, new_meta

    def get_out_paths(self, o_dir, model_baseline, baselines) -> dict:
        if model_baseline in baselines:
            b_diff_qza = ''
            b_stat = baselines[model_baseline]
            b_plot = ''
        else:
            b_diff_qza = '%s/differentials-baseline.qza' % self.out
            b_stat = '%s/differentials-stats-baseline.qza' % self.out
            b_plot = '%s/differentials-biplot-baseline.qza' % self.out
            baselines[model_baseline] = b_stat
        out_paths = {
            'diff': '%s/differentials.tsv' % o_dir,
            'diff_qza': '%s/differentials.qza' % o_dir,
            'stat': '%s/differentials-stats.qza' % o_dir,
            'plot': '%s/differentials-biplot.qza' % o_dir,
            'tens': '%s/tensorboard.qzv' % self.out,
            'html': '%s/tensorboard.html' % self.out,
            'bdiff_qza': b_diff_qza,
            'bstat': b_stat,
            'bplot': b_plot
        }
        return out_paths

    @staticmethod
    def write_new_meta(meta_pd, new_meta, meta_vars, drop, params):
        meta_cols = set(meta_pd.columns)
        if params['train'] in meta_cols:
            meta_vars.add(params['train'])
        new_meta_pd = meta_pd[
            (['sample_name'] + [x for x in meta_vars if x in meta_cols])
        ].copy()
        new_meta_pd = new_meta_pd.loc[~new_meta_pd.isna().any(1)]
        new_meta_pd = rename_duplicate_columns(new_meta_pd)
        if drop:
            to_remove = pd.concat([
                new_meta_pd[meta_var].isin(var_drop)
                # new_meta_pd[meta_var.lower()].isin(var_drop)
                for meta_var, var_drop in drop.items()
            ], axis=1).any(axis=1)
            new_meta_pd = new_meta_pd.loc[~to_remove]
        new_meta_pd.to_csv(rep(new_meta), index=False, sep='\t')
        return new_meta_pd.shape[0]

    def summarize_songbirds(self):
        q2s = []
        songbird = self.get_output()
        for root, dirs, files in os.walk(rep(songbird)):
            for fil in files:
                if fil == 'tensorboard.html':
                    path = root + '/' + fil
                    diff = '%s/differentials.tsv' % dirname(root)
                    root_split = root.split('%s/' % songbird)[-1].split('/')
                    d, pr, fr, sb, sr, ps, ml, be = root_split
                    with open(rep(path)) as f:
                        for line in f:
                            if 'Pseudo Q-squared' in line:
                                ls = line.split(
                                    'Pseudo Q-squared:</a></strong> ')
                                q2s.append([
                                    pr, d, fr, sb, ml, sr, ps, be, diff,
                                    float(ls[-1].split('<')[0])])
        if q2s:
            self.q2s_pd = pd.DataFrame(q2s, columns=[
                'pair', 'dataset', 'filter', 'subset', 'model',
                'songbird_filter', 'parameters', 'baseline', 'differentials',
                'Pseudo_Q_squared'])
            q2s_fp = '%s/songbird_q2.tsv' % songbird
            self.q2s_pd.to_csv(q2s_fp, index=False, sep='\t')
            print('Written -> %s' % q2s_fp)

    def create_songbird_feature_metadata(self):
        if self.q2s_pd.shape[0]:
            q2_pd = self.q2s_pd.copy()
            # q2_pd = self.q2s_pd.loc[(self.q2s_pd.pair == 'unpaired') &
            #                         (self.q2s_pd.Pseudo_Q_squared > 0)]
            for dat, dataset_pd in q2_pd.groupby('dataset'):
                dataset_sbs = []
                for r, row in dataset_pd.iterrows():
                    pr = 'pair=%s' % row['pair']
                    fr = 'filter=%s' % row['filter']
                    sb = 'subset=%s' % row['subset']
                    ml = 'model=%s' % row['model']
                    st = 'sb_filt=%s' % row['songbird_filter']
                    ps = 'params=%s' % row['parameters']
                    be = 'baseline=%s' % row['baseline']
                    q2 = '[Q2=%s]' % row['Pseudo_Q_squared']
                    diffs = row['differentials']
                    sb_pd = pd.read_csv(rep(diffs), index_col=0, sep='\t')
                    sb_pd.columns = ['%s %s: %s' % (
                        '__'.join([dat, pr, fr, sb, ml, st, ps, be]),
                        q2, x) for x in sb_pd.columns]
                    dataset_sbs.append(sb_pd)
                if len(dataset_sbs):
                    dataset_sbs_pd = pd.concat(dataset_sbs, axis=1, sort=False)
                    self.get_output(dat)
                    fpo_tsv = '%s/differentials_%s.tsv' % (self.out, dat)
                    self.project.datasets[dat].sb = fpo_tsv
                    fpo_qza = '%s/differentials_%s.qza' % (self.out, dat)
                    dataset_sbs_pd.index.name = 'Feature ID'
                    dataset_sbs_pd.to_csv(rep(fpo_tsv), sep='\t')
                    run_import(fpo_tsv, fpo_qza, 'FeatureData[Differential]')
                    io_update(self, i_f=fpo_tsv, o_f=fpo_qza, key=(dat, ''))

    def get_songbird_pd(self, songbird):
        self.songbird_pd = pd.DataFrame(
            songbird, columns=[
                'dataset', 'qza', 'meta', 'filter', 'params', 'subset',
                'differentials', 'baseline', 'html', 'pair'])

    def check_metadata_models(self, meta_fp, meta_pd, songbird_models):
        models = {}
        for model, formula_ in songbird_models.items():
            print(model, formula_)
            drop = {}
            levels = {}
            variables = set()
            formula = formula_.strip('"').strip("'")
            if formula.startswith('C('):
                formula_split = formula.split('C(')[-1].rsplit(')', 1)
                formula_split_c = formula_split[0].split(',')[0].strip().strip()
                formula = 'C(%s)' % formula_split[0].replace(
                    formula_split_c, formula_split_c)
                formula += formula_split[1]
                if 'Diff' in formula:
                    levels = {formula_split_c: [
                        x.strip().strip('"').strip("'")
                        for x in formula.split(
                            "levels=['")[-1].split("']")[0].split(",")]}
                elif "Treatment('" in formula:
                    levels = {formula_split_c: [
                        formula.split("Treatment('")[-1].split("')")[0]]}
                elif 'Treatment("' in formula:
                    levels = {formula_split_c: [
                        formula.split('Treatment("')[-1].split('")')[0]]}
                variables.add(formula_split_c)
                variables.update(set([x for x in re.split(
                    '[+/:*]', formula_split[1]) if x]))
            else:
                formula_split = re.split('[+/:*]', formula)
                formula = formula
                variables.update(set([x for x in formula_split]))

            common_with_md = set(meta_pd.columns.values) & variables
            if sorted(variables) != sorted(common_with_md):
                only_formula = sorted(variables ^ common_with_md)
                issue = 'Songbird formula term(s) missing in metadata:\n\t' \
                        '%s\n\t  [not used]: %s=%s' % (
                            ', '.join(only_formula), model, formula)
                self.models_issues.setdefault(issue, set()).add(meta_fp)
                continue

            if levels:
                levels_set = sorted([x for x in meta_pd[
                    formula_split_c].unique() if str(x) != 'nan'])
                if 'Diff' in formula:
                    cur_levels = levels[formula_split_c]
                    common_levels = set(levels_set) & set(cur_levels)
                    only_meta = set(levels_set) ^ common_levels
                    only_model = set(cur_levels) ^ common_levels
                    if len(only_model):
                        issue = 'Songbird formula "Diff" factors(s) missing' \
                                ' in metadata "%s": %s' % (
                                    formula_split_c, list(only_model))
                        self.models_issues.setdefault(issue, set()).add(meta_fp)
                        continue
                    if len(only_meta):
                        drop[formula_split_c] = list(only_meta)
                        issue = 'Songbird formula "Diff" factors(s) ' \
                                'incomplete for metadata "%s":\n' \
                                '\t -> skipping samples with %s' % (
                                    formula_split_c, list(only_meta))
                        self.models_issues.setdefault(issue, set()).add(meta_fp)
                elif 'Treatment(' in formula:
                    levels = {formula_split_c: formula.split(
                        "Treatment('")[-1].split("')")[0]}
                    if levels[formula_split_c] not in levels_set:
                        issue = 'Songbird formula "Treatment" factors(s)' \
                                ' missing in metadata "%s" [%s]' % (
                                    formula_split_c, levels)
                        self.models_issues.setdefault(issue, set()).add(meta_fp)
                        continue
            models[model] = [formula, variables, drop]
        return models

    def show_models_issues(self, mess):
        if mess:
            for m in sorted(mess):
                print(m)
        if self.models_issues:
            print('\n## Issues with model (will not run)')
            for model_issue, metas in self.models_issues.items():
                print('# -', model_issue)
                for meta in sorted(metas):
                    print('#\t.%s' % meta.replace(dirname(self.dir), ''))
            print('\n')

    def make_qurros(self) -> None:
        """Make qurro plots"""
        self.analysis = 'qurro'
        for r, row in self.songbird_pd.iterrows():
            dat = row['dataset']
            tax = self.project.datasets[dat].tax[-1]
            qurro_qzv = '%s_qurro.qzv' % splitext(row['differentials'])[0]
            is_qurro = to_do(qurro_qzv)
            is_diff = to_do(row['differentials'])
            if self.config.force or is_qurro and not is_diff:
                qza = row['qza']
                meta = row['meta']
                diffs = '%s.qza' % splitext(row['differentials'])[0]
                cmd = 'qiime qurro differential-plot'
                cmd += ' --i-table %s' % qza
                cmd += ' --i-ranks %s' % diffs
                cmd += ' --m-sample-metadata-file %s' % meta
                cmd += ' --m-feature-metadata-file %s' % tax
                cmd += ' --o-visualization %s' % qurro_qzv
                self.cmds.setdefault(dat, []).append(cmd)
                io_update(
                    self, i_f=[qza, meta, diffs], o_f=qurro_qzv, key=(dat, ''))
        self.register_io_command()

    def songbird(self) -> None:
        """Main script for the creation of songbird jobs.
        It iterates over the rows of the table created
        upfront and over each combination of parameters
        and collect the output info for potential reuse
        in figure generation and post-analysis.
        """
        mess = set()
        songbird = []
        self.analysis = 'songbird'
        params_pd = self.get_params_combinations()
        for r, row in self.songbirds.iterrows():
            qza, pair, meta_fp = row['qza'], row['pair'], row['meta']
            dat, filt, subset = row['dataset'], row['filter'], row['subset']
            if dat not in self.songbird_models:
                continue
            pair_dir = self.get_dat_pair_dir(dat, pair)
            meta_pd = read_meta_pd(rep(meta_fp))
            models = self.check_metadata_models(
                meta_fp, meta_pd, self.songbird_models[dat])
            row_params_pd = params_pd.copy()
            self.process_params_combinations(dat, meta_pd, row_params_pd, mess)
            for p, params in row_params_pd.iterrows():
                filt_list, params_list = self.get_filt_params(params)
                baselines, model_baselines = {}, {'1': '1'}
                for modx, model in enumerate(models.keys()):
                    formula, meta_vars, drop = models[model]
                    dat_dir, p_dir, o_dir, new_qza, new_meta = self.get_dirs(
                        pair_dir, filt, subset, filt_list, params_list, model)
                    nsams = self.write_new_meta(
                        meta_pd, new_meta, meta_vars, drop, params)
                    if dat in self.baselines and model in self.baselines[dat]:
                        model_baselines = self.baselines[dat][model]
                    for model_base in model_baselines:
                        b_formula = model_baselines[model_base]
                        self.get_output('%s/b-%s' % (dat_dir, model_base))
                        out_paths = self.get_out_paths(
                            o_dir, model_base, baselines)
                        # convergence = self.check_stats_convergence(out_paths)
                        cmd, fcmd, bcmd = write_songbird(
                            self, dat, qza, new_qza, new_meta, nsams, params,
                            formula, b_formula, out_paths)
                        songbird.append([
                            dat, new_qza, meta_fp, filt, '%s_%s' % (
                                p_dir.replace('/', '__'), model),
                            subset, out_paths['diff'], model_base,
                            out_paths['html'], pair])
                        if cmd:
                            self.cmds.setdefault(dat, []).append(cmd)
                        if fcmd:
                            self.fcmds.setdefault(dat, []).append(fcmd)
                        if bcmd:
                            self.bcmds.setdefault(dat, []).append(bcmd)
        if songbird:
            self.get_songbird_pd(songbird)
        self.show_models_issues(mess)

        self.analysis = 'songbird_filter'
        self.register_io_command()
        self.analysis = 'songbird_baselines'
        self.register_io_command()
        self.analysis = 'songbird'
        self.register_io_command()

        self.summarize_songbirds()
        self.create_songbird_feature_metadata()

    def register_io_command(self):
        if self.analysis == 'songbird_baselines':
            io = 'b'
            AnalysisPrep.analyses_commands[self.analysis] = dict(self.bcmds)
        elif self.analysis == 'songbird_filter':
            io = 'f'
            AnalysisPrep.analyses_commands[self.analysis] = dict(self.fcmds)
        else:
            io = ''
            AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)

        ios = dict((x, {i[0]: j for i, j in y.items() if i[1] == io})
                   for x, y in self.ios.items())
        AnalysisPrep.analyses_ios[self.analysis] = ios
        self.ios = {}
        self.cmds = {}
