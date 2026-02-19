# ----------------------------------------------------------------------------
# Copyright (c) 2026, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import itertools
import pandas as pd
from os.path import dirname, isdir, splitext
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import run_import, write_songbird
from microbiome_analyzer.core.metadata import (
    get_subset_pd, rename_duplicate_columns, make_random_train_test,
    get_train_perc_from_numeric, get_cat_vars_and_vc, make_train_test_from_cat)
from microbiome_analyzer.analyses.models import parse_formula
from microbiome_analyzer.analyses.filter import (
    filter_mb_table, filter_non_mb_table, skip_subset)
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
        self.model = ''
        self.formula = ''
        if config.diff_abund:
            (self.songbird_models, self.filtering, self.params, self.baselines,
             self.songbird_subsets) = self.get_songbird_dicts()
            self.models = {}
            self.models_issues = {'var': {}, 'fac': {}, 'exp': {}}
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
        if 'models' not in self.config.diff_abund:
            print('No models in %s' % self.config.diff_abund)
            sys.exit(0)
        is_mb = {}
        models = {}
        for dat in self.config.diff_abund['models'].keys():
            if dat[-1] == '*':
                is_mb[dat[:-1]] = (dat[:-1], 1)
                models[dat[:-1]] = self.config.diff_abund['models'][dat]
            else:
                is_mb[dat] = (dat, 0)
                models[dat] = self.config.diff_abund['models'][dat]
        return models, is_mb

    def get_songbird_baselines(self):
        baselines = {}
        if 'baselines' in self.config.diff_abund:
            return self.config.diff_abund['baselines']
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
        if 'params' not in self.config.diff_abund:
            print('No parameters set in %s:\nUsing defaults: %s' % (
                self.config.diff_abund,
                ', '.join(['%s: %s' % (k, v) for k, v in params.items()])))
        else:
            for param in self.config.diff_abund['params']:
                cur_param = self.config.diff_abund['params'][param]
                if not isinstance(cur_param, list):
                    print('Parameter %s should be a list (correct in %s)\n' % (
                        param, self.config.diff_abund))
                    sys.exit(0)
                params[param] = cur_param
        return params

    def get_filtering(self, models, is_mb):
        filtering = {
            '0_0': dict((is_mb[dat], ['0', '0']) for dat in models.keys())}
        if 'filtering' not in self.config.diff_abund:
            pass
        else:
            if 'global' in self.config.diff_abund['filtering']:
                for fname, p_a in self.config.diff_abund[
                    'filtering']['global'].items():
                    for dat in models.keys():
                        if fname not in filtering:
                            filtering[fname] = {}
                        filtering[fname][is_mb[dat]] = p_a
            for dat, filts in self.config.diff_abund['filtering'].items():
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
        if 'subsets' in self.config.diff_abund:
            subsets.update(self.config.diff_abund['subsets'])
        return models, filtering, params, baselines, subsets

    def merge_subsets_apply(self):
        subsets_fp = [[dataset, name, subsets, '']
                      for name, subsets in self.songbird_subsets.items()
                      for dataset in self.songbirds.dataset.unique()]
        if subsets_fp:
            subsets = pd.DataFrame(
                subsets_fp, columns=['dataset', 'subset', 'subsets', 'pair'])
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
        for rdx, row in enumerate(mmvec_pd.values):
            pair, _, subset, dat1, dat2, prev1, abun1, prev2, abun2 = row[:9]
            meta_fp, tsv1, tsv2, qza1, qza2 = row[10:15]
            filt1 = '%s-%s' % (prev1, abun1)
            filt2 = '%s-%s' % (prev2, abun2)
            mmvecs.append([
                dat1, filt1, prev1, abun1, subset, pair, tsv1, qza1, meta_fp])
            mmvecs.append([
                dat2, filt2, prev2, abun2, subset, pair, tsv2, qza2, meta_fp])
        if mmvecs and self.songbirds.shape[0]:
            self.songbirds.drop(columns=['is_mb', 'subsets'], inplace=True)
            mmvecs_pd = pd.DataFrame(mmvecs, columns=self.songbirds.columns)
            self.songbirds = pd.concat([self.songbirds, mmvecs_pd])

    def make_datasets_paths(self):
        self.get_datasets_paths()
        self.analysis = 'songbird_imports'
        if self.songbirds.shape[0]:
            for (dat, _, subset), row in self.songbirds.groupby(
                    ['dataset', 'filter', 'subset']):
                row_d = row.iloc[0, :].to_dict()
                tsv, qza, meta = row_d['tsv'], row_d['qza'], row_d['meta']
                data = self.project.datasets[dat]
                subsets = row_d['subsets']
                meta_pd = get_subset_pd(data.metadata, dat, subset, subsets)
                meta_pd.to_csv(rep(meta), index=False, sep='\t')
                if not self.config.force and not (to_do(tsv) or to_do(qza)):
                    continue
                if not data.data:
                    continue
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
                for analysis in ['songbird']:
                    self.analysis = analysis
                    self.get_output('datasets/%s/%s' % (dataset, subset))
                    rad = '%s_%s' % (dataset.replace('/', '_'), filter_)
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
                    new_meta_pd, cat_pd, vc, train_perc, meta_fp, cat_vars,
                    train_col, rep_d)
            else:
                train_samples = make_random_train_test(new_meta_pd, train_perc)
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
                vars_pd = vars_pd.loc[~vars_pd.isna().any(axis=1)]
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
    def print_message_or_not(messages, mess):
        if mess not in messages:
            messages.add(mess)

    def get_train_example(
            self,
            meta_pd: pd.DataFrame,
            params: pd.Series):
        """Get the number of examples or the "Train/Test" metadata variable
        to use for training, or if the given value for any is not allow, give
        an explicit error message to print before skipping.

        Parameters
        ----------
        meta_pd : pd.DataFrame
            Dataset metadata table.
        params : pd.Series
            Combinations of parameters (rows)

        Returns
        -------
        mess : str
            Error message indicating that the model is skipped (none by default)
        """
        ns = meta_pd.shape[0]
        cols = set(meta_pd.columns)
        mess = ''
        t = params['train']
        if t.replace('.', '').isdigit():
            if float(t) < 0.1:
                mess = '%s %s for training is too low' % (float(t) * 100, '%')
            elif float(t) > 0.95:
                mess = '%s %s for training is too high' % (float(t) * 100, '%')
            else:
                params['examples'] = int(ns * (1 - float(t)))
        else:
            if t not in set(cols):
                mess = 'Training column "%s" not in metadata' % t
            else:
                vc = meta_pd[t].value_counts()
                if {'Train', 'Test'}.issubset(set(vc.index)):
                    nt = vc['Train']
                    if ns < (1.2 * nt):
                        mess = '%s samples for %s training samples:' % (ns, nt)
                else:
                    mess = 'No TrainTest in column "%s"' % t
        return mess

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
            Output directory for a mmvec model/null pair.
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

        return dat_dir, p_dir, o_dir

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
    def new_meta(meta_pd, new_meta, variables, params):
        meta_cols = set(meta_pd.columns)
        if params['train'] in meta_cols:
            variables[params['train']] = []
        new_meta_pd = meta_pd[
            (['sample_name'] + [v for v in variables if v in meta_cols])].copy()
        new_meta_pd = new_meta_pd.loc[~new_meta_pd.isna().any(axis=1)]
        new_meta_pd = rename_duplicate_columns(new_meta_pd)
        # if drop:
        #     to_remove = pd.concat([
        #         new_meta_pd[meta_var].isin(var_drop)
        #         for meta_var, var_drop in drop.items()
        #     ], axis=1).any(axis=1)
        #     new_meta_pd = new_meta_pd.loc[~to_remove]
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
                    root_s = root.split('%s/' % rep(songbird))[-1].split('/')
                    if len(root_s) == 8:
                        d, pr, fr, sb, sr, ps, ml, be = root_s
                    elif len(root_s) == 9:
                        d, sstx, pr, fr, sb, sr, ps, ml, be = root_s
                        d = d + '/' + sstx
                    else:
                        d, ss, tx, pr, fr, sb, sr, ps, ml, be = root_s
                        d = d + '/' + ss + '/' + tx
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
                'songbird_filter', 'parameters', 'baseline',
                'differentials', 'Pseudo_Q_squared'])
            q2s_fp = '%s/songbird_q2.tsv' % songbird
            self.q2s_pd.to_csv(rep(q2s_fp), index=False, sep='\t')
            print('Written -> %s' % q2s_fp)

    def create_songbird_feature_metadata(self):
        if self.q2s_pd.shape[0]:
            q2_pd = self.q2s_pd.copy()
            for dat, dataset_pd in q2_pd.groupby('dataset'):
                if dat not in self.project.datasets:
                    continue
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
                    diffs_qza = '%s.qza' % splitext(diffs)[0]
                    if to_do(diffs) and not to_do(diffs_qza):
                        self.models_issues['exp'][diffs] = diffs_qza
                        continue
                    sb_pd = pd.read_csv(rep(diffs), sep='\t')
                    sb_pd = sb_pd.set_index('featureid')
                    if sb_pd.index[0][0]=='#':
                        sb_pd = sb_pd.iloc[1:, :]
                    sb_pd = sb_pd.astype(float)
                    sb_pd.columns = ['%s %s: %s' % (
                        '__'.join([dat, pr, fr, sb, ml, st, ps, be]),
                        q2, x) for x in sb_pd.columns]
                    dataset_sbs.append(sb_pd)
                if len(dataset_sbs):
                    sbs_pd = pd.concat(dataset_sbs, axis=1, sort=False)
                    self.get_output(dat)
                    tsv = '%s/differentials.tsv' % self.out
                    self.project.datasets[dat].sb = tsv
                    qza = '%s/differentials.qza' % self.out
                    sbs_pd.index.name = 'Feature ID'
                    sbs_pd = sbs_pd[sorted(sbs_pd.columns)]
                    sbs_pd.to_csv(rep(tsv), sep='\t')
                    cmd = run_import(tsv, qza, 'FeatureData[Differential]')
                    io_update(self, i_f=tsv, o_f=qza, key=(dat, ''))

    def get_songbird_pd(self, songbird):
        self.songbird_pd = pd.DataFrame(
            songbird, columns=[
                'dataset', 'qza', 'meta', 'filter', 'params', 'subset',
                'differentials', 'baseline', 'html', 'pair'])

    def check_metadata_variables(
            self,
            variables,
            meta_pd,
            meta_fp: str
    ):
        skip = False
        not_in_meta = set(variables).difference(set(meta_pd.columns.values))
        if not_in_meta:
            print(set(variables))
            print(set(meta_pd.columns.values))
            print(not_in_meta)
            print(self.model)
            print(self.formula)
            self.models_issues['var'].setdefault(meta_fp, []).append(
                [not_in_meta, self.model, self.formula])
            skip = True
        return skip

    def check_metadata_factors(
            self,
            variables,
            meta_pd,
            meta_fp: str
    ):
        skip = False
        for v, factors in variables.items():
            if not factors:
                continue
            not_in_v = set(factors).difference(set(meta_pd[v]))
            if not_in_v:
                # self.models_issues['fac'].setdefault(meta_fp, []).append(
                #     [not_in_v, self.model, self.formula, v])
                self.models_issues['fac'].setdefault(meta_fp, []).append(
                    [not_in_v, self.model, self.formula, v])
                skip = True
        return skip

    def check_metadata_models(self, meta_fp, meta_pd, songbird_models):
        models = {}
        for model, formula in songbird_models.items():
            self.model = model
            self.formula = formula
            variables = parse_formula(formula)
            if self.check_metadata_variables(variables, meta_pd, meta_fp):
                continue
            if self.check_metadata_factors(variables, meta_pd, meta_fp):
                continue
            models[model] = (formula, variables)
        return models

    def show_models_issues(self, messages):
        if messages:
            for mess in sorted(messages):
                print(mess)

        if self.models_issues['var']:
            print('\n## [songbird] Variables(s) missing in metadata:')
            for fp, (not_in, mod, form) in self.models_issues['var'].items():
                print("")
                print("##-----file: %s" % rep(fp))
                print("##     model: %s -> %s" % (mod, form))
                print("##     missing: %s" % ', '.join(sorted(not_in)))
            print('\n')

        if self.models_issues['fac']:
            print('\n## [songbird] Factor(s) missing in metadata:')
            for fp, issues in self.models_issues['fac'].items():
                print("")
                print("##-----file: %s" % rep(fp))
                for (not_in, mod, form, v) in issues:
                    print("")
                    print("##-----file: %s" % rep(fp))
                    print("##     model: %s -> %s" % (mod, form))
                    print("##     variable: %s" % v)
                    print("##     missing: %s" % ', '.join(sorted(not_in)))
            print('\n')

        if self.models_issues['exp']:
            print('\n## [songbird] Export(s) not completed:')
            for tsv, qza in self.models_issues['exp'].items():
                print("##     qza: %s" % rep(qza))
                print("##     tsv: %s" % rep(tsv))
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
        """
        Main script for the creation of songbird jobs. It iterates over  the
        rows of the table created upfront and over each combination of
        parameters and collect the output info for potential reuse in figure
        generation and post-analysis.
        """
        messages = set()
        songbird = []
        self.analysis = 'songbird'
        params_pd = self.get_params_combinations()
        for r, row in self.songbirds.iterrows():
            qza, pair, meta_fp = row['qza'], row['pair'], row['meta']
            dat, filt, subset = row['dataset'], row['filter'], row['subset']
            if dat not in self.songbird_models or to_do(qza):
                continue
            pair_dir = self.get_dat_pair_dir(dat, pair)
            meta_pd = read_meta_pd(rep(meta_fp))
            models = self.check_metadata_models(
                meta_fp, meta_pd, self.songbird_models[dat])
            for p, params in params_pd.iterrows():
                mess = self.get_train_example(meta_pd, params)
                if mess:
                    self.print_message_or_not(messages, mess)
                    continue
                filt_list, params_list = self.get_filt_params(params)
                baselines, model_baselines = {}, {'1': '1'}
                for model, (formula, variables) in models.items():
                    dat_dir, p_dir, o_dir = self.get_dirs(
                        pair_dir, filt, subset, filt_list, params_list, model)
                    new_qza = '%s/tab.qza' % o_dir
                    new_meta = '%s/metadata.tsv' % o_dir
                    if skip_subset(variables, meta_pd):
                        continue
                    nsams = self.new_meta(meta_pd, new_meta, variables, params)
                    if dat in self.baselines and model in self.baselines[dat]:
                        if self.baselines[dat][model]:
                            model_baselines = self.baselines[dat][model]
                    for model_base in model_baselines:
                        b_formula = model_baselines[model_base]
                        self.get_output('%s/b-%s' % (dat_dir, model_base))
                        out_paths = self.get_out_paths(
                            o_dir, model_base, baselines)
                        # convergence = self.check_stats_convergence(out_paths)
                        cmd, fcmd, bcmd, skip = write_songbird(
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

        self.analysis = 'songbird_filter'
        self.register_io_command()
        self.analysis = 'songbird_baselines'
        self.register_io_command()
        self.analysis = 'songbird'
        self.register_io_command()

        self.summarize_songbirds()
        self.create_songbird_feature_metadata()
        self.show_models_issues(messages)

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
        if self.analysis not in ['songbird_baselines', 'songbird_filter']:
            self.cmds = {}
