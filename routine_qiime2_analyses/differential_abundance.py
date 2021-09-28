# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import random
import datetime
import itertools
import pandas as pd
from os.path import dirname, isdir, isfile, splitext
from routine_qiime2_analyses.analyses_prep import AnalysisPrep
from routine_qiime2_analyses._routine_q2_cmds import (
    get_case, get_new_meta_pd, run_import, songbird_cmd,
)
from routine_qiime2_analyses._routine_q2_metadata import (
    rename_duplicate_columns,  get_train_perc_from_numeric,
    get_cat_vars_and_vc, make_train_test_from_cat
)
from routine_qiime2_analyses._routine_q2_io_utils import (
    read_meta_pd, get_analysis_folder, filter_mb_table, filter_non_mb_table,
    write_filtered_tsv
)
from routine_qiime2_analyses._routine_q2_songbird import (
    get_songbird_dicts, get_unique_filterings
)


class DiffModels(object):

    def __init__(self, config, project) -> None:
        self.cmds = {}
        self.config = config
        self.project = project
        if config.diff_models:
            songbird_dicts = get_songbird_dicts(config.diff_models)
            self.songbird_models = songbird_dicts[0]
            self.songbird_filtering = songbird_dicts[1]
            self.unique_filtering = get_unique_filterings(songbird_dicts[1])
            self.params = songbird_dicts[2]
            self.models_baselines = songbird_dicts[3]
            self.models = {}
            self.models_issues = {}
            self.songbird_datasets = songbird_dicts[4]
            self.songbird_subsets = songbird_dicts[5]
            self.songbirds = pd.DataFrame(columns=[
                'dataset', 'is_mb', 'filter', 'prevalence', 'abundance'])
            self.params_list = [
                'train', 'batches', 'learns', 'epochs', 'diff_priors',
                'thresh_feats', 'thresh_samples', 'summary_interval']
            self.q2s_pd = pd.DataFrame()
        self.songbird_pd = pd.DataFrame()

    def merge_subsets_apply(self):
        subsets_fp = [
            [dataset, var, subset, get_case(subset, var), '']
            for var, subsets in self.songbird_subsets.items()
            for subset in subsets
            for dataset in self.songbirds.dataset.unique()]
        if subsets_fp:
            subsets = pd.DataFrame(
                subsets_fp, columns=['dataset', 'variable', 'factors',
                                     'subset', 'pair'])
            self.songbirds = self.songbirds.merge(
                subsets, on=['dataset'], how='outer')

    def get_songbirds_filts(self, project):
        filts_df = []
        for (dat, is_mb), filts_dats in self.unique_filtering.items():
            if dat not in project.datasets:
                continue
            for (filt, prev, abund) in filts_dats:
                filts_df.append([dat, is_mb, filt, prev, abund])
        if filts_df:
            self.songbirds = pd.DataFrame(filts_df, columns=[
                'dataset', 'is_mb', 'filter', 'prevalence', 'abundance'])

    def prep_songbirds(self, mmvec_pd, project):
        self.get_songbirds_filts(project)
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
        cmds = {}
        self.get_datasets_paths()
        if self.songbirds.shape[0]:
            for (dataset, filter, subset), row in self.songbirds.groupby(
                    ['dataset', 'filter', 'subset']):
                row_d = row.iloc[0, :].to_dict()
                tsv, qza, meta = row_d['tsv'], row_d['qza'], row_d['meta']
                if isfile(tsv) and isfile(qza) and isfile(meta):
                    continue
                data = self.project.datasets[dataset]
                variable, factors = row_d['variable'], row_d['factors']
                meta_pd = get_new_meta_pd(
                    data.metadata[0], subset, variable, factors)
                tsv_pd = data.data[0][meta_pd.sample_name.tolist()]
                preval, abund = row_d['prevalence'], row_d['abundance']
                if row_d['is_mb']:
                    tsv_pd, res = filter_mb_table(preval, abund, tsv_pd)
                else:
                    tsv_pd, res = filter_non_mb_table(preval, abund, tsv_pd)
                meta_pd.to_csv(meta, index=False, sep='\t')
                if self.config.force or not isfile(tsv):
                    write_filtered_tsv(tsv, tsv_pd)
                if self.config.force or not isfile(qza):
                    cmd = run_import(tsv, qza, 'FeatureTable[Frequency]')
                    cmds.setdefault(dataset, []).append(cmd)
        self.register_command('songbird_imports', cmds)

    def get_datasets_paths(self):
        paths = []
        if self.songbirds.shape[0]:
            for r, row in self.songbirds.iterrows():
                dataset = row['dataset']
                filter = row['filter']
                subset = row['subset']
                for analysis in ['mmvec', 'songbird']:
                    odir = get_analysis_folder(
                        self.config.i_datasets_folder, '%s/datasets/%s/%s' % (
                            analysis, dataset, subset))
                    rad = '%s_%s' % (dataset, filter)
                    tsv = '%s/tab_%s.tsv' % (odir, rad)
                    qza = '%s.qza' % splitext(tsv)[0]
                    meta = '%s/meta_%s.tsv' % (odir, rad)
                    if isfile(tsv) and isfile(qza) and isfile(meta):
                        paths.append([tsv, qza, meta])
                        break
                    elif analysis == 'songbird':
                        paths.append([tsv, qza, meta])
        if paths:
            self.songbirds = pd.concat([
                self.songbirds,
                pd.DataFrame(paths, columns=['tsv', 'qza', 'meta'])], axis=1)

    @staticmethod
    def get_traintests(meta_fp, new_meta_pd, vars, train, train_col):
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
                               meta_pd, dat) -> dict:
        train_tests = {}
        train = train_test_d['train']
        meta_tt_pd = meta_pd.set_index('sample_name').copy()
        if 'datasets' in train_test_d and dat in train_test_d['datasets']:
            for tt, vars in train_test_d['datasets'][dat].items():
                vars_pd = meta_tt_pd[vars].copy()
                vars_pd = vars_pd.loc[~vars_pd.isna().any(1)]
                vars_pd = rename_duplicate_columns(vars_pd)
                trains = self.get_traintests(meta_fp, vars_pd, vars,
                                             str(train), tt)
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
                meta_subset = read_meta_pd(meta_fp)
                train_tests = self.make_train_test_column(
                    meta_fp, self.config.train_test_dict, meta_subset, dat)
                rewrite = False
                meta_subset_cols = set(meta_subset.columns)
                for train_col, train_samples in train_tests.items():
                    if train_col not in meta_subset_cols:
                        rewrite = True
                        meta_subset[train_col] = [
                            'Train' if x in set(train_samples) else
                            'Test' for x in meta_subset.sample_name.tolist()]
                if self.config.force or rewrite:
                    meta_subset.to_csv(meta_fp, index=False, sep='\t')

    def get_params_combinations(self):
        """Make a pandas data frame from the combinations
        of songbird run/hyper-parameters. It includes the
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
                    m = '\t[skip] "%s": train %s too low (%s)' % (
                        dataset, '%', train)
                    self.print_message_or_not(mess, m)
                elif float(train) > 0.95:
                    valid_params.append(p)
                    m = '\t[skip] "%s": train %s too high (%s)' % (
                        dataset, '%', train)
                    self.print_message_or_not(mess, m)
                else:
                    examples.append(int(nsams * (1 - float(train))))
            else:
                if train not in set(meta_cols):
                    valid_params.append(p)
                    m = '\t[skip] Training column "%s" not in metadata' % (
                        train)
                    self.print_message_or_not(mess, m)
                else:
                    train_vc = meta_pd[train].value_counts()
                    if {'Train', 'Test'}.issubset(set(train_vc.index)):
                        ntrain = train_vc['Train']
                        if nsams < (1.2 * ntrain):
                            valid_params.append(p)
                            m = '\t[skip] "%s": %s samples for %s training ' \
                                'samples:' % (dataset, nsams, ntrain)
                            self.print_message_or_not(mess, m)
                    else:
                        valid_params.append(p)
                        m = '\t[skip] "%s": no TrainTest in column "%s"' % (
                            dataset, train)
                        self.print_message_or_not(mess, m)
        if valid_params:
            params_pd.drop(index=valid_params, inplace=True)
        if examples:
            params_pd['examples'] = examples

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
            dat_pair = '%s_%s' % (dat, pair)
            pair_dir = '%s/%s' % (dat, pair)
        else:
            dat_pair = '%s/unpaired' % dat
            pair_dir = '%s/unpaired' % dat
        return dat_pair, pair_dir

    @staticmethod
    def get_main_dirs(
            pair_dir, filt, subset, params_dir, model, config) -> tuple:
        datdir = '%s/%s/%s/%s/%s' % (pair_dir, filt, subset,
                                     params_dir, model)
        odir = get_analysis_folder(config.i_datasets_folder,
                                   'songbird/%s' % datdir)
        new_qza = '%s/tab.qza' % odir
        new_meta = '%s/metadata.tsv' % odir
        return datdir, odir, new_qza, new_meta

    @staticmethod
    def get_out_paths(odir, bodir, model_baseline, baselines) -> dict:
        if model_baseline in baselines:
            bdiff_qza = ''
            bstat = baselines[model_baseline]
            bplot = ''
        else:
            bdiff_qza = '%s/differentials-baseline.qza' % bodir
            bstat = '%s/differentials-stats-baseline.qza' % bodir
            bplot = '%s/differentials-biplot-baseline.qza' % bodir
            baselines[model_baseline] = bstat
        out_paths = {
            'diff': '%s/differentials.tsv' % odir,
            'diff_qza': '%s/differentials.qza' % odir,
            'stat': '%s/differentials-stats.qza' % odir,
            'plot': '%s/differentials-biplot.qza' % odir,
            'tens': '%s/tensorboard.qzv' % bodir,
            'html': '%s/tensorboard.html' % bodir,
            'bdiff_qza': bdiff_qza,
            'bstat': bstat,
            'bplot': bplot
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
        new_meta_pd.to_csv(new_meta, index=False, sep='\t')
        return new_meta_pd.shape[0]

    def summarize_songbirds(self):
        q2s = []
        songbird = get_analysis_folder(
            self.config.i_datasets_folder, 'songbird')
        for root, dirs, files in os.walk(songbird):
            for fil in files:
                if fil == 'tensorboard.html':
                    path = root + '/' + fil
                    diff = '%s/differentials.tsv' % dirname(root)
                    root_split = root.split('%s/' % songbird)[-1].split('/')
                    d, pr, fr, sb, sr, ps, ml, be = root_split
                    with open(path) as f:
                        for line in f:
                            if 'Pseudo Q-squared' in line:
                                ls = line.split(
                                    'Pseudo Q-squared:</a></strong> ')
                                q2s.append([
                                    pr, d, fr, sb, ml, sr, ps, be, diff,
                                    float(ls[-1].split('<')[0])
                                ])
        if q2s:
            self.q2s_pd = pd.DataFrame(q2s, columns=[
                'pair', 'dataset', 'filter', 'subset', 'model',
                'songbird_filter', 'parameters', 'baseline', 'differentials',
                'Pseudo_Q_squared'])
            q2s_fp = '%s/songbird_q2.tsv' % songbird
            self.q2s_pd.to_csv(q2s_fp, index=False, sep='\t')
            print('[%s]\t\t==> Written: %s' % (datetime.datetime.now(), q2s_fp))

    def create_songbird_feature_metadata(self):
        if self.q2s_pd.shape[0]:
            q2_pd = self.q2s_pd.loc[(self.q2s_pd.pair == 'no_pair') &
                                    (self.q2s_pd.Pseudo_Q_squared > 0)]
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
                    sb_pd = pd.read_csv(diffs, index_col=0, sep='\t')
                    sb_pd.columns = ['%s %s: %s' % (
                        '__'.join([dat, pr, fr, sb, ml, st, ps, be])
                        , q2, x) for x in sb_pd.columns]
                    dataset_sbs.append(sb_pd)
                if len(dataset_sbs):
                    dataset_sbs_pd = pd.concat(dataset_sbs, axis=1, sort=False)
                    odir = get_analysis_folder(self.config.i_datasets_folder,
                                               'songbird/%s' % dat)
                    fpo_tsv = '%s/differentials_%s.tsv' % (odir, dat)
                    fpo_qza = '%s/differentials_%s.qza' % (odir, dat)
                    dataset_sbs_pd = dataset_sbs_pd.reset_index()
                    dataset_sbs_pd = dataset_sbs_pd.rename(
                        columns={
                            dataset_sbs_pd.columns.tolist()[0]: 'Feature ID'})
                    dataset_sbs_pd.to_csv(fpo_tsv, index=True, sep='\t')
                    run_import(fpo_tsv, fpo_qza, 'FeatureData[Differential]')

    def get_songbird_pd(self, songbird):
        self.songbird_pd = pd.DataFrame(
            songbird, columns=[
                'dataset', 'qza', 'meta', 'filter', 'params',
                'subset', 'differentials', 'baseline', 'html', 'pair'
            ])

    def check_metadata_models(self, meta, meta_pd, songbird_models):
        models = {}
        for model, formula_ in songbird_models.items():
            vars = set()
            drop = {}
            formula = formula_.strip('"').strip("'")
            # print()
            # print()
            # print('formula:', formula)
            if formula.startswith('C('):
                formula_split = formula.split('C(')[-1].rsplit(')', 1)
                # print('formula_split:', formula_split)
                formula_split_c = formula_split[0].split(',')[0].strip().strip()
                # print('formula_split_c:', formula_split_c)
                formula = 'C(%s)' % formula_split[0].replace(
                    formula_split_c, formula_split_c)
                # print('formula:', formula)
                formula += formula_split[1]
                # print('formula:', formula)
                if 'Diff' in formula:
                    levels = {formula_split_c: [
                        x.strip().strip('"').strip("'")
                        for x in formula.split(
                            "levels=['")[-1].split("']")[0].split(",")
                    ]}
                elif "Treatment('" in formula:
                    levels = {formula_split_c: [
                        formula.split("Treatment('")[-1].split("')")[0]
                    ]}
                    # print('levels:', levels)
                elif 'Treatment("' in formula:
                    levels = {formula_split_c: [
                        formula.split('Treatment("')[-1].split('")')[0]
                    ]}
                    # print('levels:', levels)
                vars.add(formula_split_c)
                vars.update(set([x for x in re.split(
                    '[+/:*]', formula_split[1]) if x]))
                # print("vars:", vars)
            else:
                formula_split = re.split('[+/:*]', formula)
                formula = formula
                vars.update(set([x for x in formula_split]))
                levels = {}

            common_with_md = set(meta_pd.columns.values) & vars
            if sorted(vars) != sorted(common_with_md):
                only_formula = sorted(vars ^ common_with_md)
                issue = 'Songbird formula term(s) missing in metadata:\n\t' \
                        '%s\n\t  [not used]: %s=%s' % (
                            ', '.join(only_formula), model, formula)
                self.models_issues.setdefault(issue, set()).add(meta)
                # print(issue)
                continue

            if levels:
                levels_set = sorted([x for x in meta_pd[
                    formula_split_c].unique() if str(x) != 'nan'])
                # print('levels_set:', levels_set)
                if 'Diff' in formula:
                    cur_levels = levels[formula_split_c]
                    common_levels = set(levels_set) & set(cur_levels)
                    only_meta = set(levels_set) ^ common_levels
                    only_model = set(cur_levels) ^ common_levels
                    if len(only_model):
                        issue = 'Songbird formula "Diff" factors(s) missing' \
                                ' in metadata "%s": %s' % (
                                    formula_split_c, list(only_model))
                        self.models_issues.setdefault(issue, set()).add(meta)
                        continue
                    if len(only_meta):
                        drop[formula_split_c] = list(only_meta)
                        issue = 'Songbird formula "Diff" factors(s) ' \
                                'incomplete for metadata "%s":\n' \
                                '\t -> skipping samples with %s' % (
                                    formula_split_c, list(only_meta))
                        self.models_issues.setdefault(issue, set()).add(meta)
                elif 'Treatment(' in formula:
                    levels = {formula_split_c: formula.split(
                        "Treatment('")[-1].split("')")[0]}
                    # print("levels (2):", levels)
                    if levels[formula_split_c] not in levels_set:
                        issue = 'Songbird formula "Treatment" factors(s)' \
                                ' missing in metadata "%s" [%s]' % (
                                    formula_split_c, levels)
                        self.models_issues.setdefault(issue, set()).add(meta)
                        continue

            models[model] = [formula, vars, drop]
        return models

    def show_models_issues(self, mess):
        if mess:
            for m in sorted(mess):
                print('[%s]' % datetime.datetime.now(), m)
        if self.models_issues:
            print('\n[%s] %s Issues with model (will not run) %s' % (
                datetime.datetime.now(), '#'*10, '#'*10))
            for model_issue, metas in self.models_issues.items():
                print('-', model_issue)
                for meta in metas:
                    print('\t', meta.replace(self.config.i_datasets_folder, ''))
            print('#'*60)

    def make_qurros(self) -> None:
        """Make qurro plots"""
        cmds = {}
        for r, row in self.songbird_pd.iterrows():
            dat = row['dataset']
            tax = self.project.datasets[dat].tax[-1]
            qurro_qzv = '%s_qurro.qzv' % splitext(row['differentials'])[0]
            if not isfile(qurro_qzv) and isfile(row['differentials']):
                cmd = 'qiime qurro differential-plot'
                cmd += ' --i-table %s' % row['qza']
                cmd += ' --i-ranks %s' % row['differentials']
                cmd += ' --m-sample-metadata-file %s' % row['meta']
                cmd += ' --m-feature-metadata-file %s' % tax
                cmd += ' --o-visualization %s' % qurro_qzv
                cmds.setdefault(dat, []).append(cmd)
        self.register_command('qurro', cmds)

    def songbird(self) -> None:
        """Main script for the creation of songbird jobs.
        It iterates over the rows of the table created
        upfront and over each combination of parameters
        and collect the output info for potential reuse
        in figure generation and post-analysis.

        Parameters
        ----------
        config : Class instance of AnalysesConfig
            Contains all the routine analyses config info.
        project
            Darasets.
        """
        mess = set()
        songbird = []
        dat_cmds, dat_fcmds, dat_bcmds = {}, {}, {}
        params_pd = self.get_params_combinations()
        for r, row in self.songbirds.iterrows():
            qza, pair, meta_fp = row['qza'], row['pair'], row['meta']
            dat, filt, subset = row['dataset'], row['filter'], row['subset']
            if dat not in self.songbird_models:
                continue
            dat_pair, pair_dir = self.get_dat_pair_dir(dat, pair)
            meta_pd = read_meta_pd(meta_fp)
            models = self.check_metadata_models(
                meta_fp, meta_pd, self.songbird_models[dat])
            row_params_pd = params_pd.copy()
            self.process_params_combinations(dat, meta_pd, row_params_pd, mess)
            for p, params in row_params_pd.iterrows():
                params_dir = self.get_params_dir(params)
                baselines, model_baselines = {}, {'1': '1'}
                for modx, model in enumerate(models.keys()):
                    formula, meta_vars, drop = models[model]
                    datdir, odir, new_qza, new_meta = self.get_main_dirs(
                        pair_dir, filt, subset, params_dir, model, self.config)
                    nsams = self.write_new_meta(
                        meta_pd, new_meta, meta_vars, drop, params)
                    if dat in self.models_baselines and model in \
                            self.models_baselines[dat]:
                        model_baselines = self.models_baselines[dat][model]
                    for model_baseline in model_baselines:
                        bformula = model_baselines[model_baseline]
                        bodir = get_analysis_folder(
                            self.config.i_datasets_folder,
                            'songbird/%s/b-%s' % (datdir, model_baseline))
                        out_paths = self.get_out_paths(
                            odir, bodir, model_baseline, baselines)
                        # convergence = self.check_stats_convergence(out_paths)
                        cmd, fcmd, bcmd = songbird_cmd(
                            qza, new_qza, new_meta, nsams, params, formula,
                            bformula, out_paths)
                        songbird.append([
                            dat, new_qza, meta_fp, filt, '%s_%s' % (
                                params_dir.replace('/', '__'), model),
                            subset, out_paths['diff'], model_baseline,
                            out_paths['html'], pair])
                        if cmd:
                            dat_cmds.setdefault(dat, []).append(cmd)
                        if fcmd:
                            dat_fcmds.setdefault(dat, []).append(fcmd)
                        if bcmd:
                            dat_bcmds.setdefault(dat, []).append(bcmd)
        if songbird:
            self.get_songbird_pd(songbird)
        self.show_models_issues(mess)
        self.register_command('songbird_filter', dat_fcmds)
        self.register_command('songbird_baselines', dat_bcmds)
        self.register_command('songbird', dat_cmds)
        self.summarize_songbirds()
        self.create_songbird_feature_metadata()

    @staticmethod
    def register_command(analysis, cmds):
        AnalysisPrep.analyses_commands[analysis] = cmds
