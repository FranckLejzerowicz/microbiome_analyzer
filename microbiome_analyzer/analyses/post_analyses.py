# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import pandas as pd
import pkg_resources
from scipy.stats import spearmanr
from os.path import dirname, isfile, splitext
from skbio.stats.ordination import OrdinationResults

from microbiome_analyzer.analysis import AnalysisPrep
from microbiome_analyzer._io import get_output
from microbiome_analyzer._taxonomy import get_split_taxonomy
from microbiome_analyzer._cmds import (
    get_biplot_commands, get_xmmvec_commands, get_paired_heatmaps_command)

RESOURCES = pkg_resources.resource_filename(
    "routine_qiime2_analyses", "resources/python_scripts")


class PostAnalysis(object):

    def __init__(self, config, project) -> None:
        self.config = config
        self.project = project
        self.cmds = {}
        self.mmvec_songbird_pd = pd.DataFrame()
        self.taxo_pds = {}
        self.metas = {}
        self.mmvec_res = {}
        self.mmvec_issues = set()
        self.xmmvecs = self.config.xmmvec
        self.highlights = self.config.mmvec_highlights

    def merge_mmvec_songbird(self, songbird_pd):
        rename_dic1 = dict(
            (x, '%s_omic1_songbird_common_fp' % x) for x in songbird_pd.columns)
        rename_dic1.update({'pair_omic_subset_filt': 'pair_omic_subset_filt1'})
        self.mmvec_songbird_pd = self.mmvec_songbird_pd.merge(
            songbird_pd.rename(columns=rename_dic1),
            on='pair_omic_subset_filt1', how='left')
        rename_dic2 = dict(
            (x, '%s_omic2_songbird_common_fp' % x) for x in songbird_pd.columns)
        rename_dic2.update({'pair_omic_subset_filt': 'pair_omic_subset_filt2'})
        self.mmvec_songbird_pd = self.mmvec_songbird_pd.merge(
            songbird_pd.rename(columns=rename_dic2),
            on='pair_omic_subset_filt2', how='left')

    def prep_mmvec(self, mmvec_pd):
        mmvec_pd = mmvec_pd.set_index(mmvec_pd.columns.tolist()[:-1]).unstack()
        mmvec_pd.columns = mmvec_pd.columns.droplevel()
        mmvec_pd.reset_index(inplace=True)
        mmvec_pd['filter1'] = mmvec_pd['pr1'] + '_' + mmvec_pd['ab1']
        mmvec_pd['filter2'] = mmvec_pd['pr2'] + '_' + mmvec_pd['ab2']
        mmvec_pd['omic_subset_filt1'] = mmvec_pd.apply(
            lambda x: '__'.join(x[['omic1', 'subset', 'filter1']]), axis=1)
        mmvec_pd['omic_subset_filt2'] = mmvec_pd.apply(
            lambda x: '__'.join(x[['omic1', 'subset', 'filter2']]), axis=1)
        mmvec_pd['pair_omic_subset_filt1'] = mmvec_pd.apply(
            lambda x: '__'.join(x[['pair', 'omic_subset_filt1']]), axis=1)
        mmvec_pd['pair_omic_subset_filt2'] = mmvec_pd.apply(
            lambda x: '__'.join(x[['pair', 'omic_subset_filt2']]), axis=1)
        self.mmvec_songbird_pd = mmvec_pd

    @staticmethod
    def prep_songbird(songbird_pd):
        songbird_pd['pair_omic_subset_filt'] = songbird_pd.apply(
            lambda x: '__'.join(x[['pair', 'dataset', 'subset', 'filter']]),
            axis=1)
        songbird_pd = songbird_pd[
            ['params', 'pair_omic_subset_filt', 'differentials']
        ].drop_duplicates().pivot(
            columns='params', index='pair_omic_subset_filt')
        songbird_pd.columns = songbird_pd.columns.droplevel()
        songbird_pd = songbird_pd.reset_index()
        return songbird_pd

    def get_tax_fp(self, omic):
        if isfile(self.project.datasets[omic].tax[-1]):
            omic_tax_fp = self.project.datasets[omic].tax[-1]
        else:
            omic_tax_fp = ''
        return omic_tax_fp

    def get_taxo_pds(self):
        for omicn in ['1', '2']:
            for omic in self.mmvec_songbird_pd['omic%s' % omicn].unique():
                omic_tax_fp = self.get_tax_fp(omic)
                if isfile(omic_tax_fp):
                    omic_tax_pd = pd.read_csv(
                        omic_tax_fp, header=0, sep='\t', dtype=str)
                    omic_tax_pd.rename(
                        columns={omic_tax_pd.columns[0]: 'Feature ID'},
                        inplace=True)
                else:
                    omic_tax_pd = pd.DataFrame()
                self.taxo_pds[omic] = omic_tax_pd

    def get_omics_songbirds_taxa(self):
        for omicn in ['1', '2']:
            pair_case_omics_filts = ['pair', 'subset', 'omic1',
                                     'filter1', 'omic2', 'filter2']
            all_omic_sb = [x for x in self.mmvec_songbird_pd.columns if
                           x.endswith('omic%s_songbird_common_fp' % omicn)]
            omicn_songbirds = self.mmvec_songbird_pd[
                (pair_case_omics_filts + all_omic_sb)].set_index(
                pair_case_omics_filts).T.to_dict()
            for (pair, subset, omic1, filt1, omic2,
                 filt2), sb_head_diff_fp in omicn_songbirds.items():
                if omicn == '1':
                    omic = omic1
                    omic_ = omic2
                    filt = filt1
                    filt_ = filt2
                else:
                    omic = omic2
                    omic_ = omic1
                    filt = filt2
                    filt_ = filt1
                feats_diff_cols = []
                cur_mmvec_folder = get_output(
                    self.config.i_datasets_folder,
                    'mmvec/metadata/%s/%s' % (pair, subset))
                omic_diff_list = []
                if len(sb_head_diff_fp):
                    for sb_head, diff_fp in sb_head_diff_fp.items():
                        model = sb_head.replace(
                            '_omic%s_songbird_common_fp' % omicn, '')
                        if str(diff_fp) != 'nan' and isfile(diff_fp):
                            diff_pd = pd.read_csv(diff_fp, header=0, sep='\t',
                                                  dtype=str)
                            index_header = diff_pd.columns[0]
                            if diff_pd[index_header][0] == '#q2:types':
                                diff_pd = diff_pd[1:]
                            diff_pd = diff_pd.rename(
                                columns={index_header: 'Feature ID'}).set_index(
                                'Feature ID')
                            diff_pd = diff_pd.drop(
                                columns=[x for x in diff_pd.columns if
                                         'Intercept' in x])
                            q2s = {}
                            diff_htmls = glob.glob(
                                '%s/*/tensorboard.html' % dirname(diff_fp))
                            if len(diff_htmls):
                                for diff_html in diff_htmls:
                                    baseline = diff_html.split('/')[-2]
                                    with open(diff_html) as f:
                                        for line in f:
                                            if 'Pseudo Q-squared' in line:
                                                q2 = line.split(
                                                    'Pseudo Q-squared:</a></strong> ')[
                                                    -1].split('<')[0]
                                                if float(q2) > 0.01:
                                                    q2s[baseline] = q2
                                                break
                            if q2s:
                                diff_cols = ['%s__%s__%s' % (
                                    model, x, '--'.join(
                                        ['%s-Q2=%s' % (b, q) if b else 'noQ2'
                                         for b, q in q2s.items()])
                                ) for x in diff_pd.columns]
                                diff_pd.columns = diff_cols
                                feats_diff_cols.extend(diff_cols)
                                omic_diff_list.append(diff_pd)
                if len(omic_diff_list):
                    omic_songbird_ranks = pd.concat(
                        omic_diff_list, axis=1, sort=False).reset_index()
                    omic_songbird_ranks.rename(
                        columns={omic_songbird_ranks.columns[0]: 'Feature ID'},
                        inplace=True
                    )
                else:
                    omic_common_fp = self.mmvec_songbird_pd.loc[
                        (self.mmvec_songbird_pd['pair'] == pair) &
                        (self.mmvec_songbird_pd['subset'] == subset) &
                        (self.mmvec_songbird_pd['omic1'] == omic1) &
                        (self.mmvec_songbird_pd['filter1'] == filt1) &
                        (self.mmvec_songbird_pd['omic2'] == omic2) &
                        (self.mmvec_songbird_pd['filter2'] == filt2),
                        'omic%s_common_fp' % omicn
                    ].tolist()[0]
                    omic_tax_list = []
                    if not isfile(omic_common_fp):
                        continue
                    with open(omic_common_fp) as f:
                        for ldx, line in enumerate(f):
                            if ldx:
                                omic_tax_list.append([line.split('\t')[0]])
                    omic_songbird_ranks = pd.DataFrame(omic_tax_list,
                                                       columns=['Feature ID'])
                if omic in self.taxo_pds:
                    omic_tax_pd = self.taxo_pds[omic]
                    if omic_tax_pd.shape[0]:
                        if 'Taxon' in omic_tax_pd.columns:
                            omic_split_taxa_pd = get_split_taxonomy(
                                omic_tax_pd.Taxon.tolist(), True)
                            omic_tax_pd = pd.concat(
                                [omic_tax_pd, omic_split_taxa_pd], axis=1,
                                sort=False)
                        omic_songbird_ranks = omic_songbird_ranks.merge(
                            omic_tax_pd, on='Feature ID',
                            how='left').drop_duplicates()
                meta_omic_fp = '%s/feature_metadata_%s_%s__%s_%s.tsv' % (
                    cur_mmvec_folder, omic, filt, omic_, filt_)
                drop_columns = [col for col in omic_songbird_ranks.columns if
                                omic_songbird_ranks[col].unique().size == 1]
                meta_omic_pd = omic_songbird_ranks.drop(columns=drop_columns)
                meta_omic_pd.to_csv(meta_omic_fp, index=False, sep='\t')
                meta_omic_pd.set_index('Feature ID', inplace=True)
                self.metas[(pair, subset, omic, filt, omic_, filt_)] = (
                    meta_omic_fp, meta_omic_pd, feats_diff_cols)

    @staticmethod
    def get_qzs(ordi_fp):
        qza = ordi_fp.replace('.txt', '.qza')
        qzv = ordi_fp.replace('.txt', '_emperor.qzv')
        return qza, qzv

    @staticmethod
    def edit_ordi_qzv(ordi, ordi_fp, highlight, regexes_list, meta, meta_pd):

        to_keep_feats = {}
        for regex in regexes_list:
            to_keep_feats[
                regex.lower()] = ordi.features.index.str.lower().str.contains(
                regex.lower())
        to_keep_feats_pd = pd.DataFrame(to_keep_feats)
        to_keep_feats = to_keep_feats_pd.any(axis=1)
        feats_subset_list = ordi.features.index[to_keep_feats].tolist()

        if feats_subset_list:
            ordi_edit = '%s_%s%s' % (
                splitext(ordi_fp)[0], highlight, splitext(ordi_fp)[1])
            ordi.features = ordi.features.loc[feats_subset_list, :].copy()
            ordi.write(ordi_edit)
            n_edit = ordi.features.shape[0]

            meta_edit = '%s_%s%s' % (
                splitext(meta)[0], highlight, splitext(meta)[1])
            meta_edit_pd = meta_pd.loc[feats_subset_list, :].copy()
            meta_edit_pd.to_csv(meta_edit, index=True, sep='\t')
        else:
            ordi_edit = ''
            meta_edit = ''
            n_edit = 0
        return n_edit, meta_edit, ordi_edit

    @staticmethod
    def get_pc_sb_correlations(pair, case, ordi, omic1, omic2, filt1, filt2,
                               diff_cols1, meta_pd1, diff_cols2, meta_pd2,
                               meta_fp, omic1_common_fp, omic2_common_fp,
                               ranks_fp):
        corrs = []
        max_r = 0
        for r in range(3):
            if ordi.features.shape[1] > r:
                max_r = r
                feats = ordi.features[r]
                if len(diff_cols1):
                    for model in diff_cols1:
                        x = meta_pd1.loc[
                            [x for x in meta_pd1.index if x in feats.index],
                            model
                        ].astype(float)
                        x = x[x.notnull()]
                        y = feats[x.index]
                        r2, p2 = spearmanr(x, y)
                        corrs.append([pair, case, omic1, filt1,
                                      'PC%s' % (r + 1), model, r2, p2,
                                      'spearman',
                                      meta_fp, omic1_common_fp, ranks_fp])
                sams = ordi.samples[r]
                if len(diff_cols2):
                    for model in diff_cols2:
                        x = meta_pd2.loc[
                            [x for x in meta_pd2.index if x in sams.index],
                            model
                        ].astype(float)
                        x = x[x.notnull()]
                        y = sams[x.index]
                        r2, p2 = spearmanr(x, y)
                        corrs.append([pair, case, omic2, filt2,
                                      'PC%s' % (r + 1), model, r2, p2,
                                      'spearman',
                                      meta_fp, omic2_common_fp, ranks_fp])
        corrs_pd = pd.DataFrame(corrs, columns=[
            'pair',
            'case',
            'omic',
            'filt',
            'mmvec_pc',
            'model',
            'correlation_coefficient',
            'pvalue',
            'correlation_method',
            'meta_fp',
            'features_fp',
            'ranks_fp'
        ])
        return corrs_pd, max_r

    def get_mmvec_res(self):
        mmvec_out_cols = [x for x in self.mmvec_songbird_pd.columns if
                          x.startswith('mmvec_out__')]
        # for ech row of the main table that also
        # contain the mmvec output folders
        for r, row in self.mmvec_songbird_pd.iterrows():
            pair = row['pair']
            subset = row['subset']
            omic1 = row['omic1']
            omic2 = row['omic2']
            filt1 = row['filter1']
            filt2 = row['filter2']
            omic1_common_fp = row['omic1_common_fp']
            if str(omic1_common_fp) == 'nan':
                continue
            omic2_common_fp = row['omic2_common_fp']
            n_common = row['n_common']
            meta_fp = row['meta_common_fp']
            # for each mmvec-parameters result
            for mmvec_out_col in mmvec_out_cols:
                # get the current parameters output result folder
                mmvec_out = row[mmvec_out_col]
                if str(mmvec_out) == 'nan':
                    continue
                # get the ordination file and the ranks file and skip
                # + warning if not performed
                mmvec_out_ranks = mmvec_out + '/model/ranks.tsv'
                mmvec_out_ordi = mmvec_out + '/model/ordination.txt'
                if not isfile(mmvec_out_ranks) or not isfile(mmvec_out_ordi):
                    issue = '\t\t[run mmvec first] %s (%s) %s (%s)' % (
                        omic1, filt1, omic2, filt2)
                    self.mmvec_issues.add(issue)
                    continue
                # collect the ranks + ordination
                # + songbirds for each pair of omics and parameters
                self.mmvec_res[
                    (pair, subset, omic1, omic2, filt1, filt2, n_common,
                     mmvec_out_col.replace('mmvec_out__', ''))
                ] = [mmvec_out_ranks, mmvec_out_ordi, meta_fp,
                     omic1_common_fp, omic2_common_fp]

    @staticmethod
    def get_order_omics(omic1, omic2, filt1, filt2, case, omics_pairs):
        omic_feature, omic_sample = ('feature', 'sample')
        omic_microbe, omic_metabolite = ('microbe', 'metabolite')
        omic_filt1 = '%s__%s__%s' % (omic1, case, filt1)
        omic_filt2 = '%s__%s__%s' % (omic2, case, filt2)
        if (omic_filt2, omic_filt1) in omics_pairs:
            omic_feature, omic_sample = ('sample', 'feature')
            omic_microbe, omic_metabolite = ('metabolite', 'microbe')
            omic1, omic2 = omic2, omic1
            filt1, filt2 = filt2, filt1
        return omic1, omic2, filt1, filt2, omic_feature, omic_sample, omic_microbe, omic_metabolite

    def get_pair_cmds(self, omics_pairs):
        crowdeds = [0, 1]
        pc_sb_correlations = []
        pre_paired_fp = '%s/mmvec_pre_paired-heatmaps.py' % RESOURCES
        for keys, values in self.mmvec_res.items():
            pair, case, omic1, omic2, filt1, filt2, sams, mmvec = keys
            ranks_fp, ordi_fp, meta_fp, omic1_common, omic2_common = values
            order_omics = self.get_order_omics(omic1, omic2, filt1, filt2,
                                               case, omics_pairs)
            omic1 = order_omics[0]
            omic2 = order_omics[1]
            filt1 = order_omics[2]
            filt2 = order_omics[3]
            omic_feature = order_omics[4]
            omic_sample = order_omics[5]
            omic_microbe = order_omics[6]
            omic_metabolite = order_omics[7]

            # get differentials
            meta1, meta_pd1, diff_cols1 = self.metas[(pair, case, omic1,
                                                      filt1, omic2, filt2)]
            meta2, meta_pd2, diff_cols2 = self.metas[(pair, case, omic2,
                                                      filt2, omic1, filt1)]
            # features are biplot, samples are dots
            ordi = OrdinationResults.read(ordi_fp)
            cur_pc_sb_correlations, max_r = self.get_pc_sb_correlations(
                pair, case, ordi, omic1, omic2, filt1, filt2, diff_cols1,
                meta_pd1, diff_cols2, meta_pd2, meta_fp, omic1_common,
                omic2_common, ranks_fp)
            pc_sb_correlations.append(cur_pc_sb_correlations)

            cmd = ''
            if pair in self.highlights:
                pair_highlights = self.highlights[pair]
                for highlight, regexes_list in pair_highlights.items():
                    n_edit, meta_edit, ordi_edit_fp = self.edit_ordi_qzv(
                        ordi, ordi_fp, highlight, regexes_list, meta1, meta_pd1)
                    if n_edit:
                        qza, qzv = self.get_qzs(ordi_edit_fp)
                        cmd += get_biplot_commands(
                            ordi_edit_fp, qza, qzv, omic_feature, omic_sample,
                            meta_edit, meta2, n_edit, max_r)
            ordi_edit_fp = ordi_fp
            qza, qzv = self.get_qzs(ordi_edit_fp)
            for crowded in crowdeds:
                if crowded:
                    n_ordi_feats = ordi.features.shape[0]
                    qzv = qzv.replace('.qzv', '_crowded.qzv')
                else:
                    n_ordi_feats = 15
                    # heat_qza, heat_qzv = get_heatmap_qzs(ranks_fp)
                    # cmd += get_heatmap_commands(
                    #     ranks_fp, heat_qza, heat_qzv, meta1,
                    #     meta2, meta_pd1, meta_pd2)
                cmd += get_biplot_commands(
                    ordi_edit_fp, qza, qzv, omic_feature, omic_sample,
                    meta1, meta2, n_ordi_feats, max_r)
            cmd += get_xmmvec_commands(
                ordi_edit_fp, omic1, omic2, meta1, meta2, self.xmmvecs, pair)

            topn = 5
            features_names = []
            if features_names:
                heat = '%s_paired_heatmaps_custom.qzv' % splitext(ranks_fp)[0]
            else:
                heat = '%s_paired_heatmaps_top%s.qzv' % (splitext(ranks_fp)[0],
                                                         topn)
            cmd += get_paired_heatmaps_command(
                ranks_fp, omic1_common, omic2_common, meta1, features_names,
                topn, heat, pre_paired_fp)
            self.cmds.setdefault(pair, []).append(cmd)
        return pc_sb_correlations

    def show_mmvec_issues(self):
        if self.mmvec_issues:
            for mmvec_issue in self.mmvec_issues:
                print(mmvec_issue)

    def mmbird(self, paired_datasets, differentials):
        if not paired_datasets.mmvec_pd.shape[0]:
            print('No mmvec output detected...')
            return None
        self.prep_mmvec(paired_datasets.mmvec_pd)
        if differentials.songbird_pd.shape[0]:
            songbird_pd = self.prep_songbird(differentials.songbird_pd)
            self.merge_mmvec_songbird(songbird_pd)

        self.get_taxo_pds()
        self.get_omics_songbirds_taxa()
        self.get_mmvec_res()
        self.show_mmvec_issues()

        omics_pairs = [tuple(x) for x in self.mmvec_songbird_pd[
            ['omic_subset_filt1', 'omic_subset_filt2']].values.tolist()]
        pc_sb_correlations = self.get_pair_cmds(omics_pairs)

        if len(pc_sb_correlations):
            out_folder = get_output(self.config.i_datasets_folder,
                                             'mmbird')
            out_correlations = '%s/pc_vs_songbird_correlations.tsv' % out_folder
            pc_sb_correlations_pd = pd.concat(pc_sb_correlations)
            if pc_sb_correlations_pd.shape[0]:
                pc_sb_correlations_pd.to_csv(
                    out_correlations, index=False, sep='\t')
                print('\t\t==> Written:', out_correlations)
            else:
                print('\t\t==> No good songbird model to '
                      'make correlations with mmvec PCs...')
        self.register_command('mmbird')

    def register_command(self, analysis):
        AnalysisPrep.analyses_commands[analysis] = self.cmds
