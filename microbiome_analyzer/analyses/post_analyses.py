# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import pandas as pd
import pkg_resources
from scipy.stats import spearmanr
from os.path import dirname, isdir, isfile, splitext
from skbio.stats.ordination import OrdinationResults

from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer._scratch import io_update, to_do, rep
from microbiome_analyzer.analyses.taxonomy import get_split_taxonomy
from microbiome_analyzer.core.commands import (
    get_biplot_commands, get_xmmvec_commands)

RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources/python_scripts")


class PostAnalysis(object):

    def __init__(self, config, project) -> None:
        self.config = config
        self.project = project
        self.dir = project.dir
        self.dirs = project.dirs
        self.out = ''
        self.analysis = None
        self.ios = {}
        self.cmds = {}
        self.mmvec_songbird_pd = pd.DataFrame()
        self.taxo_pds = {}
        self.metas = {}
        self.mmvec_res = {}
        self.mmvec_issues = set()
        self.xmmvecs = self.config.xmmvec
        self.highlights = self.config.mmvec_highlights

    def get_output(self, dat: str = '') -> str:
        out = '%s/%s' % (self.dir, self.analysis)
        if dat:
            out += '/%s' % dat
        self.out = out
        if not isdir(rep(self.out)):
            os.makedirs(rep(self.out))
        return out

    def get_path(self, path_):
        path = path_
        params = self.config.run_params[self.analysis]
        if not self.config.jobs or not params['scratch']:
            path = rep(path_)
        return path

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
        if isfile(rep(self.project.datasets[omic].tax[-1])):
            omic_tax_fp = self.project.datasets[omic].tax[-1]
        else:
            omic_tax_fp = ''
        return omic_tax_fp

    def get_taxo_pds(self):
        for omicn in ['1', '2']:
            for omic in self.mmvec_songbird_pd['omic%s' % omicn].unique():
                omic_tax_fp = self.get_tax_fp(omic)
                if isfile(rep(omic_tax_fp)):
                    omic_tax_pd = pd.read_csv(
                        rep(omic_tax_fp), header=0, sep='\t', dtype=str)
                    omic_tax_pd.rename(
                        columns={omic_tax_pd.columns[0]: 'Feature ID'},
                        inplace=True)
                else:
                    omic_tax_pd = pd.DataFrame()
                self.taxo_pds[omic] = omic_tax_pd

    def get_omics_songbirds_taxa(self):
        self.analysis = 'mmvec'
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
                omic_diff_list = []
                feats_diff_cols = []
                self.get_output('metadata/%s/%s' % (pair, subset))
                if len(sb_head_diff_fp):
                    for sb_head, diff_fp in sb_head_diff_fp.items():
                        if str(diff_fp) == 'nan':
                            continue
                        model = sb_head.replace(
                            '_omic%s_songbird_common_fp' % omicn, '')
                        if str(rep(diff_fp)) != 'nan' and not to_do(diff_fp):
                            diff_pd = pd.read_csv(rep(diff_fp), header=0,
                                                  sep='\t', dtype=str)
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
                                '%s/*/tensorboard.html' % dirname(rep(diff_fp)))
                            if len(diff_htmls):
                                for diff_html in diff_htmls:
                                    baseline = diff_html.split('/')[-2]
                                    with open(rep(diff_html)) as f:
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
                                omic_diff_list.append(diff_pd)
                                feats_diff_cols.extend(diff_cols)
                if len(omic_diff_list):
                    omic_songbird_ranks = pd.concat(
                        omic_diff_list, axis=1, sort=False).reset_index()
                    omic_songbird_ranks.rename(
                        columns={omic_songbird_ranks.columns[0]: 'Feature ID'},
                        inplace=True)
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
                    if to_do(omic_common_fp):
                        continue
                    with open(rep(omic_common_fp)) as f:
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
                    self.out, omic, filt, omic_, filt_)
                drop_columns = [col for col in omic_songbird_ranks.columns if
                                omic_songbird_ranks[col].unique().size == 1]
                meta_omic_pd = omic_songbird_ranks.drop(columns=drop_columns)
                meta_omic_pd.to_csv(rep(meta_omic_fp), index=False, sep='\t')
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
            ordi.write(rep(ordi_edit))
            n_edit = ordi.features.shape[0]

            meta_edit = '%s_%s%s' % (
                splitext(meta)[0], highlight, splitext(meta)[1])
            meta_edit_pd = meta_pd.loc[feats_subset_list, :].copy()
            meta_edit_pd.to_csv(rep(meta_edit), index=True, sep='\t')
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
            pair, subset = row['pair'], row['subset']
            omic1, omic2 = row['omic1'], row['omic2']
            filt1, filt2 = row['filter1'], row['filter2']
            common1_fp = row['omic1_common_fp']
            if str(common1_fp) == 'nan':
                continue
            common2_fp = row['omic2_common_fp']
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
                if to_do(mmvec_out_ranks) or to_do(mmvec_out_ordi):
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
                     common1_fp, common2_fp]

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
        order_omics = (omic1, omic2, filt1, filt2, omic_feature,
                       omic_sample, omic_microbe, omic_metabolite)
        return order_omics

    def get_paired_heatmaps_command(
            self,
            pair,
            ranks_fp: str,
            omic1_common_fp: str,
            omic2_common_fp: str,
            taxonomy_tsv: str,
            features_names: list,
            topn: int,
            paired_heatmap_qzv: str,
            pre_paired_fp: str
    ):

        cmd = ''
        # if not isfile(paired_heatmap_qzv):
        if 1:
            omic1_tmp = '%s_tmp.tsv' % splitext(omic1_common_fp)[0]
            omic2_tmp = '%s_tmp.tsv' % splitext(omic2_common_fp)[0]
            omic1_qza_tmp = '%s_tmp.qza' % splitext(omic1_common_fp)[0]
            omic2_qza_tmp = '%s_tmp.qza' % splitext(omic2_common_fp)[0]
            taxonomy_tsv_tmp = '%s_tmp.tsv' % splitext(taxonomy_tsv)[0]
            ranks_fp_tmp = '%s_tmp.tsv' % splitext(ranks_fp)[0]
            ranks_qza_tmp = '%s_tmp.qza' % splitext(ranks_fp)[0]

            py = '%s.py' % splitext(paired_heatmap_qzv)[0]
            with open(rep(py), 'w') as o, open(rep(pre_paired_fp)) as f:
                for line in f:
                    if "'OMIC1_COMMON_FP_TMP'" in line:
                        o.write(line.replace('OMIC1_COMMON_FP_TMP',
                                             self.get_path(omic1_tmp)))
                    elif "'OMIC2_COMMON_FP_TMP'" in line:
                        o.write(line.replace('OMIC2_COMMON_FP_TMP',
                                             self.get_path(omic2_tmp)))
                    elif "'OMIC1_COMMON_QZA_TMP'" in line:
                        o.write(line.replace('OMIC1_COMMON_QZA_TMP',
                                             self.get_path(omic1_qza_tmp)))
                    elif "'OMIC2_COMMON_QZA_TMP'" in line:
                        o.write(line.replace('OMIC2_COMMON_QZA_TMP',
                                             self.get_path(omic2_qza_tmp)))
                    elif "'OMIC1_COMMON_FP'" in line:
                        o.write(line.replace('OMIC1_COMMON_FP',
                                             self.get_path(omic1_common_fp)))
                    elif "'OMIC2_COMMON_FP'" in line:
                        o.write(line.replace('OMIC2_COMMON_FP',
                                             self.get_path(omic2_common_fp)))
                    elif "'TAXONOMY_TSV_TMP'" in line:
                        o.write(line.replace('TAXONOMY_TSV_TMP',
                                             self.get_path(taxonomy_tsv_tmp)))
                    elif "'TAXONOMY_TSV'" in line:
                        o.write(line.replace('TAXONOMY_TSV',
                                             self.get_path(taxonomy_tsv)))
                    elif "'RANKS_FP_TMP'" in line:
                        o.write(line.replace('RANKS_FP_TMP',
                                             self.get_path(ranks_fp_tmp)))
                    elif "'RANKS_QZA_TMP'" in line:
                        o.write(line.replace('RANKS_QZA_TMP',
                                             self.get_path(ranks_qza_tmp)))
                    elif "'RANKS_FP'" in line:
                        o.write(line.replace('RANKS_FP',
                                             self.get_path(ranks_fp)))
                    else:
                        o.write(line)

            cmd += '\npython3 %s\n' % py

            cmd += '\nqiime mmvec paired-heatmap'
            cmd += ' --i-ranks %s' % ranks_qza_tmp
            cmd += ' --i-microbes-table %s' % omic1_qza_tmp
            cmd += ' --i-metabolites-table %s' % omic2_qza_tmp
            cmd += ' --m-microbe-metadata-file %s' % taxonomy_tsv_tmp
            cmd += ' --m-microbe-metadata-column Taxon'
            if features_names:
                cmd += ' --p-top-k-microbes 0'
                for features_name in features_names:
                    cmd += ' --p-features %s' % features_name
            else:
                cmd += ' --p-top-k-microbes %s' % topn
            cmd += ' --p-normalize rel_row'
            cmd += ' --p-top-k-metabolites 100'
            cmd += ' --p-level 6'
            cmd += ' --o-visualization %s\n' % paired_heatmap_qzv

            cmd += '\nrm %s %s %s %s\n' % (
                ranks_qza_tmp, omic1_qza_tmp, omic2_qza_tmp, taxonomy_tsv_tmp)
            io_update(self, i_f=py, o_f=paired_heatmap_qzv, key=pair)

        return cmd

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
            ordi = OrdinationResults.read(rep(ordi_fp))
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
                            self, pair, ordi_edit_fp, qza, qzv, omic_feature,
                            omic_sample, meta_edit, meta2, n_edit, max_r)
            ordi_edit_fp = ordi_fp
            qza, qzv = self.get_qzs(ordi_edit_fp)
            for crowded in crowdeds:
                if crowded:
                    n_ordi_feats = ordi.features.shape[0]
                    qzv = qzv.replace('.qzv', '_crowded.qzv')
                else:
                    n_ordi_feats = 15
                cmd += get_biplot_commands(
                    self, pair, ordi_edit_fp, qza, qzv, omic_feature,
                    omic_sample, meta1, meta2, n_ordi_feats, max_r)
            cmd += get_xmmvec_commands(
                self, pair, ordi_edit_fp, omic1, omic2, meta1, meta2)

            topn = 5
            features_names = []
            if features_names:
                heat = '%s_paired_heatmaps_custom.qzv' % splitext(ranks_fp)[0]
            else:
                heat = '%s_paired_heatmaps_top%s.qzv' % (splitext(ranks_fp)[0],
                                                         topn)
            cmd += self.get_paired_heatmaps_command(
                pair, ranks_fp, omic1_common, omic2_common, meta1,
                features_names, topn, heat, pre_paired_fp)
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

        self.analysis = 'mmbird'
        if len(pc_sb_correlations):
            self.get_output()
            out_correlations = '%s/pc_vs_songbird_correlations.tsv' % self.out
            pc_sb_correlations_pd = pd.concat(pc_sb_correlations)
            if pc_sb_correlations_pd.shape[0]:
                pc_sb_correlations_pd.to_csv(
                    rep(out_correlations), index=False, sep='\t')
                print('\t\t==> Written:', out_correlations)
            else:
                print('\t\t==> No good songbird model to '
                      'make correlations with mmvec PCs...')

        self.register_io_command()

    def register_io_command(self):
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        AnalysisPrep.analyses_ios[self.analysis] = dict(self.ios)
        self.ios = {}
        self.cmds = {}
