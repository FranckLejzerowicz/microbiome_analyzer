# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import glob
import pandas as pd
from os.path import isfile, splitext

from scipy.stats import spearmanr, pearsonr
from skbio.stats.ordination import OrdinationResults

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder, get_highlights
from routine_qiime2_analyses._routine_q2_cmds import get_split_taxonomy


def get_mmvec_outputs(mmvec_outputs: list):
    mmvec_outputs_pd = pd.DataFrame(
        mmvec_outputs,
        columns=[
             'pair',
             'omic1',
             'omic2',
             'filt1',
             'filt2',
             'n_common',
             'meta_common_fp',
             'omic1_common_fp',
             'omic2_common_fp',
             'omic1_common_qza',
             'omic2_common_qza',
             'mmvec_parameters',
             'mmvec_out'
         ])
    mmvec_outputs_pd = mmvec_outputs_pd.set_index(
        mmvec_outputs_pd.columns.tolist()[:-1]).unstack()
    mmvec_outputs_pd.columns = mmvec_outputs_pd.columns.droplevel()
    mmvec_outputs_pd.reset_index(inplace=True)
    mmvec_outputs_pd['omic_filt1'] = mmvec_outputs_pd['omic1'] + '__' + mmvec_outputs_pd['filt1']
    mmvec_outputs_pd['omic_filt2'] = mmvec_outputs_pd['omic2'] + '__' + mmvec_outputs_pd['filt2']
    mmvec_outputs_pd['pair_omic_filt1'] = mmvec_outputs_pd['pair'] + '__' + mmvec_outputs_pd['omic_filt1']
    mmvec_outputs_pd['pair_omic_filt2'] = mmvec_outputs_pd['pair'] + '__' + mmvec_outputs_pd['omic_filt2']
    return mmvec_outputs_pd


def get_songbird_outputs(songbird_outputs: list):
    songbird_outputs_pd = pd.DataFrame(
        songbird_outputs,
        columns=[
            'songbird_dat',
            'songbird_filt',
            'songbird_parameters',
            'songbird_case',
            'songbird_fp',
            'pair'
        ])
    songbird_outputs_pd['pair_omic_filt'] = songbird_outputs_pd['pair'] + '__' + songbird_outputs_pd[
        'songbird_dat'] + '__' + songbird_outputs_pd['songbird_filt']
    songbird_outputs_pd['case_params'] = songbird_outputs_pd['songbird_case'] + '__' + songbird_outputs_pd[
        'songbird_parameters']
    songbird_outputs_pd = songbird_outputs_pd.drop(
        columns=['pair', 'songbird_dat', 'songbird_filt', 'songbird_case', 'songbird_parameters'])
    songbird_outputs_pd = songbird_outputs_pd.pivot(columns='case_params', index='pair_omic_filt')
    songbird_outputs_pd.columns = songbird_outputs_pd.columns.droplevel()
    songbird_outputs_pd = songbird_outputs_pd.reset_index()
    return songbird_outputs_pd


def merge_mmvec_songbird_outputs(mmvec_outputs_pd, songbird_outputs_pd):

    rename_dic1 = dict((x, '%s_omic1_songbird_common_fp' % x) for x in songbird_outputs_pd.columns)
    rename_dic1.update({'pair_omic_filt': 'pair_omic_filt1'})
    mmvec_songbird_pd = mmvec_outputs_pd.merge(
        songbird_outputs_pd.rename(columns=rename_dic1),
        on='pair_omic_filt1',
        how='left'
    )
    rename_dic2 = dict((x, '%s_omic2_songbird_common_fp' % x) for x in songbird_outputs_pd.columns)
    rename_dic2.update({'pair_omic_filt': 'pair_omic_filt2'})
    mmvec_songbird_pd = mmvec_songbird_pd.merge(
        songbird_outputs_pd.rename(columns=rename_dic2),
        on='pair_omic_filt2',
        how='left'
    )
    return mmvec_songbird_pd


def get_mmvec_res(mmvec_songbird_pd):
    mmvec_out_cols = [x for x in mmvec_songbird_pd.columns if x.startswith('mmvec_out__')]
    mmvec_res = {}
    # for ech row of the main table that also contain the mmvec output folders
    for r, row in mmvec_songbird_pd.iterrows():
        pair = row['pair']
        omic1 = row['omic1']
        omic2 = row['omic2']
        filt1 = row['filt1']
        filt2 = row['filt2']
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

            # get the ordination file and the ranks file and skip + warning if not performed
            mmvec_out_ranks = mmvec_out + '/ranks.tsv'
            mmvec_out_ordi = mmvec_out + '/ordination.txt'
            if not isfile(mmvec_out_ranks) or not isfile(mmvec_out_ordi):
                print('[ - mmvec] %s (%s) %s (%s)' % (omic1, filt1, omic2, filt2))
                continue

            # collect the ranks + ordination + songbirds for each pair of omics and parameters
            mmvec_res[
                (
                    pair, omic1, omic2,
                    filt1, filt2, n_common,
                    mmvec_out_col.replace('mmvec_out__', '')
                )
            ] = [
                mmvec_out_ranks,
                mmvec_out_ordi,
                meta_fp,
                omic1_common_fp,
                omic2_common_fp
            ]
    return mmvec_res


def get_qzs(ordi_fp):
    qza = ordi_fp.replace('.txt', '.qza')
    qzv = ordi_fp.replace('.txt', '_emperor.qzv')
    return qza, qzv


def get_heatmap_qzs(ranks_fp):
    qza = ranks_fp.replace('.tsv', '.qza')
    qzv = ranks_fp.replace('.tsv', '_heatmap.qzv')
    return qza, qzv


def get_order_omics(
        omic1, omic2,
        filt1, filt2,
        omics_pairs
):
    omic_feature, omic_sample = ('feature', 'sample')
    omic_microbe, omic_metabolite = ('microbe', 'metabolite')
    omic_filt2 = '%s__%s' % (omic1, filt1)
    omic_filt1 = '%s__%s' % (omic2, filt2)
    if (omic_filt2, omic_filt1) not in omics_pairs:
        omic_feature, omic_sample = ('sample', 'feature')
        omic_microbe, omic_metabolite = ('metabolite', 'microbe')
        omic1, omic2 = omic2, omic1
        filt1, filt2 = filt2, filt1
    return omic1, omic2, filt1, filt2, omic_feature, omic_sample, omic_microbe, omic_metabolite


# def get_tax_extended_fps(
#         omic_filt,
#         omic_common_fp,
#         omic_tax_pd,
#         all_omic_songbird_ranks,
#         ordi_fp
# ):
#     if isfile(omic_tax_fp):
#         omic_tax_pd = pd.read_csv(omic_tax_fp, header=0, sep='\t', dtype=str)
#         if 'Taxon' in omic_tax_pd.columns:
#             omic_split_taxo = get_split_taxonomy(omic_tax_pd.Taxon.tolist())
#             omic_tax_pd = omic_tax_pd.merge(omic_split_taxo, on='Taxon', how='left').drop_duplicates()
#     else:
#         omic_tax_list = []
#         with open(omic_common_fp) as f:
#             for ldx, line in enumerate(f):
#                 if ldx:
#                     omic_tax_list.append([line.split('\t')[0]])
#         omic_tax_pd = pd.DataFrame(omic_tax_list, columns=['Feature ID'])
#
#     if all_omic_songbird_ranks.shape[0]:
#         print('all_omic_songbird_ranks.shape[0] > 0 !!!!')
#         omic_tax_pd = omic_tax_pd.merge(
#             all_omic_songbird_ranks,
#             on='Feature ID',
#             how='left'
#         ).drop_duplicates()
#     metatax_omic_fp = '%s_meta-%s.tsv' % (splitext(ordi_fp)[0], omic_filt)
#     omic_tax_pd.to_csv(metatax_omic_fp, index=False, sep='\t')
#
#     return metatax_omic_fp

def get_biplot_commands(
        ordi_edit_fp, qza, qzv, omic_feature, omic_sample,
        meta1_fp, meta2_fp, n_edit
):
    cmd = '\n'
    if not isfile(qza):
        cmd += '\nqiime tools import'
        cmd += ' --input-path %s' % ordi_edit_fp
        cmd += ' --output-path %s' % qza
        cmd += ' --type "PCoAResults %s Properties([\'biplot\'])"\nsleep 3' % '%'
    cmd += '\nqiime emperor biplot'
    cmd += ' --i-biplot %s' % qza
    cmd += ' --m-%s-metadata-file %s' % (omic_feature, meta1_fp)
    cmd += ' --m-%s-metadata-file %s' % (omic_sample, meta2_fp)
    cmd += ' --p-number-of-features %s' % n_edit
    cmd += ' --o-visualization %s\n' % qzv
    return cmd


def get_heatmap_commands(
        ranks_fp, qza, qzv, meta1, meta2,
        meta_pd1, meta_pd2):

    cmd = '\n'
    if not isfile(qza):
        cmd += '\nqiime tools import'
        cmd += ' --input-path %s' % ranks_fp
        cmd += ' --output-path %s' % qza
        cmd += ' --type FeatureData[Conditional]'
        cmd += ' --input-format ConditionalFormat\n'

    cmd += '\nqiime mmvec heatmap'
    cmd += ' --i-ranks %s' % qza

    for level in 'CB':
        taxolevel_microbe = 'Taxolevel_%s' % level
        if taxolevel_microbe in meta_pd1.columns:
            break
    else:
        taxolevel_microbe = ''

    for level in 'DCB':
        taxolevel_metabolite = 'Taxolevel_%s' % level
        if taxolevel_metabolite in meta_pd2.columns:
            break
    else:
        taxolevel_metabolite = ''

    if taxolevel_microbe:
        cmd += ' --m-microbe-metadata-file %s' % meta1
        cmd += ' --m-microbe-metadata-column %s' % taxolevel_microbe
    if taxolevel_metabolite:
        cmd += ' --m-metabolite-metadata-file %s' % meta2
        cmd += ' --m-metabolite-metadata-column %s' % taxolevel_metabolite

    cmd += ' --o-visualization %s\n' % qzv
    return cmd


def edit_ordi_qzv(ordi, ordi_fp, highlight, regexes_list, meta, meta_pd):

    to_keep_feats = {}
    for regex in regexes_list:
        to_keep_feats[regex.lower()] = ordi.features.index.str.lower().str.contains(regex.lower())
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


def get_tax_fp(i_datasets_folder: str, omic: str, input_to_filtered: dict) -> str:

    tax_dir = get_analysis_folder(i_datasets_folder, 'taxonomy')

    omic_taxs = [x for x, y in input_to_filtered.items() if y == omic]
    if len(omic_taxs):
        omic_tax_ = omic_taxs[0]
        if omic_tax_.endswith('__raref'):
            omic_tax = '__raref'.join(omic_tax_.split('__raref')[:-1])
        else:
            omic_tax = omic_tax_
    else:
        print('\nNo taxonomy file for "%s"' % omic)
        return ''

    omic_tax_fps = glob.glob('%s/%s/tax_%s*.tsv' % (tax_dir, omic_tax, omic_tax))
    if len(omic_tax_fps):
        omic_tax_fp = omic_tax_fps[0]
    else:
        omic_tax_fp = ''
    return omic_tax_fp


def get_pair_cmds(mmvec_res: dict, omics_pairs_metas: dict,
                  omics_pairs: list, force: bool, highlights: dict):
    crowdeds = [0, 1]
    pc_sb_correlations = []
    # mmvec_tab = []
    pair_cmds = {}
    for keys, values in mmvec_res.items():

        pair, omic1, omic2, filt1, filt2, sams, mmvec = keys
        ranks_fp, ordi_fp, meta_fp, omic1_common_fp, omic2_common_fp = values
        if ranks_fp != '/projects/nutrition/analysis/tree_based/qiime/mmvec/paired/foods_ssu/vioscreen_foods_consumed_grams_per_day_1800s_noLiquids_0_0__16S_100nt_745s_0_0/b-50_l-1e-4_e-200_p-05_f-0_d-3_t-None_n-200_gpu-F/ranks.tsv':
            continue
        print('-'*30)
        print('ranks_fp, ordi_fp, meta_fp, omic1_common_fp, omic2_common_fp')
        print(ranks_fp)
        print(ordi_fp)
        print(meta_fp)
        print(omic1_common_fp)
        print(omic2_common_fp)
        print('-'*50)

        order_omics = get_order_omics(omic1, omic2, filt1, filt2, omics_pairs)
        omic1 = order_omics[0]
        omic2 = order_omics[1]
        filt1 = order_omics[2]
        filt2 = order_omics[3]
        omic_feature = order_omics[4]
        omic_sample = order_omics[5]
        omic_microbe = order_omics[6]
        omic_metabolite = order_omics[7]

        # get differentials
        meta1, meta_pd1, diff_cols1 = omics_pairs_metas[(pair, omic1, filt1)]
        meta2, meta_pd2, diff_cols2 = omics_pairs_metas[(pair, omic2, filt2)]
        print('-'*50)
        print('meta1')
        print(meta1)
        print('meta_pd1')
        print(meta_pd1.iloc[:3,:3])
        print('meta2')
        print(meta2)
        print('meta_pd2')
        print(meta_pd2.iloc[:3,:3])
        print('-'*50)

        # features are biplot, samples are dots
        ordi = OrdinationResults.read(ordi_fp)

        cur_pc_sb_correlations = get_pc_sb_correlations(
            pair, ordi, omic1, omic2, filt1, filt2,
            diff_cols1, meta_pd1, diff_cols2, meta_pd2,
            meta_fp,  omic1_common_fp, omic2_common_fp, ranks_fp)
        print('-'*30)
        print("cur_pc_sb_correlations")
        print(cur_pc_sb_correlations)
        print('-'*30)
        pc_sb_correlations.append(cur_pc_sb_correlations)

        cmd = ''
        if pair in highlights:
            pair_highlights = highlights[pair]
            for highlight, regexes_list in pair_highlights.items():
                n_edit, meta_edit, ordi_edit_fp = edit_ordi_qzv(
                    ordi, ordi_fp, highlight, regexes_list, meta1, meta_pd1)
                if n_edit:
                    qza, qzv = get_qzs(ordi_edit_fp)
                    cmd += get_biplot_commands(
                        ordi_edit_fp, qza, qzv,
                        omic_feature, omic_sample,
                        meta_edit, meta2, n_edit)
        ordi_edit_fp = ordi_fp
        qza, qzv = get_qzs(ordi_edit_fp)
        for crowded in crowdeds:
            if crowded:
                n_ordi_feats = ordi.features.shape[0]
                qzv = qzv.replace('.qzv', '_crowded.qzv')
            else:
                n_ordi_feats = 25
                heat_qza, heat_qzv = get_heatmap_qzs(ranks_fp)
                cmd += get_heatmap_commands(
                    ranks_fp, heat_qza, heat_qzv, meta1,
                    meta2, meta_pd1, meta_pd2)
            cmd += get_biplot_commands(
                ordi_edit_fp, qza, qzv,
                omic_feature, omic_sample,
                meta1, meta2, n_ordi_feats)

        if cmd:
            pair_cmds.setdefault(pair, []).append(cmd)

    pc_sb_correlations_pd = pd.concat(pc_sb_correlations)
    return pair_cmds, pc_sb_correlations_pd


def get_omics_songbirds_taxa(i_datasets_folder, mmvec_songbird_pd, taxo_pds):

    omics_pairs_metas = {}
    for omicn in ['1', '2']:
        pair_omic_filt = ['pair', 'omic%s' % omicn, 'filt%s' % omicn]
        all_omic_sb = [x for x in mmvec_songbird_pd.columns if x.endswith('omic%s_songbird_common_fp' % omicn)]
        omicn_songbirds = mmvec_songbird_pd[
            (pair_omic_filt + all_omic_sb)
        ].set_index(pair_omic_filt).T.to_dict()
        for (pair, omic, filt), sb_head_diff_fp in omicn_songbirds.items():
            # print()
            # print()
            # print('-' *50)
            # print(pair, omic, filt)
            # print(sb_head_diff_fp)
            feats_diff_cols = []
            cur_mmvec_folder = get_analysis_folder(i_datasets_folder, 'mmvec/metadata/%s' % pair)
            omic_diff_list = []
            if len(sb_head_diff_fp):
                for sb_head, diff_fp in sb_head_diff_fp.items():
                    model = sb_head.replace('_omic%s_songbird_common_fp' % omicn, '')
                    if str(diff_fp) != 'nan' and isfile(diff_fp):
                        diff_html = '%s-tensorboard.html' % splitext(diff_fp)[0]
                        if isfile(diff_html):
                            with open(diff_html) as f:
                                for line in f:
                                    if 'Pseudo Q-squared' in line:
                                        q2 = line.split('Pseudo Q-squared:</a></strong> ')[-1].split('<')[0]
                                        break
                        else:
                            q2 = 'Not_found'
                        diff_pd = pd.read_csv(diff_fp, header=0, sep='\t', dtype=str)
                        index_header = diff_pd.columns[0]
                        if diff_pd[index_header][0] == '#q2:types':
                            diff_pd = diff_pd[1:]
                        diff_pd = diff_pd.rename(columns={index_header: 'Feature ID'}).set_index('Feature ID')
                        diff_pd = diff_pd.drop(columns=[x for x in diff_pd.columns if 'Intercept' in x])
                        diff_cols = ['%s__%s__Q2=%s' % (model, x, q2) for x in diff_pd.columns]
                        diff_pd.columns = diff_cols
                        feats_diff_cols.extend(diff_cols)
                        omic_diff_list.append(diff_pd)
            if len(omic_diff_list):
                omic_songbird_ranks = pd.concat(omic_diff_list, axis=1, sort=False).reset_index()
                omic_songbird_ranks.rename(columns={omic_songbird_ranks.columns[0]: 'Feature ID'}, inplace=True)
                # print('1.', omic_songbird_ranks.shape)
            else:
                omic_common_fp = mmvec_songbird_pd.loc[
                    (mmvec_songbird_pd['pair'] == pair) &
                    (mmvec_songbird_pd['omic%s' % omicn] == omic) &
                    (mmvec_songbird_pd['filt%s' % omicn] == filt),
                    'omic%s_common_fp' % omicn
                ].tolist()[0]
                omic_tax_list = []
                # print('2.', omic_common_fp)
                with open(omic_common_fp) as f:
                    for ldx, line in enumerate(f):
                        if ldx:
                            omic_tax_list.append([line.split('\t')[0]])
                omic_songbird_ranks = pd.DataFrame(omic_tax_list, columns=['Feature ID'])
                # print('2.', omic_songbird_ranks.shape)
            # print("omic")
            # print(omic)
            # print("taxo_pds.keys()")
            # print(taxo_pds.keys())
            if omic in taxo_pds:
                omic_tax_pd = taxo_pds[omic]
                if 'Taxon' in omic_tax_pd.columns:
                    omic_split_taxa_pd = get_split_taxonomy(omic_tax_pd.Taxon.tolist())
                    omic_tax_pd = pd.concat([omic_tax_pd, omic_split_taxa_pd], axis=1, sort=False)
                omic_songbird_ranks = omic_songbird_ranks.merge(
                    omic_tax_pd, on='Feature ID', how='left').drop_duplicates()
            # print('3.', omic_songbird_ranks.shape)
            meta_omic_fp = '%s/feature_metadata_%s__%s.tsv' % (cur_mmvec_folder, omic, filt)
            drop_columns = [col for col in omic_songbird_ranks.columns if omic_songbird_ranks[col].unique().size == 1]
            meta_omic_pd = omic_songbird_ranks.drop(columns=drop_columns)
            meta_omic_pd.to_csv(meta_omic_fp, index=False, sep='\t')
            # print('<< written: %s >>' % meta_omic_fp)
            # print('-' *50)
            meta_omic_pd.set_index('Feature ID', inplace=True)
            omics_pairs_metas[(pair, omic, filt)] = (meta_omic_fp, meta_omic_pd, feats_diff_cols)
    return omics_pairs_metas


def get_taxo_pds(i_datasets_folder, mmvec_songbird_pd, input_to_filtered):
    taxo_pds = {}
    for omicn in ['1', '2']:
        for omic in mmvec_songbird_pd['omic%s' % omicn].unique():
            omic_tax_fp = get_tax_fp(i_datasets_folder, omic, input_to_filtered)
            if isfile(omic_tax_fp):
                omic_tax_pd = pd.read_csv(omic_tax_fp, header=0, sep='\t', dtype=str)
                omic_tax_pd.rename(columns={omic_tax_pd.columns[0]: 'Feature ID'}, inplace=True)
            else:
                omic_tax_pd = pd.DataFrame()
            taxo_pds[omic] = omic_tax_pd
    return taxo_pds


def get_pc_sb_correlations(pair, ordi, omic1, omic2, filt1, filt2,
                           diff_cols1, meta_pd1, diff_cols2, meta_pd2,
                           meta_fp, omic1_common_fp, omic2_common_fp, ranks_fp):

    corrs = []
    for r in range(3):
        print("r")
        print(r)
        print("ordi.features.iloc[:3,:]")
        print(ordi.features.iloc[:3,:])
        feats = ordi.features[r]
        print("feats[:4]")
        print(feats[:4])
        if len(diff_cols1) > 1:
            print("diff_cols1")
            print(diff_cols1)
            for model in diff_cols1:
                print(' -', model)
                x = meta_pd1.loc[ordi.features.index, model].astype(float)
                x = x[x.notnull()]
                y = feats[x.index]
                r1, p1 = pearsonr(x, y)
                r2, p2 = spearmanr(x, y)
                # corrs.append([pair, omic1, filt1, 'PC%s' % (r+1), model, r1, p1, 'pearson',
                #              meta_fp, omic1_common_fp, ranks_fp])
                corrs.append([pair, omic1, filt1, 'PC%s' % (r + 1), model, r1, p2, 'spearman',
                              meta_fp,  omic1_common_fp, ranks_fp])
        sams = ordi.samples[r]
        if len(diff_cols2) > 1:
            for model in diff_cols2:
                x = meta_pd2.loc[ordi.samples.index, model].astype(float)
                x = x[x.notnull()]
                y = sams[x.index]
                r1, p1 = pearsonr(x, y)
                r2, p2 = spearmanr(x, y)
                # corrs.append([pair, omic2, filt2, 'PC%s' % (r+1), model, r1, p1, 'pearson',
                #               meta_fp, omic2_common_fp, ranks_fp])
                corrs.append([pair, omic2, filt2, 'PC%s' % (r + 1), model, r1, p2, 'spearman',
                              meta_fp, omic2_common_fp, ranks_fp])
    corrs_pd = pd.DataFrame(corrs, columns=[
        'pair',
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
    return corrs_pd


def get_log_ratios(pc_sb_correlations_pd, correlation_threshold = 0.75):
    print(' - Checking log_ratio bins for features correlating to mmvec PCs')
    pc_sb_correlations_thresh_pd = pc_sb_correlations_pd.loc[
        ((pc_sb_correlations_pd['correlation_coefficient']).apply(abs) > correlation_threshold) &
        (pc_sb_correlations_pd['pvalue'] < 0.05)
    ]
    if pc_sb_correlations_thresh_pd.shape[0]:
        print('    > found:')
        for i in [['pair'], ['omic', 'filt'], ['model']]:
            print('       * %s for %s' % (pc_sb_correlations_thresh_pd.groupby(i).count().shape[0], i))
        for pair, pair_pd in pc_sb_correlations_thresh_pd.groupby('pair'):
            print()
            print(pair)
            for row in pair_pd.values:
                pair, omic, filt, mmvec_pc, model = row[:5]
                correlation_coefficient, pvalue, correlation_method, = row[5:8]
                meta_fp, features_fp, ranks_fp = row[-3:]


def summarize_songbirds(i_datasets_folder) -> pd.DataFrame:
    q2s = []
    songbird_folder = get_analysis_folder(i_datasets_folder, 'songbird')
    for root, dirs, files in os.walk(songbird_folder):
        for fil in files:
            if fil.endswith('.html'):
                path = root + '/' + fil
                diff = path.replace('-tensorboard.html', '.tsv')
                root_split = root.split('%s/' % songbird_folder)[-1].split('/')
                dat = root_split[0]
                if len(root_split) == 7:
                    pair, dataset_filter, subset, model, songbird_filter, parameters = root_split[1:]
                else:
                    pair = 'no_pair'
                    dataset_filter, subset, model, songbird_filter, parameters = root_split[1:]
                with open(path) as f:
                    for line in f:
                        if 'Pseudo Q-squared' in line:
                            q2s.append([
                                pair, dat, dataset_filter, subset, model, songbird_filter, parameters, diff,
                                float(line.split('Pseudo Q-squared:</a></strong> ')[-1].split('<')[0])
                            ])
    q2s_pd = pd.DataFrame(q2s, columns=['pair', 'dat', 'dataset_filter', 'subset',
                                        'model', 'songbird_filter', 'parameters',
                                        'differentials', 'Pseudo_Q_squared'])
    return q2s_pd


def run_mmbird(i_datasets_folder: str, songbird_outputs: list, p_mmvec_highlights: str,
               mmvec_outputs: list, force: bool, prjct_nm: str,
               qiime_env: str, chmod: str, noloc: bool, filt_raref: str,
               run_params: dict, input_to_filtered: dict) -> pd.DataFrame:

    if not mmvec_outputs:
        print('No mmvec output detected...')
        sys.exit(0)
    if not songbird_outputs:
        print('No songbird output detected...')
        sys.exit(0)

    print('\t-> [mmbird] Get mmvec output...', end=' ')
    mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
    print('Done.')

    print('\t-> [mmbird] Get songbird output...', end=' ')
    if len(songbird_outputs):
        songbird_outputs_pd = get_songbird_outputs(songbird_outputs)
        mmvec_songbird_pd = merge_mmvec_songbird_outputs(mmvec_outputs_pd, songbird_outputs_pd)
    else:
        mmvec_songbird_pd = mmvec_outputs_pd.copy()
    print('Done.')

    q2s_pd = summarize_songbirds(i_datasets_folder)
    out_folder = get_analysis_folder(i_datasets_folder, 'songbird')
    q2s_fp = '%s/songbird_q2.tsv' % out_folder
    q2s_pd.to_csv(q2s_fp, index=False, sep='\t')
    print(' --> Written:', q2s_fp)

    omics_pairs = [tuple(x) for x in mmvec_songbird_pd[['omic_filt1', 'omic_filt2']].values.tolist()]

    print('\t-> [mmbird] Get taxonomy (feature metadata)...', end=' ')
    taxo_pds = get_taxo_pds(
        i_datasets_folder, mmvec_songbird_pd, input_to_filtered)
    print('Done.')

    print('\t-> [mmbird] Get songbird differentials + taxonomy...', end=' ')
    omics_pairs_metas = get_omics_songbirds_taxa(
        i_datasets_folder, mmvec_songbird_pd, taxo_pds)
    print('Done.')

    print('\t-> [mmbird] Get res dict...', end=' ')
    mmvec_res = get_mmvec_res(mmvec_songbird_pd)
    print('Done.')

    print('\t-> [mmbird] Get commands...', end=' ')
    highlights = get_highlights(p_mmvec_highlights)
    pair_cmds, pc_sb_correlations_pd = get_pair_cmds(
        mmvec_res, omics_pairs_metas, omics_pairs, force, highlights)
    print('Done.')

    out_folder = get_analysis_folder(i_datasets_folder, 'mmbird')
    out_correlations = '%s/pc_vs_songbird_correlations.tsv' % out_folder
    pc_sb_correlations_pd.to_csv(out_correlations, index=False, sep='\t')
    print(' --> Written:', out_correlations)

    # get_log_ratios(pc_sb_correlations_pd)

    job_folder = get_job_folder(i_datasets_folder, 'mmbird')
    job_folder2 = get_job_folder(i_datasets_folder, 'mmbird/chunks')
    written = 0
    run_pbs = '%s/run_mmbird%s.sh' % (job_folder, filt_raref)
    with open(run_pbs, 'w') as o:
        for pair, cmds in pair_cmds.items():
            out_sh = '%s/run_mmbird_%s%s.sh' % (job_folder2, pair, filt_raref)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for cmd in cmds:
                    cur_sh.write(cmd)
                    written += 1
            run_xpbs(out_sh, out_pbs, '%s.mmbrd.%s%s' % (prjct_nm, pair, filt_raref), qiime_env,
                     run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                     run_params["mem_num"], run_params["mem_dim"],
                     chmod, written, 'single', o, noloc)
    if written:
        print_message('# Generate mmvec biplot with songbird models', 'sh', run_pbs)
    return pc_sb_correlations_pd
