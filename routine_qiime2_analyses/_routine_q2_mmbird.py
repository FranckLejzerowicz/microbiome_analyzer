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
import time
import pandas as pd
import multiprocessing as mp
from os.path import dirname, isfile, splitext

from scipy.stats import spearmanr
from skbio.stats.ordination import OrdinationResults

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder, get_highlights
from routine_qiime2_analyses._routine_q2_cmds import get_split_taxonomy


def get_mmvec_outputs(mmvec_outputs: list):
    # print('-'*30)
    # for i in mmvec_outputs:
    #     print()
    #     print(i)
    # print('-'*30)
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


def get_songbird_outputs(songbird_outputs: list) -> pd.DataFrame:
    songbird_outputs_pd = pd.DataFrame(
        songbird_outputs,
        columns=[
            'songbird_dat',
            'songbird_filt',
            'songbird_parameters',
            'songbird_case',
            'songbird_fp',
            'songbird_baseline',
            'songbird_q2',
            'pair'
        ])
    songbird_outputs_pd['pair_omic_filt'] = songbird_outputs_pd['pair'] + '__' + songbird_outputs_pd[
        'songbird_dat'] + '__' + songbird_outputs_pd['songbird_filt']
    songbird_outputs_pd['case_params'] = songbird_outputs_pd['songbird_case'] + '__' + songbird_outputs_pd[
        'songbird_parameters']
    songbird_outputs_pd['case_params_baseline'] = songbird_outputs_pd[
        'case_params'] + '__' + songbird_outputs_pd['songbird_baseline']
    songbird_outputs_drop_pd = songbird_outputs_pd.drop(
        columns=['pair', 'songbird_dat', 'songbird_filt', 'songbird_case',
                 'songbird_parameters', 'songbird_baseline'])

    songbird_outputs_pd = songbird_outputs_drop_pd[
        ['case_params', 'pair_omic_filt', 'songbird_fp']
    ].drop_duplicates().pivot(
        columns='case_params', index='pair_omic_filt'
    )
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
        if '__raref' in omic_tax_:
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
        print()
        print(pair, omic1, omic2, filt1, filt2, sams)
        print(mmvec)
        ranks_fp, ordi_fp, meta_fp, omic1_common_fp, omic2_common_fp = values

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
        meta1, meta_pd1, diff_cols1 = omics_pairs_metas[(pair, omic1, filt1, omic2, filt2)]
        meta2, meta_pd2, diff_cols2 = omics_pairs_metas[(pair, omic2, filt2, omic1, filt1)]

        # features are biplot, samples are dots
        ordi = OrdinationResults.read(ordi_fp)

        start = time.time()
        cur_pc_sb_correlations = get_pc_sb_correlations(
            pair, ordi, omic1, omic2, filt1, filt2,
            diff_cols1, meta_pd1, diff_cols2, meta_pd2,
            meta_fp,  omic1_common_fp, omic2_common_fp, ranks_fp)
        pc_sb_correlations.append(cur_pc_sb_correlations)
        end = time.time()
        print('get_pc_sb_correlations: %s' % (end-start))

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
                n_ordi_feats = 15
                # heat_qza, heat_qzv = get_heatmap_qzs(ranks_fp)
                # cmd += get_heatmap_commands(
                #     ranks_fp, heat_qza, heat_qzv, meta1,
                #     meta2, meta_pd1, meta_pd2)
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
        pair_omics_filts = ['pair', 'omic1', 'filt1', 'omic2', 'filt2']
        all_omic_sb = [x for x in mmvec_songbird_pd.columns if x.endswith('omic%s_songbird_common_fp' % omicn)]
        omicn_songbirds = mmvec_songbird_pd[(pair_omics_filts + all_omic_sb)].set_index(pair_omics_filts).T.to_dict()
        for (pair, omic1, filt1, omic2, filt2), sb_head_diff_fp in omicn_songbirds.items():
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
            cur_mmvec_folder = get_analysis_folder(i_datasets_folder, 'mmvec/metadata/%s' % pair)
            omic_diff_list = []
            if len(sb_head_diff_fp):
                for sb_head, diff_fp in sb_head_diff_fp.items():
                    model = sb_head.replace('_omic%s_songbird_common_fp' % omicn, '')
                    if str(diff_fp) != 'nan' and isfile(diff_fp):

                        diff_pd = pd.read_csv(diff_fp, header=0, sep='\t', dtype=str)
                        index_header = diff_pd.columns[0]
                        if diff_pd[index_header][0] == '#q2:types':
                            diff_pd = diff_pd[1:]
                        diff_pd = diff_pd.rename(columns={index_header: 'Feature ID'}).set_index('Feature ID')
                        diff_pd = diff_pd.drop(columns=[x for x in diff_pd.columns if 'Intercept' in x])

                        q2s = {}
                        diff_htmls = glob.glob('%s/*/tensorboard.html' % dirname(diff_fp))
                        if len(diff_htmls):
                            for diff_html in diff_htmls:
                                baseline = diff_html.split('/')[-2]
                                with open(diff_html) as f:
                                    for line in f:
                                        if 'Pseudo Q-squared' in line:
                                            q2 = line.split('Pseudo Q-squared:</a></strong> ')[-1].split('<')[0]
                                            if float(q2) > 0.1:
                                                q2s[baseline] = q2
                                                break
                        else:
                            q2s[''] = ''
                        diff_cols = ['%s__%s__%s' % (
                            model, x, '--'.join(['%s-Q2=%s' % (b, q) if b else 'noQ2' for b, q in q2s.items()])
                        ) for x in diff_pd.columns]
                        diff_pd.columns = diff_cols
                        feats_diff_cols.extend(diff_cols)
                        omic_diff_list.append(diff_pd)
            if len(omic_diff_list):
                omic_songbird_ranks = pd.concat(omic_diff_list, axis=1, sort=False).reset_index()
                omic_songbird_ranks.rename(columns={omic_songbird_ranks.columns[0]: 'Feature ID'}, inplace=True)
                # print("'\n'.join(omic_songbird_ranks.columns.tolist())")
                # print('\n'.join(omic_songbird_ranks.columns.tolist()))
                # print('1.', omic_songbird_ranks.shape)
            else:
                omic_common_fp = mmvec_songbird_pd.loc[
                    (mmvec_songbird_pd['pair'] == pair) &
                    (mmvec_songbird_pd['omic1'] == omic1) &
                    (mmvec_songbird_pd['filt1'] == filt1) &
                    (mmvec_songbird_pd['omic2'] == omic2) &
                    (mmvec_songbird_pd['filt2'] == filt2),
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
                if omic_tax_pd.shape[0]:
                    if 'Taxon' in omic_tax_pd.columns:
                        omic_split_taxa_pd = get_split_taxonomy(omic_tax_pd.Taxon.tolist())
                        omic_tax_pd = pd.concat([omic_tax_pd, omic_split_taxa_pd], axis=1, sort=False)
                    omic_songbird_ranks = omic_songbird_ranks.merge(
                        omic_tax_pd, on='Feature ID', how='left').drop_duplicates()
            # print('3.', omic_songbird_ranks.shape)
            meta_omic_fp = '%s/feature_metadata_%s_%s__%s_%s.tsv' % (cur_mmvec_folder, omic, filt, omic_, filt_)
            drop_columns = [col for col in omic_songbird_ranks.columns if omic_songbird_ranks[col].unique().size == 1]
            meta_omic_pd = omic_songbird_ranks.drop(columns=drop_columns)
            meta_omic_pd.to_csv(meta_omic_fp, index=False, sep='\t')
            # print('<< written: %s >>' % meta_omic_fp)
            # print('-' *50)
            meta_omic_pd.set_index('Feature ID', inplace=True)
            omics_pairs_metas[(pair, omic, filt, omic_, filt_)] = (meta_omic_fp, meta_omic_pd, feats_diff_cols)
    return omics_pairs_metas


def get_taxo_pds(i_datasets_folder, mmvec_songbird_pd, input_to_filtered):
    taxo_pds = {}
    for omicn in ['1', '2']:
        # print('omicn:', omicn)
        # print(mmvec_songbird_pd['omic%s' % omicn].unique())
        for omic in mmvec_songbird_pd['omic%s' % omicn].unique():
            omic_tax_fp = get_tax_fp(i_datasets_folder, omic, input_to_filtered)
            # print("omic_tax_fp")
            # print(omic_tax_fp)
            if isfile(omic_tax_fp):
                omic_tax_pd = pd.read_csv(omic_tax_fp, header=0, sep='\t', dtype=str)
                omic_tax_pd.rename(columns={omic_tax_pd.columns[0]: 'Feature ID'}, inplace=True)
            else:
                omic_tax_pd = pd.DataFrame()
            # print("omic_tax_pd.iloc[:3,:3]")
            # print(omic_tax_pd.iloc[:3, :3])
            taxo_pds[omic] = omic_tax_pd
    return taxo_pds


def get_mp_corrs(args):
    corrs, r, feats_sams, diff_cols, meta_pd, pair = args[:6]
    omic, filt, meta_fp, omic_common_fp, ranks_fp = args[6:]
    if len(diff_cols):
        for model in diff_cols:
            x = meta_pd.loc[
                [x for x in meta_pd.index if x in feats_sams.index], model
            ].astype(float)
            x = x[x.notnull()]
            y = feats_sams[x.index]
            r2, p2 = spearmanr(x, y)
            corrs.append([pair, omic, filt, 'PC%s' % (r + 1), model, r2, p2, 'spearman',
                          meta_fp, omic_common_fp, ranks_fp])


def get_pc_sb_correlations(pair, ordi, omic1, omic2, filt1, filt2,
                           diff_cols1, meta_pd1, diff_cols2, meta_pd2,
                           meta_fp, omic1_common_fp, omic2_common_fp, ranks_fp):
    do_mp = False
    if do_mp:
        manager = mp.Manager()
        corrs = manager.list()
        pool = mp.Pool(processes=6)
        pool_runs = []
    else:
        corrs = []

    for r in range(3):
        feats = ordi.features[r]
        if do_mp:
            pool_runs.append((corrs, r, feats, diff_cols1, meta_pd1, pair, omic1, filt1, meta_fp,  omic1_common_fp, ranks_fp))
        else:
            if len(diff_cols1):
                for model in diff_cols1:
                    x = meta_pd1.loc[
                        [x for x in meta_pd1.index if x in feats.index], model
                    ].astype(float)
                    x = x[x.notnull()]
                    y = feats[x.index]
                    r2, p2 = spearmanr(x, y)
                    corrs.append([pair, omic1, filt1, 'PC%s' % (r + 1), model, r2, p2, 'spearman',
                                  meta_fp,  omic1_common_fp, ranks_fp])
        sams = ordi.samples[r]
        if do_mp:
            pool_runs.append((corrs, r, sams, diff_cols2, meta_pd2, pair, omic2, filt2, meta_fp,  omic1_common_fp, ranks_fp))
        else:
            if len(diff_cols2):
                for model in diff_cols2:
                    x = meta_pd2.loc[
                        [x for x in meta_pd2.index if x in sams.index], model
                    ].astype(float)
                    x = x[x.notnull()]
                    y = sams[x.index]
                    r2, p2 = spearmanr(x, y)
                    corrs.append([pair, omic2, filt2, 'PC%s' % (r + 1), model, r2, p2, 'spearman',
                                  meta_fp, omic2_common_fp, ranks_fp])
    if do_mp:
        pool.map(get_mp_corrs, [x for x in pool_runs])
        pool.close()
        pool.join()

    corrs_pd = pd.DataFrame(list(corrs), columns=[
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
            if fil == 'tensorboard.html':
                path = root + '/' + fil
                diff = '%s/differentials.tsv' % dirname(root)
                root_split = root.split('%s/' % songbird_folder)[-1].split('/')
                if len(root_split) == 8:
                    dat, pair, dataset_filter, subset, songbird_filter, parameters, model, baseline = root_split
                else:
                    pair = 'no_pair'
                    dat, dataset_filter, subset, songbird_filter, parameters, model, baseline = root_split
                with open(path) as f:
                    for line in f:
                        if 'Pseudo Q-squared' in line:
                            q2s.append([
                                pair, dat, dataset_filter, subset, model, songbird_filter, parameters, baseline, diff,
                                float(line.split('Pseudo Q-squared:</a></strong> ')[-1].split('<')[0])
                            ])
    q2s_pd = pd.DataFrame(q2s, columns=['pair', 'dat', 'dataset_filter', 'subset', 'model',
                                        'songbird_filter', 'parameters', 'baseline',
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
        # songbird_outputs_pd = get_songbird_outputs(songbird_outputs)
        songbird_outputs_pd = get_songbird_outputs(songbird_outputs)
        mmvec_songbird_pd = merge_mmvec_songbird_outputs(mmvec_outputs_pd, songbird_outputs_pd)
    else:
        mmvec_songbird_pd = mmvec_outputs_pd.copy()
    print('Done.')

    q2s_pd = summarize_songbirds(i_datasets_folder)
    out_folder = get_analysis_folder(i_datasets_folder, 'songbird')
    q2s_fp = '%s/songbird_q2.tsv' % out_folder
    q2s_pd.to_csv(q2s_fp, index=False, sep='\t')
    print('\t\t==> Written:', q2s_fp)

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

    print('\t-> [mmbird] Get commands...')
    highlights = get_highlights(p_mmvec_highlights)
    pair_cmds, pc_sb_correlations_pd = get_pair_cmds(
        mmvec_res, omics_pairs_metas, omics_pairs, force, highlights)

    out_folder = get_analysis_folder(i_datasets_folder, 'mmbird')
    out_correlations = '%s/pc_vs_songbird_correlations.tsv' % out_folder
    pc_sb_correlations_pd.to_csv(out_correlations, index=False, sep='\t')
    print('\t\t==> Written:', out_correlations)

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
