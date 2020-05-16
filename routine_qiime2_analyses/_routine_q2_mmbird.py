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

from skbio.stats.ordination import OrdinationResults

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder
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


def get_mmvec_res(mmvec_outputs_pd):
    mmvec_out_cols = [x for x in mmvec_outputs_pd.columns if x.startswith('mmvec_out__')]

    all_omic1_sb = [x for x in mmvec_outputs_pd.columns if x.endswith('omic1_songbird_common_fp')]
    all_omic2_sb = [x for x in mmvec_outputs_pd.columns if x.endswith('omic2_songbird_common_fp')]

    mmvec_res = {}
    # for ech row of the main table that also contain the mmvec output folders
    for r, row in mmvec_outputs_pd.iterrows():
        # get the omics
        pair = row['pair']
        omic1 = row['omic1']
        omic2 = row['omic2']
        omic_filt1 = row['omic_filt1']
        omic_filt2 = row['omic_filt2']
        omic1_common_fp = row['omic1_common_fp']
        if str(omic1_common_fp) == 'nan':
            continue
        omic2_common_fp = row['omic2_common_fp']
        n_common = row['n_common']
        meta_fp = row['meta_common_fp']

        # get the songbirds
        omic1_songbird_common_fps = [(x, all_omic1_sb[idx].replace('_omic1_songbird_common_fp', '')) for idx, x in enumerate(row[all_omic1_sb])]
        omic2_songbird_common_fps = [(x, all_omic2_sb[idx].replace('_omic2_songbird_common_fp', '')) for idx, x in enumerate(row[all_omic2_sb])]

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
                print('[ - mmvec] %s %s' % (omic_filt1, omic_filt2))
                continue

            # collect the ranks + ordination + songbirds for each pair of omics and parameters
            mmvec_res[
                (
                    pair,
                    omic1,
                    omic2,
                    omic_filt1,
                    omic_filt2,
                    n_common,
                    mmvec_out_col.replace('mmvec_out__', '')
                )
            ] = [
                mmvec_out_ranks,
                mmvec_out_ordi,
                omic1_songbird_common_fps,
                omic2_songbird_common_fps,
                meta_fp,
                omic1_common_fp,
                omic2_common_fp
            ]

    return mmvec_res


def get_all_omics_songbirds(omic1_diff_fps, omic2_diff_fps):


    all_omic1_diff_list = []
    all_omic1_songbird_ranks = pd.DataFrame()
    for (omic1_diff_fp, model) in omic1_diff_fps:
        if isfile(omic1_diff_fp):
            omic1_diff_pd = pd.read_csv(omic1_diff_fp, header=0, sep='\t')
            if omic1_diff_pd['featureid'][0] == '#q2:types':
                omic1_diff_pd = omic1_diff_pd[1:]
            omic1_diff_pd = omic1_diff_pd.rename(columns={'featureid': 'Feature ID'})
            omic1_diff_pd = omic1_diff_pd.set_index('Feature ID')
            omic1_diff_pd = omic1_diff_pd.drop(columns='Intercept')
            omic1_diff_pd.columns = ['%s__%s' % (model, x) for x in omic1_diff_pd.columns]
            all_omic1_diff_list.append(omic1_diff_pd)
    if len(all_omic1_diff_list):
        all_omic1_diff_pd = pd.concat(all_omic1_diff_list, axis=1, sort=False)
        all_omic1_songbird_ranks = all_omic1_diff_pd.reset_index()
        all_omic1_songbird_ranks = all_omic1_songbird_ranks.rename(
            columns={all_omic1_songbird_ranks.columns.tolist()[0]: 'Feature ID'})

    all_omic2_diff_list = []
    all_omic2_songbird_ranks = pd.DataFrame()
    for (omic2_diff_fp, model) in omic2_diff_fps:
        if isfile(omic2_diff_fp):
            omic2_diff_pd = pd.read_csv(omic2_diff_fp, header=0, sep='\t')
            if omic2_diff_pd['featureid'][0] == '#q2:types':
                omic2_diff_pd = omic2_diff_pd[1:]
            omic2_diff_pd = omic2_diff_pd.rename(columns={'featureid': 'Feature ID'})
            omic2_diff_pd = omic2_diff_pd.set_index('Feature ID')
            omic2_diff_pd = omic2_diff_pd.drop(columns='Intercept')
            omic2_diff_pd.columns = ['%s__%s' % (model, x) for x in omic2_diff_pd.columns]
            all_omic2_diff_list.append(omic2_diff_pd)
    if len(all_omic2_diff_list):
        all_omic2_diff_pd = pd.concat(all_omic2_diff_list, axis=1, sort=False)
        all_omic2_songbird_ranks = all_omic2_diff_pd.reset_index()
        all_omic2_songbird_ranks = all_omic2_songbird_ranks.rename(
            columns={all_omic2_songbird_ranks.columns.tolist()[0]: 'Feature ID'})
    return all_omic1_songbird_ranks, all_omic2_songbird_ranks


def get_qzs(ordi_fp):
    qza = ordi_fp.replace('.txt', '.qza')
    qzv = ordi_fp.replace('.txt', '_emperor.qzv')
    return qza, qzv


def get_heatmap_qzs(ranks_fp):
    qza = ranks_fp.replace('.csv', '.qza')
    qzv = ranks_fp.replace('.csv', '_heatmap.qzv')
    return qza, qzv


# def get_heatmap_qzs(ranks_fp, metatax_omic2_pd_col):
#     qza = ranks_fp.replace('.csv', '_%s.qza' % metatax_omic2_pd_col)
#     qzv = ranks_fp.replace('.csv', '_%s_heatmap.qzv' % metatax_omic2_pd_col)
#     return qza, qzv


def get_order_omics(omic1, omic2, omic_filt1, omic_filt2, omics_pairs):
    omic_feature, omic_sample = ('feature', 'sample')
    omic_microbe, omic_metabolite = ('microbe', 'metabolite')
    if (omic_filt2, omic_filt1) not in omics_pairs:
        omic_feature, omic_sample = ('sample', 'feature')
        omic_microbe, omic_metabolite = ('metabolite', 'microbe')
        omic1, omic2 = omic2, omic1
        omic_filt1, omic_filt2 = omic_filt2, omic_filt1
    return omic1, omic2, omic_filt1, omic_filt2, omic_feature, omic_sample, omic_microbe, omic_metabolite


def get_tax_extended_fps(
        omic_filt1, omic_filt2,
        omic1_tax_fp, omic2_tax_fp,
        all_omic1_songbird_ranks,
        all_omic2_songbird_ranks,
        ordi_fp
):
    omic1_tax_pd = pd.read_csv(omic1_tax_fp, header=0, sep='\t', dtype=str)
    omic2_tax_pd = pd.read_csv(omic2_tax_fp, header=0, sep='\t', dtype=str)
    if 'Taxon' in omic1_tax_pd:
        omic1_split_taxo = get_split_taxonomy(omic1_tax_pd.Taxon.tolist())
        omic1_tax_pd = omic1_tax_pd.merge(omic1_split_taxo, on='Taxon', how='left').drop_duplicates()
    if 'Taxon' in omic2_tax_pd:
        omic2_split_taxo = get_split_taxonomy(omic2_tax_pd.Taxon.tolist())
        omic2_tax_pd = omic2_tax_pd.merge(omic2_split_taxo, on='Taxon', how='left').drop_duplicates()

    if all_omic1_songbird_ranks.shape[0]:
        omic1_tax_pd = omic1_tax_pd.merge(
            all_omic1_songbird_ranks,
            on='Feature ID',
            how='left'
        ).drop_duplicates()

    if all_omic2_songbird_ranks.shape[0]:
        omic2_tax_pd = omic2_tax_pd.merge(
            all_omic2_songbird_ranks,
            on='Feature ID',
            how='left'
        ).drop_duplicates()

    metatax_omic1_fp = ordi_fp.replace('.txt', '_meta-%s.tsv' % omic_filt1)
    metatax_omic2_fp = ordi_fp.replace('.txt', '_meta-%s.tsv' % omic_filt2)
    omic1_tax_pd.to_csv(metatax_omic1_fp, index=False, sep='\t')
    omic2_tax_pd.to_csv(metatax_omic2_fp, index=False, sep='\t')
    return metatax_omic1_fp, metatax_omic2_fp


def get_biplot_commands(ordi_edit_fp, qza, qzv,
                        omic_feature, omic_sample,
                        metatax_omic1_fp, metatax_omic2_fp,
                        edit, n_mbAnnot_CLAs, crowded, max_feats):

    if max_feats == 0:
        ordi = OrdinationResults.read(ordi_edit_fp)
        max_feats = ordi.features.shape[0]

    cmd = '\n'
    if not isfile(qza):
        cmd += '\nqiime tools import'
        cmd += ' --input-path %s' % ordi_edit_fp
        cmd += ' --output-path %s' % qza
        cmd += ' --type "PCoAResults %s Properties([\'biplot\'])"\nsleep 3' % '%'

    cmd += '\nqiime emperor biplot'
    cmd += ' --i-biplot %s' % qza
    cmd += ' --m-%s-metadata-file %s' % (omic_feature, metatax_omic1_fp)
    cmd += ' --m-%s-metadata-file %s' % (omic_sample, metatax_omic2_fp)
    if edit:
        cmd += ' --p-number-of-features %s' % n_mbAnnot_CLAs
        cmd += ' --o-visualization %s\n' % qzv
    elif crowded:
        cmd += ' --p-number-of-features %s' % max_feats
        cmd += ' --o-visualization %s\n' % qzv.replace('.qzv', '_crowded.qzv')
    else:
        cmd += ' --o-visualization %s\n' % qzv

    return cmd


# def get_heatmap_commands(ranks_edit_fp, qza, qzv,
#                          omic_microbe, omic_metabolite,
#                          # metatax_omic2_pd_col,
#                          metatax_omic1_fp, metatax_omic2_fp):
#     cmd = '\n\nqiime tools import'
#     cmd += ' --input-path %s' % ranks_edit_fp
#     cmd += ' --output-path %s' % qza
#     cmd += ' --type FeatureData[Conditional]'
#     cmd += ' --input-format ConditionalFormat\n'
#
#     cmd += '\nqiime mmvec heatmap'
#     cmd += ' --i-ranks %s' % qza
#     if omic_microbe == 'metabolite':
#         cmd += ' --m-%s-metadata-file %s' % (omic_microbe, metatax_omic2_fp)
#         cmd += ' --m-%s-metadata-column Node_level_D' % omic_microbe
#         # cmd += ' --m-%s-metadata-file %s \\ \n' % (omic_metabolite, metatax_omic2_fp)
#         # cmd += ' --m-%s-metadata-column "" \\ \n' % omic_metabolite
#         # cmd += ' --m-%s-metadata-column "%s" \\ \n' % (omic_metabolite, metatax_omic2_pd_col)
#     else:
#         # cmd += ' --m-%s-metadata-file %s \\ \n' % (omic_microbe, metatax_omic1_fp)
#         # cmd += ' --m-%s-metadata-column "%s" \\ \n' % omic_microbe
#         # cmd += ' --m-%s-metadata-column "%s" \\ \n' % (omic_microbe, metatax_omic2_pd_col)
#         cmd += ' --m-%s-metadata-file %s' % (omic_metabolite, metatax_omic1_fp)
#         cmd += ' --m-%s-metadata-column Node_level_D' % omic_metabolite
#     # cmd += ' --p-level 5 \\ \n'
#     cmd += ' --o-visualization %s\n' % qzv
#     return cmd


# def get_n_mbAnnot_CLAs_in_file(ordi_fp, mbAnnot_CLAs):
#     n_mbAnnot_CLAs_in_file = 0
#     with open(ordi_fp) as ordi_f:
#         is_mb = 0
#         for line in ordi_f:
#             if line.startswith('Species'):
#                 is_mb = 1
#                 continue
#             if line.startswith('Site'):
#                 is_mb = 0
#                 continue
#             if is_mb:
#                 for mbAnnot_CLA in mbAnnot_CLAs:
#                     if mbAnnot_CLA in line:
#                         n_mbAnnot_CLAs_in_file += 1
#                         break
#     return n_mbAnnot_CLAs_in_file
#
#
# def write_edit_ordi(edit_fp, ordi_edit_fp, n_mbAnnot_CLAs_in_file, mbAnnot_CLAs):
#     with open(ordi_edit_fp, 'w') as ordi_fo, open(ordi_fp) as ordi_f:
#         is_mb = 0
#         for line in ordi_f:
#             if line.startswith('Site\t'):
#                 ordi_fo.write('\n')
#                 is_mb = 0
#             if line.startswith('Species'):
#                 ordi_fo.write('Species\t%s\t%s\n' % (
#                     n_mbAnnot_CLAs_in_file,
#                     line.strip().split('\t')[-1]
#                 ))
#                 is_mb = 1
#                 continue
#             if is_mb:
#                 for mbAnnot_CLA in mbAnnot_CLAs:
#                     if mbAnnot_CLA in line:
#                         ordi_fo.write(line)
#                         break
#             else:
#                 ordi_fo.write(line)
#
#
# def get_ordi_edit_qzv_fp(ordi_fp, edit, omic1_common_fp, omic2_common_fp):
#     ordi_edit_fp = '%s%s%s' % (
#         splitext(ordi_fp)[0], edit,
#         splitext(ordi_fp)[1]
#     )
#     n_mbAnnot_CLAs_in_file = 0
#     if len(edit):
#         mbAnnot_fp = [x for x in [omic1_common_fp, omic2_common_fp] if 'mbAnnot' in x][0]
#         mbAnnot_pd = pd.read_csv(mbAnnot_fp, header=0, index_col=0, sep='\t')
#         mbAnnot_CLAs = [x for x in mbAnnot_pd.index if 'inoleic' in x]
#         n_mbAnnot_CLAs_in_file = get_n_mbAnnot_CLAs_in_file(ordi_fp, mbAnnot_CLAs)
#         write_edit_ordi(ordi_fp, ordi_edit_fp, n_mbAnnot_CLAs_in_file, mbAnnot_CLAs)
#     return n_mbAnnot_CLAs_in_file, ordi_edit_fp


def get_pair_cmds(mmvec_songbird_pd, mmvec_res, taxonomies, force):

    crowdeds = [0, 1]
    omics_pairs = [tuple(x) for x in mmvec_songbird_pd[['omic_filt1', 'omic_filt2']].values.tolist()]

    mmvec_tab = []
    pair_cmds = {}
    for crowded in crowdeds:
        for keys, values in mmvec_res.items():

            pair, omic1, omic2, omic_filt1, omic_filt2, sams, mmvec = keys

            ranks_fp = values[0]
            ordi_fp = values[1]
            omic1_diff_fps = values[2]
            omic2_diff_fps = values[3]
            meta_fp = values[4]
            omic1_common_fp = values[5]
            omic2_common_fp = values[6]

            omic1, omic2, omic_filt1, omic_filt2, omic_feature, omic_sample, omic_microbe, omic_metabolite = get_order_omics(
                omic1, omic2, omic_filt1, omic_filt2, omics_pairs)

            # get differentials
            all_omic1_songbird_ranks, all_omic2_songbird_ranks = get_all_omics_songbirds(omic1_diff_fps, omic2_diff_fps)
            omic1_tax_fp = '%s.tsv' % splitext(taxonomies[omic1][1])[0]
            omic2_tax_fp = '%s.tsv' % splitext(taxonomies[omic2][1])[0]
            metatax_omic1_fp, metatax_omic2_fp = get_tax_extended_fps(
                omic_filt1, omic_filt2,
                omic1_tax_fp, omic2_tax_fp,
                all_omic1_songbird_ranks,
                all_omic2_songbird_ranks,
                ordi_fp
            )

            if all_omic1_songbird_ranks.shape[0]:
                max_feats = all_omic1_songbird_ranks.shape[0]
            else:
                max_feats = 0

            mmvec_tab.append([
                omic_filt1, omic_filt2,
                '%s-%s-%s' % (omic_filt1, omic_filt2, sams),
                mmvec,
                ranks_fp,
                ordi_fp,
                meta_fp,
                omic1_tax_fp, omic2_tax_fp,
                omic1_common_fp, omic2_common_fp
            ])

            # if 'mbAnnot' in ordi_fp:
            #     for edit in ['', '_CLAs']:
            #         n_mbAnnot_CLAs_in_file, ordi_edit_fp = get_ordi_edit_qzv_fp(
            #             ordi_fp, edit,
            #             omic1_common_fp,
            #             omic2_common_fp
            #         )
            #         qza, qzv = get_qzs(ordi_edit_fp)
            #         if not isfile(qzv):
            #             cmd_import, cmd_biplot = get_biplot_commands(
            #                 ordi_edit_fp, qza, qzv,
            #                 omic_feature, omic_sample,
            #                 metatax_omic1_fp, metatax_omic2_fp,
            #                 edit, n_mbAnnot_CLAs_in_file
            #             )
            #             all_cmds.setdefault(key, []).extend([cmd_import, cmd_biplot])
            # else:
            n_mbAnnot_CLAs_in_file = 0
            ordi_edit_fp = ordi_fp
            qza, qzv = get_qzs(ordi_edit_fp)
            # if 1:
            cmd = ''
            if force or not isfile(qzv):
                cmd = get_biplot_commands(
                    ordi_edit_fp, qza, qzv,
                    omic_feature, omic_sample,
                    metatax_omic1_fp, metatax_omic2_fp,
                    '', n_mbAnnot_CLAs_in_file, crowded, max_feats
                )

            ranks_edit_fp = ranks_fp
            # if omic_metabolite == 'metabolite':
            #     metatax_omic2_pd = pd.read_csv(metatax_omic2_fp, header=0, sep='\t')
            #     metatax_omic2_pd_cols = metatax_omic2_pd.columns.tolist()[(metatax_omic2_pd.columns.tolist().index('Intercept')+1):]
            # for metatax_omic2_pd_col in metatax_omic2_pd_cols:
            #     metatax_omic2_col = ''.join([x for x in metatax_omic2_pd_col if x not in ["'", '(', ')', ' ', ',', '[', ']', '=', '"']])
            # qza, qzv = get_heatmap_qzs(ranks_edit_fp, metatax_omic2_col)
            # qza, qzv = get_heatmap_qzs(ranks_edit_fp)
            # if not crowded:
            #     cmd_heat = get_heatmap_commands(
            #         ranks_edit_fp, qza, qzv,
            #         omic_microbe, omic_metabolite,
            #         # metatax_omic2_pd_col,
            #         metatax_omic1_fp, metatax_omic2_fp
            #     )
            #     if cmd:
            #         cmd += cmd_heat
            #     else:
            #         cmd = cmd_heat

            if cmd:
                pair_cmds.setdefault(pair, []).append(cmd)
    return pair_cmds


def run_mmbird(i_datasets_folder: str, taxonomies: dict,
               songbird_outputs: list, mmvec_outputs: list, force: bool,
               prjct_nm: str, qiime_env: str, chmod: str, noloc: bool) -> None:

    if not mmvec_outputs:
        print('No mmvec output detected...')
        sys.exit(0)
    if not songbird_outputs:
        print('No songbird output detected...')
        sys.exit(0)

    songbird_outputs_pd = get_songbird_outputs(songbird_outputs)
    mmvec_outputs_pd = get_mmvec_outputs(mmvec_outputs)
    mmvec_songbird_pd = merge_mmvec_songbird_outputs(mmvec_outputs_pd, songbird_outputs_pd)
    mmvec_res = get_mmvec_res(mmvec_songbird_pd)

    pair_cmds = get_pair_cmds(mmvec_songbird_pd, mmvec_res, taxonomies, force)
    job_folder = get_job_folder(i_datasets_folder, 'mmbird')
    job_folder2 = get_job_folder(i_datasets_folder, 'mmbird/chunks')

    written = 0
    run_pbs = '%s/run_mmbird.sh' % job_folder
    with open(run_pbs, 'w') as o:
        for pair, cmds in pair_cmds.items():
            out_sh = '%s/run_mmbird_%s.sh' % (job_folder2, pair)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as cur_sh:
                for cmd in cmds:
                    cur_sh.write(cmd)
                    written += 1
            run_xpbs(out_sh, out_pbs, '%s.mmbrd.%s' % (prjct_nm, pair),
                     qiime_env, '24', '1', '1', '10', 'gb',
                     chmod, written, 'single', o, noloc)
    if written:
        print_message('# Generate mmvec biplot with sopngbird models', 'sh', run_pbs)

