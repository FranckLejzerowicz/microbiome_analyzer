# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, glob
import pandas as pd
# import numpy as np
import pkg_resources
from os.path import basename, isdir, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    read_yaml_file,
    read_meta_pd,
    write_main_sh
)
from routine_qiime2_analyses._routine_q2_cmds import (
    filter_feature_table,
    run_add_metadata,
    write_nestedness_graph,
    write_nestedness_nodfs,
    get_new_meta_pd,
    get_case, run_export
)
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict

# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# import seaborn as sns


def get_nestedness_config(nestedness_config: dict) -> (dict, list, dict, list, list):
    subsets = {'ALL': [[]]}
    if 'subsets' in nestedness_config:
        subsets.update(nestedness_config['subsets'])
    nodfs = []
    if 'nodfs' in nestedness_config:
        nodfs.extend(nestedness_config['nodfs'])
    colors = {'sample': [], 'feature': []}
    if 'sample_colors' in nestedness_config:
        colors['sample'].extend(nestedness_config['sample_colors'])
    if 'feature_colors' in nestedness_config:
        colors['feature'].extend(nestedness_config['feature_colors'])
    nulls = ['equiprobablefixed']
    if 'nulls' in nestedness_config:
        nulls = nestedness_config['nulls']
    modes = ['betweeneachpairoftypes']
    if 'modes' in nestedness_config:
        modes = nestedness_config['modes']
    return subsets, nodfs, colors, nulls, modes


def run_nestedness(i_datasets_folder: str, betas: dict, p_nestedness_groups: str,
                   datasets_rarefs: dict, force: bool, prjct_nm: str, qiime_env: str,
                   chmod: str, noloc: bool, split: bool, run_params: dict,
                   filt_raref: str, jobs: bool, chunkit: int) -> (dict, list):

    job_folder2 = get_job_folder(i_datasets_folder, 'nestedness/chunks')
    nestedness_config = read_yaml_file(p_nestedness_groups)
    if 'soft' not in nestedness_config:
        print('Must provide the path to the Nestedness soft (containing bin/Autocorrelation.jar)')
        return {}
    if nestedness_config['soft'].endswith('Autocorrelation.jar') and isfile(nestedness_config['soft']):
        binary = nestedness_config['soft']
    else:
        binary = '%s/bin/Autocorrelation.jar' % nestedness_config['soft']
        if not isfile(binary):
            print('Must provide the path to the Nestedness soft (containing bin/Autocorrelation.jar)')
            return {}
    subsets, nodfs, colors, nulls, modes = get_nestedness_config(nestedness_config)

    all_sh_pbs = {}
    nestedness_res = {}
    for dat, rarefs_metrics_groups_metas_qzas_dms_trees in betas.items():
        if not split:
            out_sh = '%s/run_nestedness_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
        nestedness_res[dat] = []
        for idx, metrics_groups_metas_qzas_dms_trees in enumerate(rarefs_metrics_groups_metas_qzas_dms_trees):
            nestedness_raref = {}
            cur_raref = datasets_rarefs[dat][idx]
            odir = get_analysis_folder(i_datasets_folder, 'nestedness/%s%s' % (dat, cur_raref))
            if split:
                out_sh = '%s/run_nestedness_%s_%s%s%s.sh' % (job_folder2, prjct_nm, dat, cur_raref, filt_raref)
            for _, groups_metas_qzas_dms_trees in metrics_groups_metas_qzas_dms_trees.items():
                for group, metas_qzas_mat_qzas_trees in groups_metas_qzas_dms_trees.items():
                    meta, qza, mat_qza, tree = metas_qzas_mat_qzas_trees[0]
                    meta_pd = read_meta_pd(meta).set_index('sample_name')
                    cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(subsets), 'nestedness')
                    for case_var, case_vals_list in cases_dict.items():
                        for case_vals in case_vals_list:
                            case = get_case(case_vals, case_var).replace(' ', '_')
                            cur_sh = '%s/run_nestedness_%s%s_%s_%s%s.sh' % (
                                job_folder2, dat, cur_raref, group, case, filt_raref)
                            cur_sh = cur_sh.replace(' ', '-')
                            all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                            res = run_single_nestedness(odir, group, meta_pd, nodfs, nulls,
                                                        modes, cur_sh, qza, case, case_var,
                                                        case_vals, binary, force)
                            nestedness_raref[(group, case)] = res
                break
            nestedness_res[dat].append(nestedness_raref)

    job_folder = get_job_folder(i_datasets_folder, 'nestedness')
    main_sh = write_main_sh(job_folder, '3_run_nestedness_%s%s' % (prjct_nm, filt_raref), all_sh_pbs,
                            '%s.prm%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        if p_nestedness_groups:
            print("# nestedness (config in %s)" % p_nestedness_groups)
        else:
            print("# nestedness")
        print_message('', 'sh', main_sh, jobs)

    return nestedness_res, colors


def run_single_nestedness(odir: str, group: str, meta_pd: pd.DataFrame, nodfs: list,
                          nulls: list, modes: list, cur_sh: str, qza: str, case: str,
                          case_var: str, case_vals: list, binary: str, force: bool) -> dict:
    res = {}
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        if group:
            cur_rad = '%s/%s_%s_%s' % (odir, splitext(basename(qza))[0], group, case)
        else:
            cur_rad = '%s/%s_%s' % (odir, splitext(basename(qza))[0], case)
        if not isdir(cur_rad):
            os.makedirs(cur_rad)

        new_meta = '%s.meta' % cur_rad
        new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
        cols = set()
        lat_lon_date = ['latitude', 'longitude', 'datetime']
        nodfs_valid = []
        for col in (nodfs + lat_lon_date):
            if col not in set(new_meta_pd.columns):
                continue
            if new_meta_pd[col].unique().size == 1:
                continue
            if col not in lat_lon_date and min(new_meta_pd[col].value_counts()) == 1:
                continue
            cols.add(col)
            if col in nodfs:
                nodfs_valid.append(col)
        new_meta_pd = new_meta_pd[sorted(cols)].reset_index()
        new_meta_pd.columns = (['#SampleID'] + sorted(cols))
        new_meta_pd.to_csv(new_meta, index=False, sep='\t')
        new_qza = '%s.qza' % cur_rad
        new_biom = '%s.biom' % cur_rad
        new_tsv = '%s.tsv' % cur_rad
        new_biom_meta = '%s_w-md.biom' % cur_rad

        if not isfile(new_biom):
            cmd = filter_feature_table(qza, new_qza, new_meta)
            cmd += run_export(new_qza, new_tsv, 'FeatureTable')
            cur_sh_o.write('echo "%s"\n' % cmd)
            cur_sh_o.write(cmd)

        cmd = run_add_metadata(new_biom, new_biom_meta, new_meta)
        cur_sh_o.write('echo "%s"\n' % cmd)
        cur_sh_o.write(cmd)

        graphs = '%s/graphs.csv' % cur_rad
        graphs_pdf = '%s/graphs.pdf' % cur_rad
        fields = '%s/fields.txt' % cur_rad
        res['graph'] = graphs
        res['graph_pdf'] = graphs_pdf
        res['fields'] = fields
        if not isfile(graphs) or not isfile(fields):
            write_nestedness_graph(
                new_biom_meta, cur_rad, graphs, binary,
                fields, nodfs_valid, cur_sh_o)
            remove = False

        for mode in modes:
            odir = '%s/%s' % (cur_rad, mode)
            if not isdir(odir):
                os.makedirs(odir)

            to_write = write_nestedness_nodfs(
                new_biom_meta, odir, binary,
                nodfs_valid, mode, nulls, cur_sh_o)
            if to_write:
                remove = False
            res.setdefault('modes', []).append(odir)
    if remove:
        os.remove(cur_sh)
    return res


# def get_comparisons_statistics_pd(mode_dir, level, cur_raref, group,
#                                   case, com_sta) -> pd.DataFrame:
#     com_sta_pds = []
#     for com_sta_fp in glob.glob('%s/*_%s.csv' % (mode_dir, com_sta)):
#         mode = mode_dir.split('/')[-1]
#         com_sta_pd = pd.read_csv(com_sta_fp)
#         if com_sta == 'simulate':
#             com_sta_pd['NULL'] = basename(com_sta_fp).split('_')[0]
#             if mode == 'overall':
#                 meta_name = np.nan
#             else:
#                 meta_name = basename(com_sta_fp).split('_simulate')[0].split('_', 1)[1]
#             com_sta_pd['METADATA'] = meta_name
#             com_sta_pd['MODE'] = mode
#             com_sta_pd['LEVEL'] = level
#             com_sta_pd['RAREF'] = cur_raref
#             com_sta_pd['FEATURE_SUBSET'] = group
#             com_sta_pd['SAMPLE_SUBSET'] = case
#             bins_lt = ['', '>', '>>', '>>>']
#             bins_gt = ['', '<', '<<', '<<<']
#             bins_vals = [0, 0.95, 0.97, 0.99]
#             ps = []
#             for (p, b) in [('PR_LT_OBSERVED', bins_lt), ('PR_GT_OBSERVED', bins_gt)]:
#                 ps.append(p)
#                 com_sta_pd[p] = [b[x-1] for x in np.digitize(com_sta_pd[p], bins=bins_vals)]
#             com_sta_pd['PVALUE'] = com_sta_pd[ps].apply(func=lambda x: ''.join(x), axis=1)
#             com_sta_pd = com_sta_pd.drop(columns=(['PR_ET_OBSERVED', 'NODF_SES'] + ps))
#         else:
#             com_sta_pd['COMPARISON'] = com_sta_pd.fillna('overall')[
#                 ['VERTEX_1_CLASSIFICATION', 'VERTEX_2_CLASSIFICATION']
#             ].apply(
#                 func=lambda x: '-'.join(list(set(x))), axis=1
#             )
#             com_sta_pd = com_sta_pd[['GRAPH_ID', 'COMPARISON']].drop_duplicates()
#         com_sta_pds.append(com_sta_pd)
#     if com_sta_pds:
#         com_sta_pd = pd.concat(com_sta_pds)
#         return com_sta_pd
#     return pd.DataFrame()


def nestedness_graphs(i_datasets_folder: str, nestedness_res: dict,
                      datasets: dict, split_taxa_pds: dict,
                      datasets_rarefs: dict, colors: dict,
                      datasets_collapsed_map: dict, collapsed: dict,
                      filt_raref: str, prjct_nm: str, qiime_env: str,
                      chmod: str, noloc: bool, split: bool, run_params: dict,
                      jobs: bool, chunkit: int) -> dict:

    RESOURCES = pkg_resources.resource_filename("routine_qiime2_analyses", "resources")
    nestedness_graphs_fp = '%s/nestedness_graphs.py' % RESOURCES

    job_folder2 = get_job_folder(i_datasets_folder, 'nestedness_figures/chunks')

    nodfs_fps = {}
    all_sh_pbs = {}
    for dat, nestedness_rarefs in nestedness_res.items():
        if not split:
            out_sh = '%s/run_nestedness_graphs_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
        stats_tax_dat = dat
        level = 'feature'
        if dat in datasets_collapsed_map:
            stats_tax_dat = datasets_collapsed_map[dat]
            level = ('%s_tx-' % stats_tax_dat).join(dat.split('%s_tx-' % stats_tax_dat)[1:])

        for idx, nestedness_raref in enumerate(nestedness_rarefs):
            for _, res in nestedness_raref.items():
                for mode_dir in res['modes']:
                    nodf_fpo = '%s/nodfs.tsv' % mode_dir
                    nodfs_fps.setdefault(stats_tax_dat, []).append(nodf_fpo)

            cur_raref = datasets_rarefs[dat][idx]
            if split:
                out_sh = '%s/run_nestedness_graphs_%s_%s%s%s.sh' % (job_folder2, prjct_nm, dat, cur_raref, filt_raref)
            out_py = out_sh.replace('.sh', '.py')

            cur_sh = '%s/run_nestedness_graphs_%s%s%s_tmp.sh' % (
                job_folder2, dat, cur_raref, filt_raref)
            cur_sh = cur_sh.replace(' ', '-')
            with open(cur_sh, 'w') as o:
                o.write('python3 %s\n' % out_py)
            all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)

            # value to edit in template
            tab_fp, meta_fp = datasets[dat][idx]
            if stats_tax_dat in split_taxa_pds:
                split_taxa_fp = split_taxa_pds[stats_tax_dat][1]
            else:
                split_taxa_fp = ''

            with open(out_py, 'w') as o, open(nestedness_graphs_fp) as f:
                for line in f:
                    line_edit = line
                    if '<DAT>' in line:
                        line_edit = line_edit.replace('<DAT>', dat)
                    if '<CUR_RAREF>' in line:
                        line_edit = line_edit.replace('<CUR_RAREF>', cur_raref)
                    if '<TAB_FP>' in line:
                        line_edit = line_edit.replace('<TAB_FP>', tab_fp)
                    if '<META_FP>' in line:
                        line_edit = line_edit.replace('<META_FP>', meta_fp)
                    if '<COLORS_SAMPLE>' in line:
                        line_edit = line_edit.replace('<COLORS_SAMPLE>', str(colors['sample']))
                    if '<COLORS_FEATURE>' in line:
                        line_edit = line_edit.replace('<COLORS_FEATURE>', str(colors['feature']))
                    if '<STATS_TAX_DAT>' in line:
                        line_edit = line_edit.replace('<STATS_TAX_DAT>', stats_tax_dat)
                    if '<SPLIT_TAXA_FP>' in line:
                        line_edit = line_edit.replace('<SPLIT_TAXA_FP>', split_taxa_fp)
                    if '<LEVEL>' in line:
                        line_edit = line_edit.replace('<LEVEL>', level)
                    if '<COLLAPSED>' in line:
                        line_edit = line_edit.replace('<COLLAPSED>', str(collapsed))
                    if '<NESTEDNESS_RAREF>' in line:
                        line_edit = line_edit.replace('<NESTEDNESS_RAREF>', str(nestedness_raref))
                    o.write(line_edit)

    job_folder = get_job_folder(i_datasets_folder, 'nestedness_figure')
    main_sh = write_main_sh(job_folder, 'run_nestedness_graphs_%s' % filt_raref, all_sh_pbs,
                            '%s.nstd.grph%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        print("# NESTEDNESS GRAPHS")
        print_message('', 'sh', main_sh, jobs)

    return nodfs_fps


def nestedness_nodfs(i_datasets_folder: str, nodfs_fps: dict,
                     collapsed: dict, filt_raref: str, prjct_nm: str,
                     qiime_env: str, chmod: str, noloc: bool, split: bool,
                     run_params: dict, jobs: bool, chunkit: int) -> None:

    RESOURCES = pkg_resources.resource_filename("routine_qiime2_analyses", "resources")
    nestedness_nodfs_fp = '%s/nestedness_nodfs.py' % RESOURCES

    job_folder2 = get_job_folder(i_datasets_folder, 'nestedness_figures/chunks')

    all_sh_pbs = {}
    for dat, nodfs in nodfs_fps.items():
        out_sh = '%s/run_nestedness_nodfs_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
        out_py = out_sh.replace('.sh', '.py')

        cur_sh = '%s/run_nestedness_nodfs_%s%s_tmp.sh' % (job_folder2, dat, filt_raref)
        cur_sh = cur_sh.replace(' ', '-')
        with open(cur_sh, 'w') as o:
            o.write('python3 %s\n' % out_py)
        all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)

        # value to edit in template
        odir = get_analysis_folder(i_datasets_folder, 'nestedness/%s%s' % (dat, filt_raref))
        with open(out_py, 'w') as o, open(nestedness_nodfs_fp) as f:
            for line in f:
                line_edit = line
                if '<DAT>' in line:
                    line_edit = line_edit.replace('<DAT>', dat)
                if '<ODIR>' in line:
                    line_edit = line_edit.replace('<ODIR>', odir)
                if '<NODFS>' in line:
                    line_edit = line_edit.replace('<NODFS>', str(nodfs))
                if '<COLLAPSED>' in line:
                    line_edit = line_edit.replace('<COLLAPSED>', str(collapsed))
                o.write(line_edit)

    job_folder = get_job_folder(i_datasets_folder, 'nestedness_figure')
    main_sh = write_main_sh(job_folder, 'run_nestedness_nodfs_%s' % filt_raref, all_sh_pbs,
                            '%s.nstd.ndf%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        print("# NESTEDNESS NODFS")
        print_message('', 'sh', main_sh, jobs)





# def nestedness_figure(i_datasets_folder: str, nestedness_res: dict,
#                       datasets_read: dict, split_taxa_pds: dict,
#                       datasets_rarefs: dict, colors: dict,
#                       datasets_collapsed_map: dict, collapsed: dict,
#                       filt_raref: str) -> None:
#
#     job_folder2 = get_job_folder(i_datasets_folder, 'nestedness_figures/chunks')
#
#     stats_tax = {}
#     for dat, nestedness_rarefs in nestedness_res.items():
#         stats_tax_dat = dat
#         level = 'feature'
#         print("datasets_collapsed_map")
#         print(datasets_collapsed_map)
#         print("dat")
#         print(dat)
#         if dat in datasets_collapsed_map:
#             stats_tax_dat = datasets_collapsed_map[dat]
#             print("stats_tax_dat")
#             print(stats_tax_dat)
#             level = ('%s_tx-' % stats_tax_dat).join(dat.split('%s_tx-' % stats_tax_dat)[1:])
#
#         nodfs_pds = []
#         for idx, nestedness_raref in enumerate(nestedness_rarefs):
#             cur_raref = datasets_rarefs[dat][idx]
#             tab_pd, meta_pd = datasets_read[dat][idx]
#             if colors['sample']:
#                 meta_pd = meta_pd.rename(columns={meta_pd.columns[0]: 'SAMPLE_ID'})
#                 meta_pd = meta_pd.set_index('SAMPLE_ID')
#                 if set(colors['sample']).issubset(meta_pd.columns):
#                     colors_in = list(set(colors['sample']) & set(meta_pd.columns))
#                     # meta_pd = meta_pd[colors_in]
#
#             tab_pd = np.log10(tab_pd + 1).stack().reset_index().rename(
#                 columns={'level_1': 'SAMPLE_ID', 0: 'log10_reads'})
#             tab_pd = tab_pd.rename(columns={tab_pd.columns[0]: 'OBSERVATION_ID'})
#
#             tax_cols = []
#             if colors['feature'] and stats_tax_dat in split_taxa_pds:
#                 tax_pd = split_taxa_pds[stats_tax_dat].reset_index()
#                 if level != 'feature':
#                     if level not in collapsed[stats_tax_dat]:
#                         continue
#                     tax_pd = split_taxa_pds[stats_tax_dat].iloc[
#                              :, :collapsed[stats_tax_dat][level]
#                         ].drop_duplicates()
#                     print("stats_tax_dat")
#                     print(stats_tax_dat)
#                     print("level")
#                     print(level)
#                     print(split_taxa_pds)
#                     print(dat, tax_pd)
#                     tax_pd['OBSERVATION_ID'] = tax_pd.apply(func=lambda x: ';'.join(x), axis=1)
#                 else:
#                     tax_pd = tax_pd.rename(columns={tax_pd.columns[0]: 'OBSERVATION_ID'})
#
#                 tax_pd = tax_pd.set_index('OBSERVATION_ID')
#                 tax_cols = [tax_pd.columns.tolist()[(x-1)] if isinstance(x, int)
#                                 and x not in tax_pd.columns and tax_pd.shape[1] >= x
#                             else x for x in colors['feature']]
#                 tax_cols = [x for x in tax_cols if x in tax_pd.columns]
#                 if tax_cols:
#                     tax_pd = tax_pd[tax_cols]
#
#             for (group, case), res in nestedness_raref.items():
#
#                 fields_fp = res['fields']
#                 graphs_fp = res['graph']
#                 if not isfile(fields_fp) or not isfile(graphs_fp):
#                     continue
#
#                 graphs = pd.read_csv(graphs_fp, header=0, sep=',', dtype={'SAMPLE_ID': str})
#                 samples_order = [y for x, y in sorted(dict(graphs[['SAMPLE_RANK', 'SAMPLE_ID']].values).items())][::-1]
#                 features_order = [y for x, y in sorted(dict(graphs[['OBSERVATION_RANK', 'OBSERVATION_ID']].values).items())][::-1]
#                 graphs = graphs.merge(tab_pd, on=['SAMPLE_ID', 'OBSERVATION_ID'], how='left')
#
#                 matrix = graphs[['OBSERVATION_ID', 'SAMPLE_ID', 'log10_reads']].pivot_table(
#                     values='log10_reads', index='OBSERVATION_ID', columns='SAMPLE_ID')
#                 matrix = matrix.loc[features_order, samples_order]
#
#                 fields = [x.strip() for x in open(fields_fp).readlines()]
#                 new_meta = 'metadata'
#                 cur_meta_pd = meta_pd.loc[matrix.columns.tolist(), fields].copy()
#
#                 if new_meta in cur_meta_pd.columns:
#                     new_meta = 'passed_metadata'
#                 if len(fields) > 1:
#                     cur_meta_pd[new_meta] = cur_meta_pd[fields].apply(func=lambda x: '-'.join(x), axis=1)
#                 cur_meta_leg_hex = cur_meta_pd.apply(func=lambda x: dict(
#                     zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size).as_hex())
#                 )).to_dict()
#                 cur_meta_leg_rgb = cur_meta_pd.apply(func=lambda x: dict(
#                     zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size))
#                 )).to_dict()
#                 cur_meta_pd_hex = cur_meta_pd.apply(func=lambda x: x.map(dict(
#                     zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size).as_hex())
#                 )))
#                 cur_meta_pd_rgb = cur_meta_pd.apply(func=lambda x: x.map(dict(
#                     zip(x.unique(), sns.color_palette(palette='Set1', n_colors=x.unique().size))
#                 )))
#
#                 if tax_cols:
#                     cur_tax_pd = tax_pd.loc[matrix.index.tolist(), :].copy()
#                     cur_tax_leg_hex = cur_tax_pd.apply(func=lambda x: dict(
#                         zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size).as_hex())
#                     )).to_dict()
#                     cur_tax_leg_rgb = cur_tax_pd.apply(func=lambda x: dict(
#                         zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size))
#                     )).to_dict()
#                     cur_tax_pd_hex = cur_tax_pd.apply(func=lambda x: x.map(dict(
#                         zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size).as_hex())
#                     )))
#                     cur_tax_pd_rgb = cur_tax_pd.apply(func=lambda x: x.map(dict(
#                         zip(x.unique(), sns.color_palette(palette='Paired', n_colors=x.unique().size))
#                     )))
#
#                 X, Y = matrix.shape
#                 graphs_pdf = res['graph_pdf']
#                 graphs_pdf_complex = '%s_complex.pdf' % splitext(graphs_pdf)[0]
#                 # if not isfile(graphs_pdf_complex):
#                 if 1:
#                     with PdfPages(graphs_pdf_complex) as pdf:
#                         fig, ax = plt.subplots(figsize=(20, (20*Y/X)))
#                         if tax_cols:
#                             g = sns.clustermap(
#                                 matrix, col_cluster=False, row_cluster=False,
#                                 linewidths=0.1, cmap='coolwarm',
#                                 row_colors=cur_tax_pd_hex, col_colors=cur_meta_pd_hex,
#                                 yticklabels=False, xticklabels=False
#                             )
#                             simples = [('feature', cur_tax_pd_rgb, cur_tax_leg_rgb),
#                                        ('sample', cur_meta_pd_rgb, cur_meta_leg_rgb)]
#                         else:
#                             g = sns.clustermap(
#                                 matrix, col_cluster=False, row_cluster=False,
#                                 linewidths=0.1, cmap='coolwarm',
#                                 col_colors=cur_meta_pd_hex,
#                                 yticklabels=False, xticklabels=False
#                             )
#                             simples = [('sample', cur_meta_pd_rgb, cur_meta_leg_rgb)]
#
#                         n = 0
#                         N = 1 / cur_meta_pd_hex.columns.size
#                         for cdx, col in enumerate(cur_meta_pd_hex.columns):
#                             n += N
#                             if not cdx:
#                                 first_leg = g.ax_col_dendrogram.legend(handles=[
#                                     mpatches.Patch(label=x, color=y) for x, y in cur_meta_leg_hex[col].items()
#                                 ], loc="upper left")
#                                 g.ax_col_dendrogram.add_artist(first_leg)
#                             else:
#                                 g.ax_col_dendrogram.legend(handles=[
#                                     mpatches.Patch(label=x, color=y) for x, y in cur_meta_leg_hex[col].items()
#                                 ], loc="upper left", bbox_to_anchor=(n, 1))
#
#                         if tax_cols:
#                             n = 0
#                             N = 1 / cur_tax_pd_hex.columns.size
#                             for cdx, col in enumerate(cur_tax_pd_hex.columns):
#                                 n += N
#                                 if not cdx:
#                                     first_leg = g.ax_row_dendrogram.legend(handles=[
#                                         mpatches.Patch(label=x, color=y) for x, y in cur_tax_leg_hex[col].items()
#                                     ], loc="lower left")
#                                     g.ax_row_dendrogram.add_artist(first_leg)
#                                 else:
#                                     g.ax_row_dendrogram.legend(handles=[
#                                         mpatches.Patch(label=x, color=y) for x, y in cur_tax_leg_hex[col].items()
#                                     ], loc="lower left", bbox_to_anchor=(0, n))
#
#                         g.ax_heatmap.set_xlabel('Samples (sorted by richness)')
#                         g.ax_heatmap.set_ylabel('Taxon (sorted by prevalence)')
#
#                         suptitle = '%s%s' % (dat, cur_raref)
#                         if case != 'ALL':
#                             suptitle = suptitle + '\nsamples subset: %s' % case
#                         if group != '':
#                             suptitle = suptitle + '\nfeatures subset: %s' % group
#                         plt.suptitle(suptitle, fontsize=20)
#                         plt.subplots_adjust(top=.95, hspace=0.3)
#                         pdf.savefig(bbox_inches='tight', dpi=300)
#                         plt.close()
#
#                 graphs_pdf_simples = '%s_simples.pdf' % splitext(graphs_pdf)[0]
#                 # if not isfile(graphs_pdf_simples):
#                 if 1:
#                     with PdfPages(graphs_pdf_simples) as pdf:
#                         for (typ, tab, leg) in simples:
#                             num = tab.columns.size
#                             fig, axes = plt.subplots(num, 1, figsize=(9, (tab.columns.size*4)))
#                             for cdx, col in enumerate(tab.columns):
#                                 leg_col2rgb = dict((x, leg[col][x]) for idx, x in enumerate(leg[col]))
#                                 cols_d = tab[col].to_dict()
#                                 if typ == 'sample':
#                                     mat = [
#                                         [
#                                             [np.nan, np.nan, np.nan, 0.] if str(x[1]) == 'nan' else
#                                             ([y for y in cols_d[x[0]]] + [1.]) for x in row.reset_index().values
#                                         ] for r, row in matrix.iterrows()
#                                     ]
#                                 if typ == 'feature':
#                                     mat = [
#                                         [
#                                             [np.nan, np.nan, np.nan, 0.] if str(x) == 'nan' else
#                                             ([y for y in cols_d[r]] + [1.]) for x in row
#                                         ] for r, row in matrix.iterrows()
#                                     ]
#                                 mat = np.array(mat)
#                                 hands = [mpatches.Patch(label=x, color=y) for x, y in leg_col2rgb.items()]
#                                 if num == 1:
#                                     axes.imshow(mat, aspect="auto")
#                                     axes.legend(handles=hands, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#                                     axes.set_title(col)
#                                     axes.set_xlabel('Samples (sorted by richness)')
#                                     axes.set_ylabel('Phyla (sorted by prevalence)')
#                                 else:
#                                     axes[cdx].imshow(mat, aspect="auto")
#                                     axes[cdx].legend(handles=hands, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#                                     axes[cdx].set_title(col)
#                                     axes[cdx].set_xlabel('Samples (sorted by richness)')
#                                     axes[cdx].set_ylabel('Phyla (sorted by prevalence)')
#                             suptitle = '%s%s (%s coloring)' % (dat, cur_raref, typ)
#                             top = .93
#                             if case != 'ALL':
#                                 top -= 0.03
#                                 suptitle = suptitle + '\nsamples subset: %s' % case
#                             if group != '':
#                                 top -= 0.03
#                                 suptitle = suptitle + '\nfeatures subset: %s' % group
#                             plt.suptitle(suptitle, fontsize=15)
#                             plt.subplots_adjust(top=top, hspace=0.3)
#                             pdf.savefig(bbox_inches='tight', dpi=300)
#                             plt.close()
#
#                 for mode_dir in res['modes']:
#                     comparisons_pd = get_comparisons_statistics_pd(
#                         mode_dir, level, cur_raref, group, case, 'comparisons')
#                     if not comparisons_pd.shape[0]:
#                         continue
#                     # statistics_pd = get_comparisons_statistics_pd(
#                     #     mode_dir, level, raref, group, case, 'statistics')
#                     simulate_pd = get_comparisons_statistics_pd(
#                         mode_dir, level, cur_raref, group, case, 'simulate')
#                     if not simulate_pd.shape[0]:
#                         continue
#                     # nodf_pd = statistics_pd[statistics_pd['GRAPH_EDGE_COUNT'] > 10].merge(
#                     #     comparisons_pd, on='GRAPH_ID', how='left')
#                     nodf_pd = simulate_pd[simulate_pd['GRAPH_EDGE_COUNT'] > 5].merge(
#                         comparisons_pd, on='GRAPH_ID', how='left')
#                     nodfs_pds.append(nodf_pd)
#         if nodfs_pds:
#             nodfs_pd = pd.concat(nodfs_pds)
#             stats_tax.setdefault(stats_tax_dat, []).append(nodfs_pd)
#
#     if stats_tax:
#         print('Writing nodfs pdfs...', end=' ')
#         for dat, nodfs_pds in stats_tax.items():
#             if dat not in collapsed:
#                 continue
#             levels = collapsed[dat]
#             levels['feature'] = max(levels.values())+1
#             nodfs = pd.concat(nodfs_pds).drop(columns=['GRAPH_ID'])
#             nodfs_cols = ['NODF_OBSERVED', 'NODF_NULL_MEAN']
#             nodfs = nodfs.set_index([x for x in nodfs.columns if x not in nodfs_cols])
#             nodfs = nodfs.stack().reset_index().rename(columns={'level_11': 'NODF_TYPE', 0: 'NODF'})
#             nodfs_obs = nodfs.loc[nodfs.NODF_TYPE == 'NODF_OBSERVED'].drop(columns=['NODF_NULL_STDEV'])
#             nodfs_null = nodfs.loc[nodfs.NODF_TYPE == 'NODF_NULL_MEAN'].copy()
#
#             nodfs_null_up = nodfs_null.copy()
#             nodfs_null_down = nodfs_null.copy()
#             nodfs_null_up["NODF"] = nodfs_null_up.NODF + nodfs_null_up.NODF_NULL_STDEV
#             nodfs_null_down["NODF"] = nodfs_null_up.NODF - nodfs_null_up.NODF_NULL_STDEV
#             nodfs_null = pd.concat([nodfs_obs, nodfs_null_up, nodfs_null_down])
#
#             nodfs_metas = nodfs_null.loc[nodfs_null.COMPARISON != 'overall']
#             nodfs_leg = dict(zip(nodfs_metas.COMPARISON.unique(), sns.color_palette(
#                         palette='Set1', n_colors=nodfs_metas.COMPARISON.unique().size).as_hex()))
#             nodfs_leg['overall'] = '#000000'
#
#             nodfs_null['SUBSET'] = nodfs_null[
#                 ['FEATURE_SUBSET', 'SAMPLE_SUBSET']].apply(func=lambda x: ''.join(x), axis=1)
#             comparison_raref = 'COMPARISON_RAREF'
#             nodfs_null[comparison_raref] = nodfs_null[
#                 ['COMPARISON', 'RAREF']].apply(func=lambda x: ''.join(x), axis=1)
#             if not sum(nodfs_null.COMPARISON_RAREF.str.contains('raref')):
#                 comparison_raref = 'COMPARISON'
#
#             nodfs_null['MODE_METADATA'] = nodfs_null[
#                 ['MODE', 'METADATA']].fillna('').apply(func=lambda x: ' - '.join([y for y in x if y]), axis=1)
#             nodfs_null['LEVEL_SORT'] = [levels[x] for x in nodfs_null['LEVEL']]
#             nodfs_null = nodfs_null.sort_values(['LEVEL_SORT', 'NODF_TYPE'], ascending=True)
#
#             nodfs_null = nodfs_null.drop(columns=[
#                 'FEATURE_SUBSET', 'SAMPLE_SUBSET', 'RAREF', 'MODE', 'METADATA'
#             ])
#
#             odir = get_analysis_folder(i_datasets_folder, 'nestedness/%s%s' % (dat, filt_raref))
#             table_fpo = '%s/nodfs.tsv' % odir
#             nodfs_null.to_csv(table_fpo, index=False, sep='\t')
#             nodfs_pdf = '%s/nodfs.pdf' % odir
#             with PdfPages(nodfs_pdf) as pdf:
#                 for (mode_metadata), nodfs_null_mode_ in nodfs_null.groupby('MODE_METADATA'):
#                     nodfs_null_mode = nodfs_null_mode_.drop(columns=['MODE_METADATA'])
#                     SUBSETS = ['ALL'] + sorted([x for x in nodfs_null_mode.SUBSET.unique().tolist() if x != 'ALL'])
#                     g = sns.relplot(
#                         data=nodfs_null_mode, x="LEVEL", y="NODF",
#                         col="SUBSET", col_wrap=2, col_order=SUBSETS,
#                         hue=comparison_raref, palette=nodfs_leg,
#                         style='NODF_TYPE', style_order=['NODF_OBSERVED', 'NODF_NULL_MEAN'],
#                         kind="line", err_style="bars", ci='sd')
#
#                     axes = g.axes
#                     X = 0
#                     for SUBSET in SUBSETS:
#                         SUBSET_pd = nodfs_null_mode.loc[
#                             (nodfs_null_mode.SUBSET == SUBSET) &
#                             (nodfs_null_mode.NODF_TYPE == 'NODF_OBSERVED')
#                         ].drop(columns=['SUBSET'])
#
#                         for (LEVEL_SORT, NODF, PVALUE, COMP) in SUBSET_pd[
#                             ['LEVEL_SORT', 'NODF', 'PVALUE', comparison_raref]].values:
#                             if PVALUE:
#                                 axes[X].text(
#                                     (LEVEL_SORT-1), NODF+0.005, PVALUE,
#                                     color=nodfs_leg[COMP], fontsize=5
#                                 )
#                         X += 1
#
#                     suptitle = '%s - %s' % (dat, mode_metadata)
#                     plt.suptitle(suptitle, fontsize=12)
#                     plt.subplots_adjust(top=.90, hspace=0.3)
#                     pdf.savefig(bbox_inches='tight', dpi=300)
#                     plt.close()
#
#         print('Done!')

