# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import yaml
import pandas as pd
from os.path import basename, isfile, isdir, splitext
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_job_folder, get_analysis_folder, run_import


def run_permanova(i_folder: str, datasets: dict, betas: dict,
                  beta_metrics: list, testing_groups: tuple, p_perm_groups: str,
                  force: bool, prjct_nm: str, qiime_env: str):

    print("# PERMANOVA (groups config in %s)" % p_perm_groups)
    def run_multi_perm(
            job_folder2, odir, dat,
            mat, mat_pd,
            tsv, tsv_pd,
            meta_pd,
            case_, case_var, case_vals,
            testing_groups):

        written = 0
        out_sh = '%s/run_beta_group_significance_%s.sh' % (job_folder2, case_)
        out_pbs = out_sh.replace('.sh', '.pbs')
        with open(out_sh, 'w') as sh:

            for testing_group in testing_groups:

                if testing_group == 'ALL':
                    continue
                else:
                    case = '%s__%s' % (case_, testing_group)

                cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s' % case)
                new_tsv = '%s.tsv' % cur_rad
                new_meta = new_tsv.replace('/tab_', '/meta_')
                new_qza = '%s.qza' % cur_rad
                new_qzv = '%s_permanova.qzv' % cur_rad

                if force or not isfile(new_qzv):
                    new_mat = '%s/%s.tsv' % (odir, basename(mat).replace('.tsv', '_%s' % case))
                    new_mat_qza = new_mat.replace('.tsv', '.qza')

                    if 'ALL' in case:
                        new_meta_pd = meta_pd.copy()
                    elif len([x for x in case_vals if '>' in x or '<' in x]):
                        new_meta_pd = meta_pd.copy()
                        for case_val in case_vals:
                            if case_val[0] == '>':
                                new_meta_pd = new_meta_pd[new_meta_pd[case_var] >= float(case_val[1:])].copy()
                            elif case_val[0] == '<':
                                new_meta_pd = new_meta_pd[new_meta_pd[case_var] <= float(case_val[1:])].copy()
                    else:
                        new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()

                    if new_meta_pd[testing_group].unique().size > 1:
                        q2types = pd.DataFrame(
                            [(['#q2:types'] + ['categorical'] * new_meta_pd.shape[1])],
                            columns=new_meta_pd.reset_index().columns.tolist(),
                        )
                        q2types = q2types.set_index('#SampleID')
                        new_meta_pd = pd.concat([q2types, new_meta_pd])
                        new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')

                        new_mat_pd = mat_pd.loc[
                            new_meta_pd.index.tolist()[1:],
                            new_meta_pd.index.tolist()[1:]
                        ].copy()
                        new_mat_pd.to_csv(new_mat, index=True, sep='\t')

                        new_tsv_pd = tsv_pd[new_meta_pd.index.tolist()[1:]].copy()
                        new_tsv_pd = new_tsv_pd.loc[new_tsv_pd.sum(1) > 0, :]
                        new_tsv_pd = new_tsv_pd.loc[:, new_tsv_pd.sum(0) > 0]
                        new_tsv_pd.reset_index().to_csv(new_tsv, index=False, sep='\t')

                        cmd = run_import(new_tsv, new_qza, "FeatureTable[Frequency]")
                        sh.write(cmd)
                        cmd = run_import(new_mat, new_mat_qza, "DistanceMatrix")
                        sh.write(cmd)

                        cmd = [
                            'qiime', 'diversity', 'beta-group-significance',
                            '--i-distance-matrix', new_mat_qza,
                            '--m-metadata-file', new_meta,
                            '--m-metadata-column', testing_group,
                            '--p-permutations', '2999',
                            '--o-visualization', new_qzv]
                        sh.write('echo "%s"\n' % ' '.join(cmd))
                        sh.write('%s\n' % ' '.join(cmd))
                        written += 1
        if written:
            xpbs_call(out_sh, out_pbs, '%s.perm.%s' % (prjct_nm, case_), qiime_env,
                      '2', '1', '1', '2', 'gb')
        else:
            os.remove(out_sh)

    job_folder = get_job_folder(i_folder, 'permanova')
    job_folder2 = get_job_folder(i_folder, 'permanova/chunks')

    with open(p_perm_groups) as handle:
        cases_dict = yaml.load(handle)
        # cases_dict = yaml.load(handle, Loader=yaml.FullLoader)
    cases_dict.update({'ALL': [[]]})

    jobs = []
    first_print = 0
    main_sh = '%s/3_run_beta_group_significance.sh' % job_folder
    with open(main_sh, 'w') as o:

        for dat, tsvs_metas_list in datasets.items():
            odir = get_analysis_folder(i_folder, 'permanova/%s' % dat)
            for tsv_meta in tsvs_metas_list:
                tsv, meta = tsv_meta
                tsv_pd = pd.read_csv(tsv, header=0, index_col=0, sep='\t')
                meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')

                for mat_qza in betas[dat][meta]:
                    mat = mat_qza.replace('.qza', '.tsv')
                    for metric in beta_metrics:
                        if metric in mat:
                            break

                    if not isfile(mat):
                        if not first_print:
                            print('Beta diversity, distances matrices must be generated already to automatise PERMANOVA\n'
                                  '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)')
                            first_print += 1
                        continue
                    mat_pd = pd.read_csv(mat, header=0, index_col=0, sep='\t')

                    for case_var, case_vals_list in cases_dict.items():
                        testing_groups = testing_groups + [case_var]
                        for case_vals in case_vals_list:
                            if len(case_vals):
                                case = '%s_%s_%s' % (metric, case_var, '-'.join(
                                    [x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
                            else:
                                case = '%s_%s' % (metric, case_var)
                            out_pbs = '%s/run_beta_group_significance_%s.pbs' % (job_folder2, case)
                            o.write('qsub %s\n' % out_pbs)
                            p = multiprocessing.Process(
                                target=run_multi_perm,
                                args=(
                                    job_folder2, odir, dat,
                                    mat, mat_pd,
                                    tsv, tsv_pd,
                                    meta_pd,
                                    case, case_var, case_vals,
                                    testing_groups
                                )
                            )
                            p.start()
                            jobs.append(p)
    for j in jobs:
        j.join()

    print('[TO RUN] sh', main_sh)


def run_adonis(p_formulas: str, i_folder: str, datasets: dict, betas: dict,
               beta_metrics: list, p_perm_groups: str,
               force: bool, prjct_nm: str, qiime_env: str):

    print("# Run Adonis (groups config in %s)" % p_perm_groups)
    def run_multi_adonis(formula, odir, dat, mat, mat_pd, tsv, tsv_pd,
                         meta_pd, case, case_var, case_vals, out_sh):

        with open(out_sh, 'w') as cur_sh:
            cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s' % case)
            new_tsv = '%s.tsv' % cur_rad
            new_meta = new_tsv.replace('/tab_', '/meta_')
            new_qza = '%s.qza' % cur_rad
            new_qzv = '%s_adonis.qzv' % cur_rad
            print(new_qzv)

            if force or not isfile(new_qzv):
                new_mat = '%s/%s.tsv' % (odir, basename(mat).replace('.tsv', '_%s' % case))
                new_mat_qza = new_mat.replace('.tsv', '.qza')

                if 'ALL' in case:
                    new_meta_pd = meta_pd.copy()
                elif len([x for x in case_vals if '>' in x or '<' in x]):
                    new_meta_pd = meta_pd.copy()
                    for case_val in case_vals:
                        if case_val[0] == '>':
                            new_meta_pd = new_meta_pd[new_meta_pd[case_var] >= float(case_val[1:])].copy()
                        elif case_val[0] == '<':
                            new_meta_pd = new_meta_pd[new_meta_pd[case_var] <= float(case_val[1:])].copy()
                else:
                    new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()
                new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')
                new_mat_pd = mat_pd.loc[
                    new_meta_pd.index.tolist()[1:],
                    new_meta_pd.index.tolist()[1:]
                ].copy()
                new_mat_pd.to_csv(new_mat, index=True, sep='\t')

                new_tsv_pd = tsv_pd[new_meta_pd.index.tolist()[1:]].copy()
                new_tsv_pd = new_tsv_pd.loc[new_tsv_pd.sum(1) > 0, :]
                new_tsv_pd = new_tsv_pd.loc[:, new_tsv_pd.sum(0) > 0]
                new_tsv_pd.reset_index().to_csv(new_tsv, index=False, sep='\t')

                cmd = run_import(new_tsv, new_qza, "FeatureTable[Frequency]")
                cur_sh.write(cmd)
                cmd = run_import(new_mat, new_mat_qza, "DistanceMatrix")
                cur_sh.write(cmd)
                cmd = ['qiime', 'diversity', 'adonis',
                       '--i-distance-matrix', new_mat_qza,
                       '--m-metadata-file', new_meta,
                       '--p-formula', "%s" % formula,
                       '--p-permutations', '2999',
                       '--p-n-jobs', '6',
                       '--o-visualization', new_qzv]
                cur_sh.write('echo "%s"\n' % ' '.join(cmd))
                cur_sh.write('%s\n' % ' '.join(cmd))

    job_folder = get_job_folder(i_folder, 'adonis')
    job_folder2 = get_job_folder(i_folder, 'adonis/chunks')

    with open(p_perm_groups) as handle:
        cases_dict = yaml.load(handle)
        # cases_dict = yaml.load(handle, Loader=yaml.FullLoader)
    cases_dict.update({'ALL': [[]]})

    with open(p_formulas) as handle:
        formulas = yaml.load(handle)
        # formulas = yaml.load(handle, Loader=yaml.FullLoader)

    jobs = []
    first_print = 0
    main_sh = '%s/3_run_adonis.sh' % job_folder
    with open(main_sh, 'w') as main_o:

        for dat, tsvs_metas_list in datasets.items():
            odir = get_analysis_folder(i_folder, 'adonis/%s' % dat)
            out_sh = '%s/run_adonis_%s.sh' % (job_folder2, dat)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                for tsv_meta in tsvs_metas_list:
                    tsv, meta = tsv_meta
                    tsv_pd = pd.read_csv(tsv, header=0, index_col=0, sep='\t')
                    meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')

                    for mat_qza in betas[dat][meta]:
                        mat = mat_qza.replace('.qza', '.tsv')
                        for metric in beta_metrics:
                            if metric in mat:
                                break
                        if not isfile(mat):
                            if not first_print:
                                print('Beta diversity, distances matrices must be generated already to automatise adonis\n'
                                      '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)')
                                first_print += 1
                            continue
                        mat_pd = pd.read_csv(mat, header=0, index_col=0, sep='\t')
                        for form, formula in formulas.items():
                            for case_var, case_vals_list in cases_dict.items():
                                for case_vals in case_vals_list:
                                    if len(case_vals):
                                        case = '%s_%s_%s_%s' % (metric, case_var, '-'.join(
                                            [x.replace('<', 'below').replace('>', 'above') for x in case_vals]), form)
                                    else:
                                        case = '%s_%s_%s' % (metric, case_var, form)
                                    out_sh = '%s/run_adonis_%s_%s.sh' % (job_folder2, dat, case)
                                    sh.write('sh %s\n' % out_sh)
                                    print(form)
                                    print(formula)
                                    print(case_var)
                                    print(case_vals)
                                    print(case)
                                    print(out_sh)
                                    print(out_shd)
                                    p = multiprocessing.Process(
                                        target=run_multi_adonis,
                                        args=(formula, odir, dat, mat, mat_pd, tsv, tsv_pd,
                                              meta_pd, case, case_var, case_vals, out_sh)
                                    )
                                    p.start()
                                    jobs.append(p)
            xpbs_call(out_sh, out_pbs, '%s.dns.%s' % (prjct_nm, dat), qiime_env,
                      '2', '1', '6', '2', 'gb')
            main_o.write('qsub %s\n' % out_pbs)

    for j in jobs:
        j.join()

    print('[TO RUN] sh', main_sh)
