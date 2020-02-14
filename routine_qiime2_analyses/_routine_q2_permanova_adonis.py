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
from os.path import basename, isfile, splitext
import multiprocessing

from routine_qiime2_analyses._routine_q2_xpbs import xpbs_call
from routine_qiime2_analyses._routine_q2_io_utils import get_metrics, get_job_folder, get_analysis_folder


def run_permanova(i_folder: str, datasets: dict, betas: dict,
                  testing_groups: tuple, p_perm_groups: str,
                  force: bool, prjct_nm: str, qiime_env: str):

    beta_metrics = get_metrics('beta_metrics')

    def run_multi_perm(out_sh, out_pbs, odir, mat_qza, tsv, meta_pd,
                       case_, case_var, case_vals, testing_groups):

        qza = '%s.qza' % splitext(tsv)[0]
        written = 0
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
                # new_mat = '%s/%s.tsv' % (odir, basename(mat).replace('.tsv', '_%s' % case))
                # new_mat_qza = new_mat.replace('.tsv', '.qza')
                new_mat_qza = '%s/%s.tsv' % (odir, basename(mat_qza).replace('.qza', '_%s.qza' % case))

                if force or not isfile(new_qzv):

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
                        q2types.rename(columns={q2types.columns.tolist()[0]: new_meta_pd.index.name}, inplace=True)
                        q2types.set_index(new_meta_pd.index.name, inplace=True)
                        new_meta_pd = pd.concat([q2types, new_meta_pd])
                        new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')

                        cmd = 'qiime diversity filter-distance-matrix \\ \n'
                        cmd += '--m-metadata-file %s \\ \n' % new_meta
                        cmd += '--i-distance-matrix %s \\ \n' % mat_qza
                        cmd += '--o-filtered-distance-matrix %s\n' % new_mat_qza
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n' % cmd)

                        cmd = 'qiime feature-table filter-samples \\ \n'
                        cmd += '--i-table %s \\ \n' % qza
                        cmd += '--m-metadata-file %s \\ \n' % new_meta
                        cmd += '--o-filtered-table %s\n' % new_qza
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n' % cmd)

                        cmd = 'qiime diversity beta-group-significance \\ \n'
                        cmd += '--i-distance-matrix %s \\ \n' % new_mat_qza
                        cmd += '--m-metadata-file %s \\ \n' % new_meta
                        cmd += '--m-metadata-column %s \\ \n' % testing_group
                        cmd += '--p-permutations 2999 \\ \n'
                        cmd += '--o-visualization %s\n' % new_qzv
                        sh.write('echo "%s"\n' % cmd)
                        sh.write('%s\n' % cmd)
                        written += 1
        if written:
            xpbs_call(out_sh, out_pbs, '%s.perm.%s' % (prjct_nm, case_), qiime_env,
                      '2', '1', '1', '2', 'gb')
        else:
            os.remove(out_sh)

    job_folder = get_job_folder(i_folder, 'permanova')
    job_folder2 = get_job_folder(i_folder, 'permanova/chunks')

    cases_dict = {'ALL': [[]]}
    if p_perm_groups:
        with open(p_perm_groups) as handle:
            # cases_dict = yaml.load(handle)
            cases_dict.update(yaml.load(handle, Loader=yaml.FullLoader))

    jobs = []
    all_shs = []
    first_print = 0
    main_sh = '%s/3_run_beta_group_significance.sh' % job_folder
    with open(main_sh, 'w') as o:

        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            odir = get_analysis_folder(i_folder, 'permanova/%s' % dat)
            meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')

            for mat_qza in betas[dat][meta]:
                for metric in beta_metrics:
                    if metric in mat_qza:
                        break

                if not isfile(mat_qza):
                    if not first_print:
                        print('Beta diversity, distances matrices must be generated already to automatise PERMANOVA\n'
                              '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)')
                        first_print += 1
                    continue

                for case_var, case_vals_list in cases_dict.items():
                    testing_groups = list(testing_groups) + [case_var]
                    for case_vals in case_vals_list:
                        if len(case_vals):
                            case = '%s_%s_%s' % (metric, case_var, '-'.join(
                                [x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
                        else:
                            case = '%s_%s' % (metric, case_var)
                        out_sh = '%s/run_beta_group_significance_%s.sh' % (job_folder2, case)
                        out_pbs = '%s.pbs' % splitext(out_sh)[0]
                        p = multiprocessing.Process(
                            target=run_multi_perm,
                            args=(out_sh, out_pbs, odir, mat_qza, tsv, meta_pd,
                                  case, case_var, case_vals, testing_groups)
                        )
                        p.start()
                        jobs.append(p)
                        if isfile(out_sh):
                            all_shs.append(out_sh)
                            o.write('qsub %s\n' % out_pbs)
    for j in jobs:
        j.join()

    if len([1 for sh in all_shs for line in open(sh).readlines() if len(line.strip())]):
        if p_perm_groups:
            print("# PERMANOVA (groups config in %s)" % p_perm_groups)
        else:
            print("# PERMANOVA")
        print('[TO RUN] sh', main_sh)


def run_adonis(p_formulas: str, i_folder: str, datasets: dict, betas: dict,
               p_perm_groups: str, force: bool, prjct_nm: str, qiime_env: str):

    beta_metrics = get_metrics('beta_metrics')

    def run_multi_adonis(formula, odir, mat_qza, tsv,
                         meta_pd, case, case_var, case_vals, out_sh):

        qza = '%s.qza' % splitext(tsv)[0]
        with open(out_sh, 'w') as cur_sh:
            cur_rad = odir + '/' + basename(tsv).replace('.tsv', '_%s' % case)
            new_tsv = '%s.tsv' % cur_rad
            new_meta = new_tsv.replace('/tab_', '/meta_')
            new_qza = '%s.qza' % cur_rad
            new_qzv = '%s_adonis.qzv' % cur_rad
            new_mat_qza = '%s/%s.tsv' % (odir, basename(mat_qza).replace('.qza', '_%s.qza' % case))

            if force or not isfile(new_qzv):

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

                cmd = 'qiime diversity filter-distance-matrix \\ \n'
                cmd += '--m-metadata-file %s \\ \n' % new_meta
                cmd += '--i-distance-matrix %s \\ \n' % mat_qza
                cmd += '--o-filtered-distance-matrix %s\n' % new_mat_qza
                cur_sh.write('echo "%s"\n' % cmd)
                cur_sh.write('%s\n' % cmd)

                cmd = 'qiime feature-table filter-samples \\ \n'
                cmd += '--i-table %s \\ \n' % qza
                cmd += '--m-metadata-file %s \\ \n' % new_meta
                cmd += '--o-filtered-table %s\n' % new_qza
                cur_sh.write('echo "%s"\n' % cmd)
                cur_sh.write('%s\n' % cmd)

                cmd = 'qiime diversity adonis \\ \n'
                cmd += '--i-distance-matrix %s \\ \n' % new_mat_qza
                cmd += '--m-metadata-file %s \\ \n' % new_meta
                cmd += '--p-formula "%s" \\ \n' % formula
                cmd += '--p-permutations 2999 \\ \n'
                cmd += '--p-n-jobs 6 \\ \n'
                cmd += '--o-visualization %s\n' % new_qzv
                cur_sh.write('echo "%s"\n' % cmd)
                cur_sh.write('%s\n' % cmd)

    job_folder = get_job_folder(i_folder, 'adonis')
    job_folder2 = get_job_folder(i_folder, 'adonis/chunks')

    cases_dict = {'ALL': [[]]}
    if p_perm_groups:
        with open(p_perm_groups) as handle:
            # cases_dict = yaml.load(handle)
            cases_dict.update(yaml.load(handle, Loader=yaml.FullLoader))

    with open(p_formulas) as handle:
        # formulas = yaml.load(handle)
        formulas = yaml.load(handle, Loader=yaml.FullLoader)

    jobs = []
    all_shs = []
    first_print = 0
    main_sh = '%s/3_run_adonis.sh' % job_folder
    with open(main_sh, 'w') as main_o:

        for dat, tsv_meta_pds in datasets.items():
            tsv, meta = tsv_meta_pds
            odir = get_analysis_folder(i_folder, 'adonis/%s' % dat)
            out_sh = '%s/run_adonis_%s.sh' % (job_folder2, dat)
            all_shs.append(out_sh)
            out_pbs = '%s.pbs' % splitext(out_sh)[0]
            with open(out_sh, 'w') as sh:
                # tsv_pd = pd.read_csv(tsv, header=0, index_col=0, sep='\t')
                meta_pd = pd.read_csv(meta, header=0, index_col=0, sep='\t')

                for mat_qza in betas[dat][meta]:
                    # mat = mat_qza.replace('.qza', '.tsv')
                    for metric in beta_metrics:
                        if metric in mat_qza:
                            break
                    if not isfile(mat_qza):
                        if not first_print:
                            print('Beta diversity, distances matrices must be generated already to automatise adonis\n'
                                  '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)')
                            first_print += 1
                        continue
                    # mat_pd = pd.read_csv(mat, header=0, index_col=0, sep='\t')
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
                                p = multiprocessing.Process(
                                    target=run_multi_adonis,
                                    args=(formula, odir, mat_qza, tsv,
                                          meta_pd, case, case_var, case_vals, out_sh)
                                )
                                p.start()
                                jobs.append(p)
            xpbs_call(out_sh, out_pbs, '%s.dns.%s' % (prjct_nm, dat), qiime_env,
                      '2', '1', '6', '2', 'gb')
            main_o.write('qsub %s\n' % out_pbs)

    for j in jobs:
        j.join()

    if len([1 for sh in all_shs for line in open(sh).readlines() if len(line.strip())]):
        if p_perm_groups:
            print("# Run Adonis (groups config in %s)" % p_perm_groups)
        else:
            print("# Run Adonis")
        print('[TO RUN] sh', main_sh)
