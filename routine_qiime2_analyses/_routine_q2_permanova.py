# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import pkg_resources
from os.path import basename, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_main_cases_dict,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict,
    check_metadata_testing_groups
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_diversity_beta_group_significance,
    add_q2_types_to_meta
)


def run_single_perm(odir: str, subset: str, meta_pd: pd.DataFrame,
                    cur_sh: str, metric: str, case_: str, testing_group: str,
                    p_beta_type: tuple, qza: str, mat_qza: str, case_var: str,
                    case_vals: list, force: bool) -> None:
    """
    Run beta-group-significance: Beta diversity group significance.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-group-significance/
    (in-loop function).

    :param odir: output analysis directory.
    :param tsv: features table input to the beta diversity matrix.
    :param meta_pd: metadata table.
    :param cur_sh: input bash script file.
    :param case_:
    :param testing_group:
    :param mat_qza:
    :param case_var:
    :param case_vals:
    :param force: Force the re-writing of scripts for all commands.
    """
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        case = '%s__%s__%s' % (metric, case_, testing_group)
        case = case.replace(' ', '_')
        if subset:
            cur_rad = '%s/%s_%s_%s' % (odir, splitext(basename(qza))[0], subset, case)
        else:
            cur_rad = '%s/%s_%s' % (odir, splitext(basename(qza))[0], case)
        new_meta = '%s.meta' % cur_rad
        for beta_type in p_beta_type:
            new_qzv = '%s_%s.qzv' % (cur_rad, beta_type)
            new_html = '%s_%s.html' % (cur_rad, beta_type)
            new_cv = '%s_%s.cv' % (cur_rad, beta_type)
            new_mat_qza = odir + '/' + basename(mat_qza).replace('.qza', '_%s_DM.qza' % case)
            new_meta_pd = get_new_meta_pd(meta_pd, case, case_var, case_vals)
            add_q2_types_to_meta(new_meta_pd, new_meta, testing_group, new_cv)
            if force or not isfile(new_html):
                if len([x for x in new_meta_pd[testing_group].unique() if str(x) != 'nan']) > 1:
                    write_diversity_beta_group_significance(new_meta, mat_qza, new_mat_qza, testing_group,
                                                            beta_type, new_qzv, new_html, cur_sh_o)
                    remove = False
    if remove:
        os.remove(cur_sh)


def run_permanova(i_datasets_folder: str, betas: dict, main_testing_groups: tuple,
                  p_beta_type: tuple, datasets_rarefs: dict, p_perm_groups: str,
                  force: bool, prjct_nm: str, qiime_env: str, chmod: str,
                  noloc: bool, split: bool, run_params: dict,
                  filt_raref: str,  jobs: bool, chunkit: int) -> dict:
    """
    Run beta-group-significance: Beta diversity group significance.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-group-significance/
    Main per-dataset looper for the PERMANOVA tests on beta diversity matrices.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of datasets.
    :param betas: beta diversity matrices.
    :param main_testing_groups: groups to test.
    :param p_perm_groups: groups to subset.
    :param force: Force the re-writing of scripts for all commands.
    :param prjct_nm: Nick name for your project.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    """
    permanovas = {}
    job_folder2 = get_job_folder(i_datasets_folder, 'permanova/chunks')
    main_cases_dict = get_main_cases_dict(p_perm_groups)

    metric_check = set()
    all_sh_pbs = {}
    first_print = 0
    for dat, metric_groups_metas_qzas_dms_trees_ in betas.items():
        permanovas[dat] = []
        if not split:
            out_sh = '%s/run_beta_group_significance_%s_%s%s.sh' % (job_folder2, prjct_nm, dat, filt_raref)
        for idx, metric_groups_metas_qzas_dms_trees in enumerate(metric_groups_metas_qzas_dms_trees_):
            cur_depth = datasets_rarefs[dat][idx]
            odir = get_analysis_folder(i_datasets_folder, 'permanova/%s%s' % (dat, cur_depth))
            for metric, subset_files in metric_groups_metas_qzas_dms_trees.items():
                permanovas.setdefault(dat, []).append(metric)
                if split:
                    out_sh = '%s/run_beta_group_significance_%s_%s_%s%s.sh' % (job_folder2, prjct_nm,
                                                                               dat, metric, filt_raref)
                for subset, metas_qzas_mat_qzas_trees in subset_files.items():
                    for meta, qza, mat_qza, tree in metas_qzas_mat_qzas_trees:
                        if not isfile(mat_qza):
                            if not first_print:
                                print('Beta diversity, distances matrices must be generated already to automatise PERMANOVA\n'
                                      '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)')
                                first_print += 1
                            continue
                        if (dat, subset) not in metric_check:
                            meta_pd = read_meta_pd(meta)
                            meta_pd = meta_pd.set_index('sample_name')
                            cases_dict = check_metadata_cases_dict(meta, meta_pd, dict(main_cases_dict), 'PERMANOVA')
                            testing_groups = check_metadata_testing_groups(meta, meta_pd, main_testing_groups, 'PERMANOVA')
                            metric_check.add((dat, subset))

                        for case_var, case_vals_list in cases_dict.items():
                            testing_groups_case_var = list(set(testing_groups + [case_var]))
                            for case_vals in case_vals_list:
                                case = get_case(case_vals, case_var).replace(' ', '_')
                                for testing_group in testing_groups_case_var:
                                    if testing_group == 'ALL':
                                        continue
                                    cur_sh = '%s/run_beta_group_significance_%s%s_%s_%s_%s_%s%s.sh' % (
                                        job_folder2, dat, cur_depth, metric, subset, case, testing_group, filt_raref)
                                    cur_sh = cur_sh.replace(' ', '-')
                                    all_sh_pbs.setdefault((dat, out_sh), []).append(cur_sh)
                                    run_single_perm(odir, subset, meta_pd, cur_sh, metric, case,
                                                    testing_group, p_beta_type, qza, mat_qza,
                                                    case_var, case_vals, force)

    job_folder = get_job_folder(i_datasets_folder, 'permanova')
    main_sh = write_main_sh(job_folder, '3_run_beta_group_significance_%s%s' % (prjct_nm, filt_raref), all_sh_pbs,
                            '%s.prm%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        if p_perm_groups:
            print("# PERMANOVA (groups config in %s)" % p_perm_groups)
        else:
            print("# PERMANOVA")
        print_message('', 'sh', main_sh, jobs)

    return permanovas


def summarize_permanova(i_datasets_folder: str, permanovas: dict,
                        prjct_nm: str, qiime_env: str, chmod: str,
                        noloc: bool, split: bool, run_params: dict,
                        filt_raref: str,  jobs: bool, chunkit: int) -> dict:

    RESOURCES = pkg_resources.resource_filename("routine_qiime2_analyses", "resources")
    summarize_fp = '%s/summarize_permanovas.py' % RESOURCES

    all_sh_pbs = {}
    job_folder2 = get_job_folder(i_datasets_folder, 'permanova_summarize/chunks')
    for dat, metrics in permanovas.items():
        metrics = [x for x in [
            'aitchison',
            'jaccard',
            'braycurtis',
            'unweighted_unifrac',
            'weighted_unifrac'
        ] if x in metrics]
        permanovas[dat] = []
        out_sh = '%s/run_permanova_summarize_%s%s.sh' % (job_folder2, dat, filt_raref)
        out_py = '%s/run_permanova_summarize_%s%s.py' % (job_folder2, dat, filt_raref)
        with open(out_py, 'w') as o, open(summarize_fp) as f:
            for line in f:
                if 'ROUTINE_FOLDER' in line:
                    o.write(line.replace('ROUTINE_FOLDER', i_datasets_folder))
                elif 'METRICS' in line:
                    o.write(line.replace('METRICS', str(metrics)))
                elif 'DATASET' in line:
                    o.write(line.replace('DATASET', dat))
                else:
                    o.write(line)
        cur_sh = '%s/run_permanova_summarize_%s%s_tmp.sh' % (job_folder2, dat, filt_raref)
        with open(cur_sh, 'w') as o:
            o.write('python3 %s\n' % out_py)
        all_sh_pbs[(dat, out_sh)] = [cur_sh]

    job_folder = get_job_folder(i_datasets_folder, 'permanova_summarize')
    main_sh = write_main_sh(job_folder, '3_run_permanova_summarize%s' % filt_raref, all_sh_pbs,
                            '%s.prm%s' % (prjct_nm, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        print("# SUMMARIZE PERMANOVAS")
        print_message('', 'sh', main_sh, jobs)

    return permanovas
