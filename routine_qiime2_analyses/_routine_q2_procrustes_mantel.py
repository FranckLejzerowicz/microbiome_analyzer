import os
import glob
import itertools
import pandas as pd
from os.path import dirname, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_procrustes_mantel_dicts,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_procrustes_mantel
)


def run_single_procrustes_mantel(procrustes_mantel: str, odir: str, dm1: str, dm2: str,
                                 meta_pd: pd.DataFrame, dm_out1: str, dm_out2: str,
                                 output: str, cur_sh: str, cur: str, case_var: str,
                                 case_vals: list, force: bool) -> None:
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(output):
            new_meta_pd = get_new_meta_pd(meta_pd, cur, case_var, case_vals)
            common_meta_fp = '%s/meta_%s.tsv' % (odir, cur)
            new_meta_pd.to_csv(common_meta_fp, index=False, sep='\t')
            if new_meta_pd.shape[0]:
                write_procrustes_mantel(
                    procrustes_mantel, common_meta_fp, dm1, dm2,
                    dm_out1, dm_out2, output, cur_sh_o)
                remove = False
    if remove:
        os.remove(cur_sh)



def get_dat_idx(dat_) -> (str, str):
    if '__raref' in dat_:
        split = dat_.split('__raref')
        dat = '__raref'.join(split[:-1])
        raref = '_raref%s' % '__raref'.join(split[-1:])
    else:
        dat = dat_
        raref = ''
    return dat, raref


def get_dm_meta(dat, dm, meta, raref, metric, i_datasets_folder, skip):
    dm_rgx = '%s%s*/*%s_DM.qza' % (dirname(dm), raref, metric)
    dm_rgx_glob = glob.glob(dm_rgx)
    if len(dm_rgx_glob) >= 1:
        dm = sorted(dm_rgx_glob)[0]
    else:
        skip += 1
    meta_dir = get_analysis_folder(i_datasets_folder, 'rarefy/%s' % dat)
    meta_rgx = '%s/meta_%s%s*tsv' % (meta_dir, dat, raref)
    meta_rgx_glob = glob.glob(meta_rgx)
    if len(meta_rgx_glob) >= 1:
        meta = sorted(meta_rgx_glob)[0]
    else:
        skip += 1
    return dm, meta


def run_procrustes(i_datasets_folder: str, p_procrustes: str, betas: dict,
                   force: bool, prjct_nm: str, qiime_env: str, chmod: str,
                   noloc: bool, split: bool, run_params: dict, filt_raref: str,
                   eval_depths: dict) -> None:
    """
    """
    evaluation = ''
    if eval_depths:
        evaluation = '_eval'
        procrustes_pairs = {}
        for dat, depths in eval_depths.items():
            sorted_depths = sorted(depths, key=lambda x: int(x.split('_')[-1]))
            for idx, x in enumerate(sorted_depths[:-1]):
                y = sorted_depths[(idx+1)]
                n0 = x.split('_')[-1]
                n1 = y.split('_')[-1]
                procrustes_pairs['%s_%s' % (n0, n1)] = [x, y]
        procrustes_subsets = {'ALL': [[]]}
    else:
        procrustes_pairs, procrustes_subsets = get_procrustes_mantel_dicts(p_procrustes)

    get_job_folder(i_datasets_folder, 'procrustes%s' % evaluation)

    dms_tab = []
    all_sh_pbs = {}
    for pair, (dat1_, dat2_) in procrustes_pairs.items():

        if evaluation:
            dat1, dat2 = dat1_, dat2_
            metrics_groups_metas_qzas_dms_trees1 = betas[dat1]
            metrics_groups_metas_qzas_dms_trees2 = betas[dat2]
        else:
            dat1, raref1 = get_dat_idx(dat1_)
            dat2, raref2 = get_dat_idx(dat2_)
            metrics_groups_metas_qzas_dms_trees1 = betas[dat1][0]
            metrics_groups_metas_qzas_dms_trees2 = betas[dat2][0]

        job_folder2 = get_job_folder(i_datasets_folder, 'procrustes%s/chunks/%s' % (evaluation, pair))
        if not split:
            out_sh = '%s/run_procrustes%s_%s%s.sh' % (job_folder2, evaluation, pair, filt_raref)

        for metric, groups_metas_qzas_dms_trees1 in metrics_groups_metas_qzas_dms_trees1.items():
            if split:
                out_sh = '%s/run_procrustes%s_%s_%s%s.sh' % (job_folder2, evaluation, pair, metric, filt_raref)
            if metric not in metrics_groups_metas_qzas_dms_trees2:
                continue
            groups_metas_qzas_dms_trees2 = metrics_groups_metas_qzas_dms_trees2[metric]
            groups1 = sorted(groups_metas_qzas_dms_trees1.keys())
            groups2 = sorted(groups_metas_qzas_dms_trees2.keys())
            for (group1_, group2_) in itertools.product(*[groups1, groups2]):
                if group1_ == '':
                    group1 = 'full'
                else:
                    group1 = group1_
                if group2_ == '':
                    group2 = 'full'
                else:
                    group2 = group2_

                meta1, qza1, dm1, tree1 = groups_metas_qzas_dms_trees1[group1_]
                meta2, qza2, dm2, tree2 = groups_metas_qzas_dms_trees2[group2_]

                skip = 0
                if not evaluation:
                    if '__raref' in dat1_:
                        dm1, meta1 = get_dm_meta(dat1, dm1, meta1, raref1, metric, i_datasets_folder, skip)
                    if '__raref' in dat2_:
                        dm2, meta2 = get_dm_meta(dat2, dm2, meta2, raref2, metric, i_datasets_folder, skip)
                if skip:
                    print('[Proscustes] One desired rarefaction depth not run (pair %s)' % pair)
                    continue

                meta_pd1 = read_meta_pd(meta1)
                meta_pd2 = read_meta_pd(meta2)
                common_sams = list(
                    set(meta_pd1.sample_name) &
                    set(meta_pd2.sample_name)
                )
                meta_pd = meta_pd1.loc[meta_pd1.sample_name.isin(common_sams)]
                cases_dict = check_metadata_cases_dict(
                    meta1, meta_pd, dict(procrustes_subsets), 'procrustes')
                odir = get_analysis_folder(i_datasets_folder,
                                           'procrustes%s/%s/%s_vs_%s' % (evaluation, pair, group1, group2))
                job_folder3 = get_job_folder(i_datasets_folder,
                                             'procrustes%s/chunks/%s/%s_vs_%s' % (evaluation, pair, group1, group2))
                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        case_ = get_case(case_vals, case_var).replace(' ', '_')
                        cur = '%s__%s' % (metric, case_)
                        cur_sh = '%s/run_procrustes%s_%s%s.sh' % (job_folder3, evaluation, cur, filt_raref)
                        cur_sh = cur_sh.replace(' ', '-')
                        all_sh_pbs.setdefault((pair, out_sh), []).append(cur_sh)

                        dm_out1 = '%s/dm_%s__%s_DM.qza' % (odir, dat1_, cur)
                        dm_out2 = '%s/dm_%s__%s_DM.qza' % (odir, dat2_, cur)
                        dm_out1_tsv = '%s.tsv' % splitext(dm_out1)[0]
                        dm_out2_tsv = '%s.tsv' % splitext(dm_out2)[0]
                        biplot = '%s/procrustes%s_%s__%s__%s.qzv' % (odir, evaluation, dat1_, dat2_, cur)
                        run_single_procrustes_mantel('procrustes', odir, dm1, dm2, meta_pd, dm_out1, dm_out2,
                                                     biplot, cur_sh, cur, case_var, case_vals, force)
                        dms_tab.append([pair, dat1_, dat2_,
                                        group1, group2, case_, metric,
                                        dm_out1_tsv, dm_out2_tsv])

    job_folder = get_job_folder(i_datasets_folder, 'procrustes%s' % evaluation)
    main_sh = write_main_sh(job_folder, '4_run_procrustes%s%s' % (evaluation, filt_raref), all_sh_pbs,
                            '%s.prcst%s%s' % (prjct_nm, evaluation, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_procrustes and p_procrustes != 1:
            if p_procrustes.startswith('/panfs'):
                p_procrustes = p_procrustes.replace(os.getcwd(), '')
            print('# Procrustes (pairs and samples subsets config in %s)' % p_procrustes)
        else:
            print('# Procrustes')
        print_message('', 'sh', main_sh)

    dms_tab_pd = pd.DataFrame(
        dms_tab, columns=[
            'pair', 'dat1', 'dat2', 'metric',
            'group1', 'group2', 'case',
            'dm_out1', 'dm_out2',
        ]
    )
    odir = get_analysis_folder(i_datasets_folder, 'procrustes%s/R' % evaluation)
    dms_tab_fp = '%s/pairs%s.tsv' % (odir, evaluation)
    dms_tab_pd.to_csv(dms_tab_fp, index=False, sep='\t')

    out_R = '%s/pairs_proscrustes_results%s.tsv' % (odir, evaluation)
    if not isfile(out_R):
        job_folder = get_job_folder(i_datasets_folder, 'procrustes/R')
        R_script = '%s/4_run_procrustes%s.R' % (job_folder, filt_raref)
        with open(R_script, 'w') as o:
            o.write("library(vegan)\n")
            o.write("dms_files <- read.table('%s', h=T)\n" % dms_tab_fp)
            o.write("cols <- c('comparison', 'd1', 'd2', 'g1', 'g2', 'case', 'metric', 'f1', 'f2', 'M2', 'signif')\n")
            o.write("res <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), cols)\n")
            o.write("for (i in seq(1, dim(dms_files)[1])) {\n")
            o.write("    row <- as.vector(unlist(dms_files[i,]))\n")
            o.write("    com <- row[1]\n")
            o.write("    d1 <- row[2]\n")
            o.write("    d2 <- row[3]\n")
            o.write("    group1 <- row[4]\n")
            o.write("    group2 <- row[5]\n")
            o.write("    case <- row[6]\n")
            o.write("    metric <- row[7]\n")
            o.write("    f1 <- row[8]\n")
            o.write("    f2 <- row[9]\n")
            o.write("    filin_tsv_pd1 <- read.csv(f1, header = TRUE, check.names=FALSE,\n")
            o.write("                              row.names = 1, colClasses = 'character', sep = '\\t')\n")
            o.write("    filin_tsv_pd2 <- read.csv(f2, header = TRUE, check.names=FALSE,\n")
            o.write("                              row.names = 1, colClasses = 'character', sep = '\\t')\n")
            o.write("    filin_tsv_pd1 <- data.matrix(filin_tsv_pd1)\n")
            o.write("    filin_tsv_pd2 <- data.matrix(filin_tsv_pd2)\n")
            o.write("    filin_tsv_pd1 <- filin_tsv_pd1[rownames(filin_tsv_pd2), rownames(filin_tsv_pd2)]\n")
            o.write("    # procrustes12 <- procrustes(filin_tsv_pd1, filin_tsv_pd2, kind=2, permutations=999)\n")
            o.write("    prtst <- protest(filin_tsv_pd1, filin_tsv_pd2, permutations = 999)\n")
            o.write("    res[i,] <- c(com, d1, d2, group1, group2, case, metric, f1, f2, prtst$ss, prtst$signif)\n")
            o.write("}\n")
            o.write("write.table(x = res, file = '%s')\n" % out_R)

        out_sh = '%s/4_run_procrustes%s_R%s.sh' % (job_folder, evaluation, filt_raref)
        out_pbs = '%s.pbs' % splitext(out_sh)[0]
        with open(out_sh, 'w') as o:
            o.write('R -f %s --vanilla\n' % R_script)

        run_xpbs(out_sh, out_pbs, '%s.prcrt%s.R%s' % (prjct_nm, evaluation, filt_raref), 'renv',
                 run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                 run_params["mem_num"], run_params["mem_dim"], chmod, 1,
                 '# Procrustes for stats in R (pairs and samples subsets config in %s)' % p_procrustes,
                 None, False)
    else:
        print('%s exists (remove to re-run)' % out_R)


def run_mantel(i_datasets_folder: str, p_mantel: str, betas: dict,
               force: bool, prjct_nm: str, qiime_env: str, chmod: str,
               noloc: bool, split: bool, run_params: dict, filt_raref: str,
               eval_depths: dict) -> None:
    """
    """
    evaluation = ''
    if eval_depths:
        evaluation = '_eval'
        mantel_pairs = {}
        for dat, depths in eval_depths.items():
            sorted_depths = sorted(depths, key=lambda x: int(x.split('_')[-1]))
            for idx, x in enumerate(sorted_depths[:-1]):
                y = sorted_depths[(idx+1)]
                n0 = x.split('_')[-1]
                n1 = y.split('_')[-1]
                mantel_pairs['%s_%s' % (n0, n1)] = [x, y]
        mantel_subsets = {'ALL': [[]]}
    else:
        mantel_pairs, mantel_subsets = get_procrustes_mantel_dicts(p_mantel)

    get_job_folder(i_datasets_folder, 'mantel%s' % evaluation)

    all_sh_pbs = {}
    for pair, (dat1_, dat2_) in mantel_pairs.items():

        if evaluation:
            dat1, dat2 = dat1_, dat2_
            metrics_groups_metas_qzas1 = betas[dat1]
            metrics_groups_metas_qzas2 = betas[dat2]
        else:
            dat1, raref1 = get_dat_idx(dat1_)
            dat2, raref2 = get_dat_idx(dat2_)
            metrics_groups_metas_qzas1 = betas[dat1][0]
            metrics_groups_metas_qzas2 = betas[dat2][0]

        job_folder2 = get_job_folder(i_datasets_folder, 'mantel%s/chunks/%s' % (evaluation, pair))
        if not split:
            out_sh = '%s/run_mantel%s_%s%s.sh' % (job_folder2, evaluation, pair, filt_raref)

        for metric, groups_metas_qzas1 in metrics_groups_metas_qzas1.items():
            if split:
                out_sh = '%s/run_mantel%s_%s_%s%s.sh' % (job_folder2, evaluation, pair, metric, filt_raref)
            if metric not in metrics_groups_metas_qzas2:
                continue
            groups_metas_qzas2 = metrics_groups_metas_qzas2[metric]
            groups1 = sorted(groups_metas_qzas1.keys())
            groups2 = sorted(groups_metas_qzas2.keys())
            for (group1_, group2_) in itertools.product(*[groups1, groups2]):
                if group1_ == '':
                    group1 = 'full'
                else:
                    group1 = group1_
                if group2_ == '':
                    group2 = 'full'
                else:
                    group2 = group2_

                meta1, qza1, dm1 = groups_metas_qzas1[group1_]
                meta2, qza2, dm2 = groups_metas_qzas2[group2_]

                skip = 0
                if not evaluation:
                    if '__raref' in dat1_:
                        dm1, meta1 = get_dm_meta(dat1, dm1, meta1, raref1, metric, i_datasets_folder, skip)
                    if '__raref' in dat2_:
                        dm2, meta2 = get_dm_meta(dat2, dm2, meta2, raref2, metric, i_datasets_folder, skip)
                if skip:
                    print('[Mantels] One desired rarefaction depth not run (pair %s)' % pair)
                    continue

                meta_pd1 = read_meta_pd(meta1)
                meta_pd2 = read_meta_pd(meta2)
                common_sams = list(
                    set(meta_pd1.sample_name) &
                    set(meta_pd2.sample_name)
                )
                meta_pd = meta_pd1.loc[meta_pd1.sample_name.isin(common_sams)]
                cases_dict = check_metadata_cases_dict(
                    meta1, meta_pd, dict(mantel_subsets), 'mantel')
                odir = get_analysis_folder(i_datasets_folder,
                                           'mantel%s/%s/%s_vs_%s' % (evaluation, pair, group1, group2))
                job_folder3 = get_job_folder(i_datasets_folder,
                                             'mantel%s/chunks/%s/%s_vs_%s' % (evaluation, pair, group1, group2))

                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        case_ = get_case(case_vals, case_var).replace(' ', '_')
                        cur = '%s__%s' % (metric, case_)
                        cur_sh = '%s/run_mantel%s_%s%s.sh' % (job_folder3, evaluation, cur, filt_raref)
                        cur_sh = cur_sh.replace(' ', '-')
                        all_sh_pbs.setdefault((pair, out_sh), []).append(cur_sh)

                        dm_out1 = '%s/dm_%s__%s_DM.qza' % (odir, dat1_, cur)
                        dm_out2 = '%s/dm_%s__%s_DM.qza' % (odir, dat2_, cur)
                        mantel_out = '%s/mantel%s_%s__%s__%s.qzv' % (odir, evaluation, dat1_, dat2_, cur)
                        run_single_procrustes_mantel('mantel', odir, dm1, dm2, meta_pd, dm_out1, dm_out2,
                                                     mantel_out, cur_sh, cur, case_var, case_vals, force)

    job_folder = get_job_folder(i_datasets_folder, 'mantel%s' % evaluation)
    main_sh = write_main_sh(job_folder, '4_run_mantel%s%s' % (evaluation, filt_raref), all_sh_pbs,
                            '%s.mntl%s%s' % (prjct_nm, evaluation, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_mantel and p_mantel != 1:
            if p_mantel.startswith('/panfs'):
                p_mantel = p_mantel.replace(os.getcwd(), '')
            print('# Mantels (pairs and samples subsets config in %s)' % p_mantel)
        else:
            print('# Mantels')
        print_message('', 'sh', main_sh)
