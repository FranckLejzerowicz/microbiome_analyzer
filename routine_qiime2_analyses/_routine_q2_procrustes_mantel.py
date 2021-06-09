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


def get_dat_idx(dat__, evaluation, datasets_filt, filt_only) -> (str, str):
    if evaluation:
        return dat__, ''
    else:
        if '__raref' in dat__:
            split = dat__.split('__raref')
            dat_ = '__raref'.join(split[:-1])
            raref = '_raref%s' % '__raref'.join(split[-1:])
        else:
            dat_ = dat__
            raref = ''
    if dat_ in datasets_filt and filt_only:
        dat = datasets_filt[dat_]
    else:
        dat = dat_
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


def check_dat_exists(betas, dat, missing_dats):
    if dat not in betas:
        if dat not in missing_dats:
            print('Dataset "%s" do not exist' % dat)
            missing_dats.add(dat)
            return 1
    return 0


def get_betas_raref(betas, dat, raref):
    metrics_groups_metas_qzas_dms_trees = betas[dat][0]
    if raref:
        if len(betas[dat]) > 1:
            if raref == '_raref':
                metrics_groups_metas_qzas_dms_trees = betas[dat][1]
            else:
                metrics_groups_metas_qzas_dms_trees = [
                    x for x in betas[dat] if x[list(x.keys())[0]][''][0].endswith('%s.tsv' % raref)][0]
    return metrics_groups_metas_qzas_dms_trees,


def get_meta_qza_dm_trees_d(meta_qza_dm_trees):
    meta_qza_dm_tree1s_d = {'': meta_qza_dm_trees[0]}
    if len(meta_qza_dm_trees) > 1:
        for (meta, qza, dm, tree) in meta_qza_dm_trees[1:]:
            meta_qza_dm_tree1s_d[dm.split('/')[-2]] = (meta, qza, dm, tree)
    return meta_qza_dm_tree1s_d


def run_procrustes(i_datasets_folder: str, datasets_filt: dict, p_procrustes: str,
                   betas: dict, force: bool, prjct_nm: str, qiime_env: str, chmod: str,
                   noloc: bool, split: bool, run_params: dict, filt_raref: str,
                   filt_only: bool, eval_depths: dict, jobs: bool, chunkit: int) -> None:
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
    missing_dats = set()
    for pair, (dat1_, dat2_) in procrustes_pairs.items():

        dat1, raref1 = get_dat_idx(dat1_, evaluation, datasets_filt, filt_only)
        dat2, raref2 = get_dat_idx(dat2_, evaluation, datasets_filt, filt_only)

        if check_dat_exists(betas, dat1, missing_dats) or check_dat_exists(betas, dat2, missing_dats):
            continue

        if evaluation:
            metrics_groups_metas_qzas_dms_trees1 = betas[dat1]
            metrics_groups_metas_qzas_dms_trees2 = betas[dat2]
        else:
            metrics_groups_metas_qzas_dms_trees1 = betas[dat1][0]
            metrics_groups_metas_qzas_dms_trees2 = betas[dat2][0]

        job_folder2 = get_job_folder(i_datasets_folder, 'procrustes%s/chunks/%s%s' % (evaluation, pair, filt_raref))
        if not split:
            out_sh = '%s/run_procrustes_%s%s_%s%s.sh' % (job_folder2, prjct_nm, evaluation, pair, filt_raref)

        for metric, groups_metas_qzas_dms_trees1 in metrics_groups_metas_qzas_dms_trees1.items():
            if split:
                out_sh = '%s/run_procrustes_%s%s_%s_%s%s.sh' % (job_folder2, prjct_nm, evaluation,
                                                                pair, metric, filt_raref)
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

                meta1, qza1, dm1, tree1 = groups_metas_qzas_dms_trees1[group1_][0]
                meta2, qza2, dm2, tree2 = groups_metas_qzas_dms_trees2[group2_][0]

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
                if len(common_sams) < 3:
                    continue

                meta_pd = meta_pd1.loc[meta_pd1.sample_name.isin(common_sams)]
                cases_dict = check_metadata_cases_dict(
                    meta1, meta_pd, dict(procrustes_subsets), 'procrustes')
                odir = get_analysis_folder(
                    i_datasets_folder, 'procrustes%s/%s%s/%s_vs_%s' % (
                        evaluation, pair, filt_raref, group1, group2))
                job_folder3 = get_job_folder(
                    i_datasets_folder, 'procrustes%s/chunks/%s%s/%s_vs_%s' % (
                        evaluation, pair, filt_raref, group1, group2))
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
    main_sh = write_main_sh(job_folder, '4_run_procrustes_%s%s%s' % (prjct_nm, evaluation, filt_raref), all_sh_pbs,
                            '%s.prcst%s%s' % (prjct_nm, evaluation, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        if p_procrustes and p_procrustes != 1:
            if p_procrustes.startswith('/panfs'):
                p_procrustes = p_procrustes.replace(os.getcwd(), '')
            print('# Procrustes (pairs and samples subsets config in %s)' % p_procrustes)
        else:
            print('# Procrustes')
        print_message('', 'sh', main_sh, jobs)

    dms_tab_pd = pd.DataFrame(
        dms_tab, columns=[
            'pair', 'dat1', 'dat2', 'metric',
            'group1', 'group2', 'case',
            'dm_out1', 'dm_out2',
        ]
    )

    odir = get_analysis_folder(i_datasets_folder, 'procrustes%s/R' % evaluation)
    out_Rs = glob.glob('%s/pairs_proscrustes_results%s%s*.tsv' % (odir, evaluation, filt_raref))
    if len(out_Rs):
        done_R = pd.concat([pd.read_table(x, sep=' ') for x in out_Rs])
        dms_tab_pd = dms_tab_pd.loc[
            ~dms_tab_pd[['dm_out1', 'dm_out2']].sum(1).isin(done_R[['f1', 'f2']].sum(1))
        ]

    if dms_tab_pd.shape[0]:
        fp_num = 0
        if len(out_Rs):
            print("out_Rs")
            print(out_Rs)
            last = sorted(out_Rs, key=lambda fp: int(fp.split('.tsv')[0].split('_')[-1]))
            fp_num = int(last[-1].split('.tsv')[0].split('_')[-1]) + 1

        dms_tab_fp = '%s/pairs%s%s_%s.tsv' % (odir, evaluation, filt_raref, fp_num)
        dms_tab_pd.to_csv(dms_tab_fp, index=False, sep='\t')
        out_R = '%s/pairs_proscrustes_results%s%s_%s.tsv' % (odir, evaluation, filt_raref, fp_num)
        job_folder = get_job_folder(i_datasets_folder, 'procrustes/R')
        R_script = '%s/4_run_procrustes_%s%s.R' % (job_folder, prjct_nm, filt_raref)
        with open(R_script, 'w') as o:
            o.write("library(vegan)\n")
            o.write("dms_files <- read.table('%s', h=T)\n" % dms_tab_fp)
            o.write("cols <- c('pair', 'd1', 'd2', 'g1', 'g2', 'case', 'metric', 'f1', 'f2', 'samples', 'M2', 'p-value')\n")
            o.write("res <- setNames(data.frame(matrix(ncol = 12, nrow = 0)), cols)\n")
            o.write("for (i in seq(1, dim(dms_files)[1])) {\n")
            o.write("    row <- as.vector(unlist(dms_files[i,]))\n")
            o.write("    pair <- row[1]\n")
            o.write("    d1 <- row[2]\n")
            o.write("    d2 <- row[3]\n")
            o.write("    group1 <- row[4]\n")
            o.write("    group2 <- row[5]\n")
            o.write("    case <- row[6]\n")
            o.write("    metric <- row[7]\n")
            o.write("    f1 <- row[8]\n")
            o.write("    f2 <- row[9]\n")
            o.write("    if (sum(file.exists(f1, f2)) == 2) {\n")
            o.write("        filin_tsv_pd1 <- read.csv(f1, header = TRUE, check.names=FALSE,\n")
            o.write("                                  row.names = 1, colClasses = 'character', sep = '\\t')\n")
            o.write("        filin_tsv_pd2 <- read.csv(f2, header = TRUE, check.names=FALSE,\n")
            o.write("                                  row.names = 1, colClasses = 'character', sep = '\\t')\n")
            o.write("        filin_tsv_pd1 <- data.matrix(filin_tsv_pd1)\n")
            o.write("        filin_tsv_pd2 <- data.matrix(filin_tsv_pd2)\n")
            o.write("        filin_tsv_pd1 <- filin_tsv_pd1[rownames(filin_tsv_pd2), rownames(filin_tsv_pd2)]\n")
            o.write("        # procrustes12 <- procrustes(filin_tsv_pd1, filin_tsv_pd2, kind=2, permutations=999)\n")
            o.write("        prtst <- protest(filin_tsv_pd1, filin_tsv_pd2, permutations = 999)\n")
            o.write("        n <- dim(filin_tsv_pd1)[1]\n")
            o.write("        res[i,] <- c(pair, d1, d2, group1, group2, case, metric, f1, f2, n, prtst$ss, prtst$signif)\n")
            o.write("    }\n")
            o.write("}\n")
            o.write("write.table(x = res, file = '%s')\n" % out_R)

        out_sh = '%s/4_run_procrustes_%s%s_R%s.sh' % (job_folder, prjct_nm, evaluation, filt_raref)
        out_pbs = '%s.pbs' % splitext(out_sh)[0]
        with open(out_sh, 'w') as o:
            o.write('R -f %s --vanilla\n' % R_script)

        run_xpbs(out_sh, out_pbs, '%s.prcrt%s.R%s' % (prjct_nm, evaluation, filt_raref), 'renv',
                 run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                 run_params["mem_num"], run_params["mem_dim"], chmod, 1,
                 '# Procrustes for stats in R (pairs and samples subsets config in %s)' % p_procrustes,
                 None, False, jobs)


def run_mantel(i_datasets_folder: str, datasets_filt: dict, p_mantel: str,
               betas: dict, force: bool, prjct_nm: str, qiime_env: str, chmod: str,
               noloc: bool, split: bool, run_params: dict, filt_raref: str, filt_only: bool,
               eval_depths: dict, jobs: bool, chunkit: int) -> None:
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
    missing_dats = set()
    for pair, (dat1_, dat2_) in mantel_pairs.items():

        dat1, raref1 = get_dat_idx(dat1_, evaluation, datasets_filt, filt_only)
        dat2, raref2 = get_dat_idx(dat2_, evaluation, datasets_filt, filt_only)

        if check_dat_exists(betas, dat1, missing_dats) or check_dat_exists(betas, dat2, missing_dats):
            continue

        if evaluation:
            metrics_groups_metas_qzas_dms_trees1 = betas[dat1]
            metrics_groups_metas_qzas_dms_trees2 = betas[dat2]
        else:
            metrics_groups_metas_qzas_dms_trees1 = betas[dat1][0]
            metrics_groups_metas_qzas_dms_trees2 = betas[dat2][0]

        job_folder2 = get_job_folder(i_datasets_folder, 'mantel%s/chunks/%s%s' % (evaluation, pair, filt_raref))
        if not split:
            out_sh = '%s/run_mantel_%s%s_%s%s.sh' % (job_folder2, prjct_nm, evaluation, pair, filt_raref)

        for metric, groups_metas_qzas_dms_trees1 in metrics_groups_metas_qzas_dms_trees1.items():
            if split:
                out_sh = '%s/run_mantel_%s%s_%s_%s%s.sh' % (job_folder2, prjct_nm, evaluation,
                                                            pair, metric, filt_raref)
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

                meta1, qza1, dm1, tree1 = groups_metas_qzas_dms_trees1[group1_][0]
                meta2, qza2, dm2, tree2 = groups_metas_qzas_dms_trees2[group2_][0]

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
                if len(common_sams) < 3:
                    continue

                meta_pd = meta_pd1.loc[meta_pd1.sample_name.isin(common_sams)]
                cases_dict = check_metadata_cases_dict(
                    meta1, meta_pd, dict(mantel_subsets), 'mantel')
                odir = get_analysis_folder(i_datasets_folder,
                                           'mantel%s/%s%s/%s_vs_%s' % (evaluation, pair, filt_raref, group1, group2))
                job_folder3 = get_job_folder(i_datasets_folder,
                                             'mantel%s/chunks/%s%s/%s_vs_%s' % (evaluation, pair, filt_raref, group1, group2))

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
    main_sh = write_main_sh(job_folder, '4_run_mantel_%s%s%s' % (prjct_nm, evaluation, filt_raref), all_sh_pbs,
                            '%s.mntl%s%s' % (prjct_nm, evaluation, filt_raref),
                            run_params["time"], run_params["n_nodes"], run_params["n_procs"],
                            run_params["mem_num"], run_params["mem_dim"],
                            qiime_env, chmod, noloc, jobs, chunkit)
    if main_sh:
        if p_mantel and p_mantel != 1:
            if p_mantel.startswith('/panfs'):
                p_mantel = p_mantel.replace(os.getcwd(), '')
            print('# Mantels (pairs and samples subsets config in %s)' % p_mantel)
        else:
            print('# Mantels')
        print_message('', 'sh', main_sh, jobs)
