import os
import glob
import itertools
import pandas as pd
from os.path import dirname, isfile, splitext

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs, print_message
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_job_folder,
    get_analysis_folder,
    get_procrustes_dicts,
    write_main_sh,
    read_meta_pd
)
from routine_qiime2_analyses._routine_q2_metadata import (
    check_metadata_cases_dict
)
from routine_qiime2_analyses._routine_q2_cmds import (
    get_new_meta_pd, get_case,
    write_procrustes
)


def run_single_procrustes(odir: str, dm1: str, dm2: str, meta_pd: pd.DataFrame,
                          dm_out1: str, dm_out2: str, biplot: str, cur_sh: str,
                          cur: str, case_var: str, case_vals: list, force: bool) -> None:
    remove = True
    with open(cur_sh, 'w') as cur_sh_o:
        if force or not isfile(biplot):
            new_meta_pd = get_new_meta_pd(meta_pd, cur, case_var, case_vals)
            common_meta_fp = '%s/meta_%s.tsv' % (odir, cur)
            new_meta_pd.to_csv(common_meta_fp, index=False, sep='\t')
            if new_meta_pd.shape[0]:
                write_procrustes(common_meta_fp, dm1, dm2,
                                 dm_out1, dm_out2,
                                 biplot, cur_sh_o)
                remove = False
    if remove:
        os.remove(cur_sh)


def run_procrustes(i_datasets_folder: str, datasets: dict, datasets_filt: dict,
                   p_procrustes: str, betas: dict, force: bool, prjct_nm: str,
                   qiime_env: str, chmod: str, noloc: bool, split: bool) -> None:
    """
    """
    job_folder = get_job_folder(i_datasets_folder, 'procrustes')
    procrustes_pairs, procrustes_subsets = get_procrustes_dicts(p_procrustes)

    dms_tab = []
    all_sh_pbs = {}
    run_pbs = '%s/4_run_procrustes.sh' % job_folder
    for pair, (dat1_, dat2_) in procrustes_pairs.items():

        if dat1_.endswith('__raref') and dat1_ not in datasets:
            dat1 = dat1_.split('__raref')[0]
        else:
            dat1 = dat1_
        if dat2_.endswith('__raref') and dat2_ not in datasets:
            dat2 = dat2_.split('__raref')[0]
        else:
            dat2 = dat2_

        job_folder2 = get_job_folder(i_datasets_folder, 'procrustes/chunks/%s' % pair)
        if not split:
            out_sh = '%s/run_procrustes_%s.sh' % (job_folder2, pair)
        metrics_groups_metas_qzas1 = betas[dat1]
        for metric, groups_metas_qzas1 in metrics_groups_metas_qzas1.items():
            if split:
                out_sh = '%s/run_procrustes_%s_%s.sh' % (job_folder2, pair, metric)
            if metric not in betas[dat2] or metric not in betas[dat1]:
                continue
            groups_metas_qzas2 = betas[dat2][metric]
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
                meta1, qza1, dm1 = betas[dat1][metric][group1_]
                meta2, qza2, dm2 = betas[dat2][metric][group2_]

                if dat1_.endswith('__raref'):
                    dm1_rgx = glob.glob('%s/tab_%s_raref*_%s_DM.qza' % (dirname(dm1), dat1, metric))
                    if len(dm1_rgx) == 1:
                        dm1 = dm1_rgx[0]
                    meta_dir = get_analysis_folder(i_datasets_folder, 'rarefy/%s' % dat1)
                    meta1_rgx = glob.glob('%s/meta_%s_raref*tsv' % (meta_dir, dat1))
                    if len(meta1_rgx) >= 1:
                        meta1 = meta1_rgx[0]
                if dat2_.endswith('__raref'):
                    dm2_rgx = glob.glob('%s/tab_%s_raref*_%s_DM.qza' % (dirname(dm2), dat2, metric))
                    if len(dm2_rgx) == 1:
                        dm2 = dm2_rgx[0]
                    meta_dir = get_analysis_folder(i_datasets_folder, 'rarefy/%s' % dat2)
                    meta2_rgx = glob.glob('%s/meta_%s_raref*tsv' % (meta_dir, dat2))
                    if len(meta2_rgx) >= 1:
                        meta2 = meta2_rgx[0]

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
                                           'procrustes/%s/%s_vs_%s' % (pair, group1, group2))
                job_folder3 = get_job_folder(i_datasets_folder,
                                             'procrustes/chunks/%s/%s_vs_%s' % (pair, group1, group2))
                for case_var, case_vals_list in cases_dict.items():
                    for case_vals in case_vals_list:
                        case_ = get_case(case_vals, case_var).replace(' ', '_')
                        cur = '%s__%s' % (metric, case_)
                        cur_sh = '%s/run_procrustes_%s.sh' % (job_folder3, cur)
                        cur_sh = cur_sh.replace(' ', '-')
                        all_sh_pbs.setdefault((pair, out_sh), []).append(cur_sh)

                        dm_out1 = '%s/dm_%s__%s_DM.qza' % (odir, dat1_, cur)
                        dm_out2 = '%s/dm_%s__%s_DM.qza' % (odir, dat2_, cur)
                        dm_out1_tsv = '%s.tsv' % splitext(dm_out1)[0]
                        dm_out2_tsv = '%s.tsv' % splitext(dm_out2)[0]
                        biplot = '%s/procrustes_%s__%s__%s.qzv' % (odir, dat1_, dat2_, cur)
                        run_single_procrustes(odir, dm1, dm2, meta_pd, dm_out1, dm_out2,
                                              biplot, cur_sh, cur, case_var, case_vals, force)
                        dms_tab.append([pair, dat1_, dat2_,
                                        group1, group2, case_, metric,
                                        dm_out1_tsv, dm_out2_tsv])

    job_folder = get_job_folder(i_datasets_folder, 'procrustes')
    main_sh = write_main_sh(job_folder, 'run_procrustes', all_sh_pbs,
                            '%s.prcst' % prjct_nm, '2', '1', '1', '1', 'gb',
                            qiime_env, chmod, noloc)
    if main_sh:
        if p_procrustes:
            if p_procrustes.startswith('/panfs'):
                p_procrustes = p_procrustes.replace(os.getcwd(), '')
            print('# Procrustes (pairs and samples subsets config in %s)' % p_procrustes)
        else:
            print('# Procrustes')
        print_message('', 'sh', main_sh)

    dms_tab_pd = pd.DataFrame(
        dms_tab, columns = [
            'pair', 'dat1', 'dat2', 'metric',
            'group1', 'group2', 'case',
            'dm_out1', 'dm_out2',
        ]
    )
    odir = get_analysis_folder(i_datasets_folder, 'procrustes/R')
    dms_tab_fp = '%s/pairs.tsv' % odir
    dms_tab_pd.to_csv(dms_tab_fp, index=False, sep='\t')

    out_R = '%s/pairs_proscrustes_results.tsv' % odir
    if not isfile(out_R):
        job_folder = get_job_folder(i_datasets_folder, 'procrustes/R')
        R_script = '%s/4_run_procrustes.R' % job_folder
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
            o.write("    filin_tsv_pd1 <- read.csv(f1, header = TRUE, check.names=FALSE, row.names = 1, sep = '\\t')\n")
            o.write("    filin_tsv_pd2 <- read.csv(f2, header = TRUE, check.names=FALSE, row.names = 1, sep = '\\t')\n")
            o.write("    filin_tsv_pd1 <- as.matrix(filin_tsv_pd1)\n")
            o.write("    filin_tsv_pd2 <- as.matrix(filin_tsv_pd2)\n")
            o.write("    # procrustes12 <- procrustes(filin_tsv_pd1, filin_tsv_pd2, kind=2, permutations=999)\n")
            o.write("    prtst <- protest(filin_tsv_pd1, filin_tsv_pd2, permutations = 999)\n")
            o.write("    res[i,] <- c(com, d1, d2, group1, group2, case, metric, f1, f2, prtst$ss, prtst$signif)\n")
            o.write("}\n")
            o.write("write.table(x = res, file = '%s')\n" % out_R)

        out_sh = '%s/4_run_procrustes_R.sh' % job_folder
        out_pbs = '%s.pbs' % splitext(out_sh)[0]
        with open(out_sh, 'w') as o:
            o.write('R -f %s --vanilla\n' % R_script)

        run_xpbs(out_sh, out_pbs, '%s.prcrt.R' % prjct_nm,
                 'renv', '4', '1', '1', '1', 'gb', chmod, 1,
                 '# Procrustes for stats in R (pairs and samples subsets config in %s)' % p_procrustes,
                 None, False)




