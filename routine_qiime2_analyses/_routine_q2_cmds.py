

import pandas as pd
from typing import TextIO
from os.path import isfile

from routine_qiime2_analyses._routine_q2_io_utils import run_import


def write_diversity_beta(out_fp: str, datasets_phylo: dict, trees: dict, dat: str,
                         qza: str, metric: str, cur_sh: TextIO) -> bool:
    if 'unifrac' in metric:
        if not datasets_phylo[dat][0]:
            return True
        cmd = 'qiime diversity beta-phylogenetic \\ \n'
        if datasets_phylo[dat][1]:
            cmd += '--i-table %s \\ \n' % trees[dat][0]
        else:
            cmd += '--i-table %s \\ \n' % qza
        cmd += '--i-phylogeny %s \\ \n' % trees[dat][1]
    else:
        cmd = 'qiime diversity beta \\ \n'
        cmd += '--i-table %s \\ \n' % qza
    cmd += '--p-metric %s \\ \n' % metric
    cmd += '--p-n-jobs 1 \\ \n'
    cmd += '--o-distance-matrix %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)
    return False


def write_diversity_pcoa(DM: str, out_pcoa: str, cur_sh: TextIO) -> None:
    cmd = 'qiime diversity pcoa \\ \n'
    cmd += '--i-distance-matrix %s \\ \n' % DM
    cmd += '--o-pcoa %s\n' % out_pcoa
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_emperor(meta: str, pcoa: str, out_plot: str, cur_sh: TextIO) -> None:
    cmd = 'qiime emperor plot \\ \n'
    cmd += '--m-metadata-file %s \\ \n' % meta
    cmd += '--i-pcoa %s \\ \n' % pcoa
    cmd += '--o-visualization %s\n' % out_plot
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_seqs_fasta(out_fp_seqs_fasta: str, out_fp_seqs_qza: str,
                     tsv_pd: pd.DataFrame, cur_sh: TextIO) -> None:
    with open(out_fp_seqs_fasta, 'w') as fas_o:
        for seq in tsv_pd.index:
            fas_o.write('>%s\n%s\n' % (seq.strip(), seq.strip()))
    cmd = run_import(out_fp_seqs_fasta, out_fp_seqs_qza, 'FeatureData[Sequence]')
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_fragment_insertion(out_fp_seqs_qza: str, ref_tree_qza: str,
                             out_fp_sepp_tree: str, out_fp_sepp_plac: str,
                             qza: str, qza_in: str, qza_out: str,
                             cur_sh: TextIO) -> None:
    cmd = 'qiime fragment-insertion sepp \\ \n'
    cmd += '--i-representative-sequences %s \\ \n' % out_fp_seqs_qza
    cmd += '--i-reference-database %s \\ \n' % ref_tree_qza
    cmd += '--o-tree %s \\ \n' % out_fp_sepp_tree
    cmd += '--o-placements %s \\ \n' % out_fp_sepp_plac
    cmd += '--p-threads 24\n'
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)
    cmd = 'qiime fragment-insertion filter-features \\ \n'
    cmd += '--i-table %s \\ \n' % qza
    cmd += '--i-tree %s \\ \n' % out_fp_sepp_tree
    cmd += '--o-filtered-table %s \\ \n' % qza_in
    cmd += '--o-removed-table %s\n' % qza_out
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_deicode_biplot(qza: str, new_meta: str, new_qza: str, ordi_qza: str,
                         new_mat_qza: str, ordi_qzv: str, cur_sh: TextIO) -> None:
    cmd = 'qiime feature-table filter-samples \\ \n'
    cmd += '--i-table %s \\ \n' % qza
    cmd += '--m-metadata-file %s \\ \n' % new_meta
    cmd += '--o-filtered-table %s \\ \n' % new_qza
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)
    cmd = 'qiime deicode rpca \\ \n'
    cmd += '--i-table %s \\ \n' % new_qza
    cmd += '--p-min-feature-count 10 \\ \n'
    cmd += '--p-min-sample-count 500 \\ \n'
    cmd += '--o-biplot %s \\ \n' % ordi_qza
    cmd += '--o-distance-matrix %s\n' % new_mat_qza
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)
    cmd = 'qiime emperor biplot \\ \n'
    cmd += '--i-biplot %s \\ \n' % ordi_qza
    cmd += '--m-sample-metadata-file %s \\ \n' % new_meta
    cmd += '--o-visualization %s \\ \n' % ordi_qzv
    cmd += '--p-number-of-features 20\n'
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def add_q2_types_to_meta(new_meta_pd: pd.DataFrame, new_meta: str) -> None:
    q2types = pd.DataFrame(
        [(['#q2:types'] + ['categorical'] * new_meta_pd.shape[1])],
        columns=new_meta_pd.reset_index().columns.tolist(),
    )
    q2types.rename(columns={q2types.columns.tolist()[0]: new_meta_pd.index.name}, inplace=True)
    q2types.set_index(new_meta_pd.index.name, inplace=True)
    new_meta_pd = pd.concat([q2types, new_meta_pd])
    new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')


def check_absence_mat(mat_qzas: list, first_print: int, analysis: str) -> bool:
    presence_mat = [mat_qza for mat_qza in mat_qzas if isfile(mat_qza)]
    if not presence_mat:
        if not first_print:
            print('Beta diversity, distances matrices must be generated already to automatise %s\n'
                  '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)' % analysis)
            first_print += 1
        return True
    return False


def write_diversity_beta_group_significance(new_meta: str, mat_qza: str, new_mat_qza: str,
                                            qza: str, new_qza: str, testing_group: str,
                                            new_qzv: str, cur_sh: TextIO) -> None:
    cmd = 'qiime diversity filter-distance-matrix \\ \n'
    cmd += '--m-metadata-file %s \\ \n' % new_meta
    cmd += '--i-distance-matrix %s \\ \n' % mat_qza
    cmd += '--o-filtered-distance-matrix %s\n' % new_mat_qza
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)
    cmd = 'qiime feature-table filter-samples \\ \n'
    cmd += '--i-table %s \\ \n' % qza
    cmd += '--m-metadata-file %s \\ \n' % new_meta
    cmd += '--o-filtered-table %s\n' % new_qza
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)
    cmd = 'qiime diversity beta-group-significance \\ \n'
    cmd += '--i-distance-matrix %s \\ \n' % new_mat_qza
    cmd += '--m-metadata-file %s \\ \n' % new_meta
    cmd += '--m-metadata-column %s \\ \n' % testing_group
    cmd += '--p-permutations 2999 \\ \n'
    cmd += '--o-visualization %s\n' % new_qzv
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def write_diversity_adonis(new_meta: str, mat_qza: str, new_mat_qza: str,
                           qza: str, new_qza: str, formula: str,
                           new_qzv: str, cur_sh: TextIO) -> None:
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


def get_metric(metrics: list, file_name: str) -> str:
    for metric in metrics:
        if metric in file_name:
            return metric
    return ''


def get_case(case_vals: list, metric: str, case_var: str, form: str = None) -> str:
    if len(case_vals):
        case = '%s_%s_%s' % (metric, case_var, '-'.join(
            [x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
    else:
        case = '%s_%s' % (metric, case_var)
    if form:
        case = '%s_%s' % (case, form)
    case = case.replace('__', '_')
    return case


def write_longitudinal_volatility(out_fp: str, meta_alphas: str,
                                  time_point: str, cur_sh: TextIO) -> None:
    cmd = 'qiime longitudinal volatility \\ \n'
    cmd += '--m-metadata-file %s \\ \n' % meta_alphas
    cmd += '--p-state-column "%s" \\ \n' % time_point
    cmd += '--p-individual-id-column "host_subject_id"'
    cmd += '--o-visualization %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_diversity_alpha_correlation(out_fp: str, qza: str, method: str,
                                      meta: str, cur_sh: TextIO) -> None:
    cmd = 'qiime diversity alpha-correlation \\ \n'
    cmd += '--i-alpha-diversity %s \\ \n' % qza
    cmd += '--p-method %s \\ \n' % method
    cmd += '--m-metadata-file %s \\ \n' % meta
    cmd += '--o-visualization %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_diversity_alpha(out_fp: str, datasets_phylo: dict, trees: dict, dat: str,
                          qza: str, metric: str, cur_sh: TextIO) -> bool:
    if metric in ['faith_pd']:
        if not datasets_phylo[dat][0]:
            return True
        cmd = 'qiime diversity alpha-phylogenetic \\ \n'
        if datasets_phylo[dat][1]:
            cmd += '--i-table %s \\ \n' % trees[dat][0]
        else:
            cmd += '--i-table %s \\ \n' % qza
        cmd += '--i-phylogeny %s \\ \n' % trees[dat][1]
        cmd += '--p-metric %s \\ \n' % metric
        cmd += '--o-alpha-diversity %s\n' % out_fp
    else:
        cmd = 'qiime diversity alpha \\ \n'
        cmd += '--i-table %s \\ \n' % qza
        cmd += '--p-metric %s \\ \n' % metric
        cmd += '--o-alpha-diversity %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)
    return False


def write_metadata_tabulate(out_fp: str, divs: list, meta: str, cur_sh: TextIO) -> None:
    cmd = 'qiime metadata tabulate \\ \n'
    cmd += '--o-visualization %s \\ \n' % out_fp
    for div in divs:
        cmd += '--m-input-file %s \\ \n' % div
    cmd += '--m-input-file %s\n' % meta
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_alpha_group_significance_cmd(alpha: str, metadata: str, visu: str, cur_sh: TextIO) -> None:
    cmd = 'qiime diversity alpha-group-significance \\ \n'
    cmd += '--i-alpha-diversity %s \\ \n' % alpha
    cmd += '--m-metadata-file %s \\ \n' % metadata
    cmd += '--o-visualization %s\n' % visu
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def get_new_meta_pd(meta_pd: pd.DataFrame, case: str, case_var: str,
                    case_vals: list) -> pd.DataFrame:
    if 'ALL' in case:
        new_meta_pd = meta_pd.copy()
    elif len([x for x in case_vals if x[0] == '>' or x[0] == '<']):
        new_meta_pd = meta_pd.copy()
        for case_val in case_vals:
            if case_val[0] == '>':
                new_meta_pd = new_meta_pd[new_meta_pd[case_var] >= float(case_val[1:])].copy()
            elif case_val[0] == '<':
                new_meta_pd = new_meta_pd[new_meta_pd[case_var] <= float(case_val[1:])].copy()
    else:
        new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()
    return new_meta_pd


def get_new_alpha_div(case: str, div_qza: str, cur_rad: str,
                      new_meta_pd: pd.DataFrame, cur_sh: TextIO) -> str:
    new_div = '%s.qza' % cur_rad
    if 'ALL' in case:
        cur_sh.write('echo "cp %s %s"\n' % (div_qza, new_div))
        cur_sh.write('cp %s %s\n' % (div_qza, new_div))
    else:
        new_tsv = '%s.tsv' % cur_rad
        new_tsv_pd = pd.read_csv(div_qza.replace('.qza', '.tsv'), header=0, sep='\t')
        new_tsv_pd.rename(columns={new_tsv_pd.columns.tolist()[0]: 'Feature ID'}, inplace=True)
        new_tsv_pd.set_index('Feature ID', inplace=True)
        new_tsv_pd = new_tsv_pd.loc[new_meta_pd.index.tolist(), :]
        new_tsv_pd.reset_index().to_csv(new_tsv, index=False, sep='\t')
        cmd = run_import(new_tsv, new_div, 'SampleData[AlphaDiversity]')
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)
    return new_div
