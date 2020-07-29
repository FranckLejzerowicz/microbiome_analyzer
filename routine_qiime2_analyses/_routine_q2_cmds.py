# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import pandas as pd
from typing import TextIO
from os.path import dirname, isdir, isfile, splitext
from skbio.stats.ordination import OrdinationResults


def get_subset(tsv_pd: pd.DataFrame, subset_regex: list) -> list:
    """
    Make a feature metadata from the regex
    to get the names of the features to keep.

    :param tsv_pd: Table containing the features to subset.
    :param feats_subset: Feature metadata to create.
    :param subset_regex: subsetting regex.
    """
    to_keep_feats = {}
    for regex in subset_regex:
        to_keep_feats[regex.lower()] = tsv_pd.index.str.lower().str.contains(regex.lower())
    to_keep_feats_pd = pd.DataFrame(to_keep_feats)
    to_keep_feats = to_keep_feats_pd.any(axis=1)
    feats_subset_list = tsv_pd.index[to_keep_feats].tolist()
    return feats_subset_list


def write_filter_features(tsv_pd: pd.DataFrame, feats: list, qza: str,
                          qza_subset: str, meta_subset: str,
                          cur_sh: TextIO, dropout: bool) -> None:
    """
    filter-features: Filter features from table¶
    https://docs.qiime2.org/2020.2/plugins/available/feature-table/filter-features/

    Filter features from table based on frequency and/or metadata. Any samples
    with a frequency of zero after feature filtering will also be removed. See
    the filtering tutorial on https://docs.qiime2.org for additional details.

    :param qza: The features table to be filtered.
    :param qza_subset: The .
    :param meta_subset: Feature metadata to write.
    :param cur_sh: writing file handle.
    """

    if dropout:
        cmd = 'qiime feature-table filter-features \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--m-metadata-file %s \\\n' % meta_subset
        cmd += '--o-filtered-table %s\n' % qza_subset
    else:
        tsv_subset = '%s.tsv' % splitext(qza_subset)[0]
        tsv_nodrop = tsv_pd.loc[list(set(tsv_pd.index) & set(feats)), :].copy()
        tsv_nodrop.to_csv(tsv_subset, index=True, sep='\t')
        cmd = run_import(tsv_subset, qza_subset, "FeatureTable[Frequency]")
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_qemistree(feature_data: str, classyfire_qza: str, classyfire_tsv: str,
                    qemistree: str, cur_sh: TextIO) -> None:
    """
    :param feature_data:
    :param classyfire_qza:
    :param classyfire_tsv:
    :param qemistree:
    :param cur_sh: writing file handle.
    """
    cmd = ''
    if not isfile(classyfire_qza):
        cmd += 'qiime qemistree get-classyfire-taxonomy \\\n'
        cmd += '--i-feature-data %s \\\n' % feature_data
        cmd += '--o-classified-feature-data %s\n' % classyfire_qza
    classyfire_tsv = '%s.tsv' % splitext(classyfire_qza)[0]

    if not isfile(classyfire_tsv):
        cmd += run_export(classyfire_qza, classyfire_tsv, '')

    classyfire_level_qza = '%s-subclass.qza' % (splitext(classyfire_qza)[0])
    if not isfile(classyfire_level_qza):
        cmd += '\nqiime qemistree prune-hierarchy \\\n'
        cmd += '--i-feature-data %s \\\n' % classyfire_qza
        cmd += '--p-column subclass \\\n'
        cmd += '--i-tree %s \\\n' % qemistree
        cmd += '--o-pruned-tree %s\n' % classyfire_level_qza
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_barplots(out_qzv: str, qza: str, meta: str,
                   tax_qza: str, cur_sh: TextIO) -> None:
    """
    barplot: Visualize taxonomy with an interactive bar plot¶
    https://docs.qiime2.org/2020.2/plugins/available/taxa/barplot/

    This visualizer produces an interactive barplot visualization of
    taxonomies. Interactive features include multi-level sorting, plot
    recoloring, sample relabeling, and SVG figure export.

    :param cur_sh: writing file handle.
    """
    cmd = 'qiime taxa barplot \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--i-taxonomy %s \\\n' % tax_qza
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s \\\n' % out_qzv
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_taxonomy_sklearn(out_qza: str, out_fp_seqs_qza: str,
                           ref_classifier_qza: str) -> str:
    """
    Classify reads by taxon using a fitted classifier.
    https://docs.qiime2.org/2020.2/plugins/available/feature-classifier/classify-sklearn

    :param out_qza: Taxonomic classifications.
    :param out_fp_seqs_qza: The features sequences that should be classified.
    :param ref_classifier_qza: Taxonomic classifier.
    """
    cmd = 'qiime feature-classifier classify-sklearn \\\n'
    cmd += '--i-reads %s \\\n' % out_fp_seqs_qza
    cmd += '--i-classifier %s \\\n' % ref_classifier_qza
    cmd += '--p-n-jobs %s \\\n' % '4'
    cmd += '--o-classification %s\n' % out_qza
    return cmd


def write_rarefy(qza: str, qza_out: str, depth: str, cur_sh: TextIO) -> None:
    """
    Subsample frequencies from all samples so that the sum of frequencies in
    each sample is equal to sampling-depth.
    https://docs.qiime2.org/2019.10/plugins/available/feature-table/rarefy/

    :param qza: The feature table from which samples should be rarefied.
    :param qza_out: The resulting rarefied feature table in qza format.
    :param depth: The rarefaction depth.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime feature-table rarefy \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--p-sampling-depth %s \\\n' % depth
    cmd += '--o-rarefied-table %s\n' % qza_out
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_mmvec_cmd(meta_fp: str, qza1: str, qza2: str, res_dir: str,
                    odir: str, ranks_tsv: str, ordination_tsv: str,
                    batch: str, learn: str, epoch: str, prior: str,
                    thresh_feat: str, latent_dim: str,
                    train_column: str, n_example: str, gpu: bool,
                    standalone: bool, cur_sh: TextIO, qiime_env: str) -> None:
    """
    Performs bi-loglinear multinomial regression and calculates the
    conditional probability ranks of metabolite co-occurence given the microbe
    presence.

    :param meta_fp:
    :param qza1:
    :param qza2:
    :param res_dir:
    :param conditionals_tsv:
    :param ordination_tsv:
    :param batch:
    :param learn:
    :param epoch:
    :param prior:
    :param thresh_feat:
    :param latent_dim:
    :param train_column:
    :param n_example:
    :param gpu:
    :param standalone:
    :param cur_sh:
    :return:
    """

    cmd = ''
    if gpu or standalone:
        biom1 = '%s.biom' % splitext(qza1)[0]
        biom2 = '%s.biom' % splitext(qza2)[0]
        cmd += '\nmmvec paired-omics \\\n'
        if gpu:
            cmd += '--arm-the-gpu \\\n'
        cmd += '--microbe-file %s \\\n' % biom1
        cmd += '--metabolite-file %s \\\n' % biom2
        cmd += '--min-feature-count %s \\\n' % thresh_feat
        cmd += '--epochs %s \\\n' % epoch
        cmd += '--batch-size %s \\\n' % batch
        cmd += '--latent-dim %s \\\n' % latent_dim
        cmd += '--input-prior %s \\\n' % prior
        cmd += '--learning-rate %s \\\n' % learn
        cmd += '--beta1 0.85 \\\n'
        cmd += '--beta2 0.90 \\\n'
        cmd += '--checkpoint-interval 30 \\\n'
        cmd += '--summary-interval 30 \\\n'
        cmd += '--summary-dir %s \\\n' % res_dir
        cmd += '--ranks-file %s\n' % ranks_tsv
    else:
        ranks_qza = '%s.qza' % splitext(ranks_tsv)[0]
        ordination_qza = '%s.qza' % splitext(ordination_tsv)[0]
        if not isfile(ranks_qza) or not isfile(ordination_qza):
            cmd += 'current_time=$(date "+%y%m%d_%H%M%S")\n'
            cmd += '\ncd %s\n' % odir
            cmd += '\nqiime mmvec paired-omics \\\n'
            cmd += '--i-microbes %s \\\n' % qza1
            cmd += '--i-metabolites %s \\\n' % qza2
            cmd += '--m-metadata-file %s \\\n' % meta_fp
            if str(train_column) != 'None':
                cmd += '--p-training-column %s \\\n' % train_column
            cmd += '--p-num-testing-examples %s \\\n' % n_example
            cmd += '--p-min-feature-count %s \\\n' % thresh_feat
            cmd += '--p-epochs %s \\\n' % epoch
            cmd += '--p-batch-size %s \\\n' % batch
            cmd += '--p-latent-dim %s \\\n' % latent_dim
            cmd += '--p-input-prior %s \\\n' % prior
            cmd += '--p-learning-rate %s \\\n' % learn
            cmd += '--p-summary-interval 30 \\\n'
            if qiime_env == 'qiime2-2020.2':
                cmd += '--p-equalize-biplot \\\n'
            cmd += '--o-conditionals %s \\\n' % ranks_qza
            cmd += '--o-conditional-biplot %s \\\n' % ordination_qza
            cmd += '--output-dir %s/logdir_${current_time}\n' % dirname(ordination_qza)
        if not isfile(ranks_tsv):
            cmd += run_export(ranks_qza, ranks_tsv, '')
        if not isfile(ordination_tsv):
            cmd += run_export(ordination_qza, ordination_tsv, 'mmvec')
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def filter_feature_table(qza: str, new_qza: str, meta: str) -> str:
    """
    :param qza:
    :param new_qza:
    :param meta:
    :return:
    """
    cmd = '\nqiime feature-table filter-samples \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-filtered-table %s\n' % new_qza
    return cmd


def write_songbird_cmd(qza: str, new_qza: str, new_meta: str, formula: str,
                       epoch: str, batch: str, diff_prior: str, learn: str,
                       thresh_sample: str, thresh_feat: str, train_column: str,
                       metadatas: dict, diffs: str, diffs_qza: str, stats: str,
                       plot: str, base_diff_qza: str, base_stats: str, base_plot: str,
                       baseline_formula: str, tensor: str, tensor_html: str,
                       cur_sh: TextIO) -> None:
    """
    :param qza:
    :param new_qza:
    :param new_meta:
    :param formula:
    :param epoch:
    :param batch:
    :param diff_prior:
    :param learn:
    :param thresh_sample:
    :param thresh_feat:
    :param train:
    :param diffs:
    :param diffs_qza:
    :param stats:
    :param plot:
    :param base_diff_qza:
    :param base_stats:
    :param base_plot:
    :param tensor:
    :param cur_sh:
    """

    if new_qza not in metadatas:
        if not isfile(new_qza):
            cmd = filter_feature_table(qza, new_qza, new_meta)
            # cur_sh.write('echo "%s"\n' % cmd)
            cur_sh.write('%s\n' % cmd)
        metadatas[new_qza] = new_meta

        # if not isfile(diffs_qza):
        cmd = '\nqiime songbird multinomial \\\n'
        cmd += ' --i-table %s \\\n' % new_qza
        cmd += ' --m-metadata-file %s \\\n' % new_meta
        cmd += ' --p-formula "%s" \\\n' % formula
        cmd += ' --p-epochs %s \\\n' % epoch
        cmd += ' --p-batch-size %s \\\n' % batch
        cmd += ' --p-differential-prior %s \\\n' % diff_prior
        cmd += ' --p-learning-rate %s \\\n' % learn
        cmd += ' --p-min-sample-count %s \\\n' % thresh_sample
        cmd += ' --p-min-feature-count %s \\\n' % thresh_feat
        cmd += ' --p-training-column %s \\\n' % train_column
        cmd += ' --p-summary-interval 2 \\\n'
        cmd += ' --o-differentials %s \\\n' % diffs_qza
        cmd += ' --o-regression-stats %s \\\n' % stats
        cmd += ' --o-regression-biplot %s\n' % plot
        cur_sh.write('%s\n' % cmd)

    if not isfile(diffs):
        cmd = run_export(diffs_qza, diffs, '')
        # cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if len(base_diff_qza) and not isfile(base_diff_qza):
        cmd = '\nqiime songbird multinomial \\\n'
        cmd += ' --i-table %s \\\n' % new_qza
        cmd += ' --m-metadata-file %s \\\n' % new_meta
        cmd += ' --p-formula "%s" \\\n' % baseline_formula
        cmd += ' --p-epochs %s \\\n' % epoch
        cmd += ' --p-batch-size %s \\\n' % batch
        cmd += ' --p-differential-prior %s \\\n' % diff_prior
        cmd += ' --p-learning-rate %s \\\n' % learn
        cmd += ' --p-min-sample-count %s \\\n' % thresh_sample
        cmd += ' --p-min-feature-count %s \\\n' % thresh_feat
        cmd += ' --p-training-column %s \\\n' % train_column
        cmd += ' --p-summary-interval 2 \\\n'
        cmd += ' --o-differentials %s \\\n' % base_diff_qza
        cmd += ' --o-regression-stats %s \\\n' % base_stats
        cmd += ' --o-regression-biplot %s\n' % base_plot
        cur_sh.write('%s\n' % cmd)

    if not isfile(tensor):
        cmd = '\n\nqiime songbird summarize-paired \\\n'
        cmd += ' --i-regression-stats %s \\\n' % stats
        cmd += ' --i-baseline-stats %s \\\n' % base_stats
        cmd += ' --o-visualization %s\n' % tensor
        # cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if not isdir(tensor_html):
        cmd = run_export(tensor, tensor_html, 'songbird')
        # cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)


def write_phate_cmd(tsv: str,
                    cur_sh: TextIO) -> None:
    """
    :param qza:
    :param new_qza:
    :param new_meta:
    :param formula:
    :param epoch:
    :param batch:
    :param diff_prior:
    :param learn:
    :param thresh_sample:
    :param thresh_feat:
    :param train:
    :param diffs:
    :param diffs_qza:
    :param stats:
    :param plot:
    :param base_diff_qza:
    :param base_stats:
    :param base_plot:
    :param tensor:
    :param cur_sh:
    """
    pass


def run_import(input_path: str, output_path: str, typ: str) -> str:
    """
    Return the import qiime2 command.

    :param input_path: input file path.
    :param output_path: output file path.
    :param typ: qiime2 type.
    :return: command to qiime2.
    """
    cmd = ''
    if typ.startswith("FeatureTable"):
        if not input_path.endswith('biom'):
            cur_biom = '%s.biom' % splitext(input_path)[0]
            cmd += 'biom convert \\\n'
            cmd += '  -i %s \\\n' % input_path
            cmd += '  -o %s \\\n' % cur_biom
            cmd += '  --table-type="OTU table" \\\n'
            cmd += '  --to-hdf5\n\n'
            cmd += 'qiime tools import \\\n'
            cmd += '  --input-path %s \\\n' % cur_biom
            cmd += '  --output-path %s \\\n' % output_path
            cmd += '  --type "FeatureTable[Frequency]"\n'
        else:
            cmd += 'qiime tools import \\\n'
            cmd += '  --input-path %s \\\n' % input_path
            cmd += '  --output-path %s \\\n' % output_path
            cmd += '  --type "FeatureTable[Frequency]"\n'
    else:
        cmd += 'qiime tools import \\\n'
        cmd += '  --input-path %s \\\n' % input_path
        cmd += '  --output-path %s \\\n' % output_path
        cmd += '  --type "%s"\n' % typ
    return cmd


def run_export(input_path: str, output_path: str, typ: str) -> str:
    """
    Return the export qiime2 command.

    :param input_path: input file path.
    :param output_path: output file path.
    :param typ: qiime2 type.
    :return: command to qiime2.
    """
    cmd = ''
    if typ.startswith("FeatureTable"):
        if not output_path.endswith('biom'):
            cur_biom = '%s.biom' % splitext(output_path)[0]
            cmd += 'qiime tools export \\\n'
            cmd += '  --input-path %s \\\n' % input_path
            cmd += '  --output-path %s\n' % splitext(output_path)[0]
            cmd += 'mv %s/*.biom %s\n' % (splitext(output_path)[0], cur_biom)
            cmd += 'biom convert'
            cmd += '  -i %s \\\n' % cur_biom
            cmd += '  -o %s.tmp \\\n' % output_path
            cmd += '  --to-tsv\n\n'
            cmd += 'tail -n +2 %s.tmp > %s\n\n' % (output_path, output_path)
            cmd += 'rm -rf %s %s.tmp\n' % (splitext(output_path)[0], output_path)
        else:
            cmd += 'qiime tools export \\\n'
            cmd += '  --input-path %s \\\n' % input_path
            cmd += '  --output-path %s\n' % splitext(output_path)[0]
            cmd += 'mv %s/*.biom %s\n' % (splitext(input_path)[0], output_path)
            cmd += 'rm -rf %s\n' % splitext(input_path)[0]
    else:
        cmd += 'qiime tools export \\\n'
        cmd += '  --input-path %s \\\n' % input_path
        cmd += '  --output-path %s\n' % splitext(output_path)[0]
        if 'Phylogeny' in typ:
            cmd += 'mv %s/*.nwk %s\n' % (splitext(output_path)[0], output_path)
        elif 'biplot' in typ:
            cmd += 'mv %s/*.txt %s\n' % (splitext(output_path)[0], output_path)
        elif 'mmvec' in typ:
            cmd += 'mv %s/*.txt %s\n' % (splitext(output_path)[0], output_path)
        elif 'pcoa' in typ:
            cmd += 'mv %s/*.txt %s\n' % (splitext(output_path)[0], output_path)
        elif 'perms' in typ:
            cmd += 'mv %s/index.html %s\n' % (splitext(output_path)[0], output_path)
        elif 'songbird' in typ:
            cmd += 'mv %s/index.html %s\n' % (splitext(output_path)[0], output_path)
        else:
            cmd += 'mv %s/*.tsv %s\n' % (splitext(output_path)[0], output_path)
        cmd += 'rm -rf %s\n' % splitext(output_path)[0]
    return cmd


def write_diversity_beta(out_fp: str, datasets_phylo: dict, trees: dict, dat: str,
                         qza: str, metric: str, cur_sh: TextIO, subset: bool) -> None:
    """
    Computes a user-specified beta diversity metric for all pairs of samples
    in a feature table.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta/
    and
    Computes a user-specified phylogenetic beta diversity metric for all pairs
    of samples in a feature table.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-phylogenetic/

    :param out_fp: The resulting distance matrix.
    :param datasets_phylo: pholygenetic decision.
    :param trees: pholygenetic trees.
    :param dat: current dataset.
    :param qza: The feature table containing the samples over which beta diversity should be computed.
    :param metric: The beta diversity metric to be computed.
    :param cur_sh: writing file handle.
    :return: whether the command is to be skipped or not.
    """
    if 'unifrac' in metric:
        # if not datasets_phylo[dat][0] or dat not in trees:
        #     return True
        cmd = 'qiime diversity beta-phylogenetic \\\n'
        if subset:
            cmd += '--i-table %s \\\n' % qza
        else:
            if datasets_phylo[dat][1]:
                cmd += '--i-table %s \\\n' % trees[dat][0]
            else:
                cmd += '--i-table %s \\\n' % qza
        cmd += '--i-phylogeny %s \\\n' % trees[dat][1]
    else:
        cmd = 'qiime diversity beta \\\n'
        cmd += '--i-table %s \\\n' % qza
    cmd += '--p-metric %s \\\n' % metric
    cmd += '--p-n-jobs 1 \\\n'
    cmd += '--o-distance-matrix %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)
    # return False


def write_diversity_pcoa(DM: str, out_pcoa: str, out_tsv: str, cur_sh: TextIO) -> None:
    """
    Apply principal coordinate analysis.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/

    :param DM: The distance matrix on which PCoA should be computed.
    :param out_pcoa: The resulting PCoA matrix.
    :param cur_sh: writing file handle.
    """
    if not isfile(out_pcoa):
        cmd = 'qiime diversity pcoa \\\n'
        cmd += '--i-distance-matrix %s \\\n' % DM
        cmd += '--o-pcoa %s\n' % out_pcoa
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n\n' % cmd)
    if not isfile(out_tsv):
        cmd = run_export(out_pcoa, out_tsv, 'pcoa')
        cur_sh.write('%s\n\n' % cmd)


def write_diversity_biplot(tsv: str, qza: str, out_pcoa: str,
                           out_biplot: str, tax_qza: str,
                           tsv_tax: str, cur_sh: TextIO) -> None:
    """
    pcoa-biplot: Principal Coordinate Analysis Biplot.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa-biplot/

    :param DM: The distance matrix on which PCoA should be computed.
    :param out_pcoa: The resulting PCoA matrix.
    :param cur_sh: writing file handle.
    """
    cmd = ''
    if tax_qza != 'missing':
        biplot_tab_tsv = '%s_table.tsv' % splitext(out_biplot)[0]
        biplot_tab_qza = '%s.qza' % splitext(biplot_tab_tsv)[0]
        tax_dict = {}
        with open('%s.tsv' % splitext(tax_qza)[0]) as f, open(tsv_tax, 'w') as o_tax:
            o_tax.write('Feature ID\tTaxon\tPrevious ID\n')
            n = 0
            for ldx, line in enumerate(f):
                if ldx and not line.startswith('#q2:types'):
                    new = 'x__%s;%s' % (n, line.strip().split('\t')[1])
                    tax_dict[line.split('\t')[0]] = new
                    o_tax.write('%s\t%s\t%s\n' % (new, new, line.split('\t')[0]))
                    n += 1
        with open(tsv) as f, open(biplot_tab_tsv, 'w') as o_tab:
            for ldx, line in enumerate(f):
                t = line.strip().split('\t')
                if t[0] in tax_dict:
                    o_tab.write('%s\t%s\n' % (tax_dict[t[0]], '\t'.join(t[1:])))
                else:
                    o_tab.write(line)
        cmd += run_import(biplot_tab_tsv, biplot_tab_qza, 'FeatureTable[Frequency]')
    else:
        biplot_tab_qza = qza

    biplot_tab_rel_qza = '%s_rel.qza' % splitext(biplot_tab_qza)[0]

    cmd += '\nqiime feature-table relative-frequency \\\n'
    cmd += '--i-table %s \\\n' % biplot_tab_qza
    cmd += '--o-relative-frequency-table %s\n' % biplot_tab_rel_qza

    cmd += 'qiime diversity pcoa-biplot \\\n'
    cmd += '--i-pcoa %s \\\n' % out_pcoa
    cmd += '--i-features %s \\\n' % biplot_tab_rel_qza
    cmd += '--o-biplot %s\n' % out_biplot

    cmd += 'rm %s\n' % biplot_tab_rel_qza
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)

    out_biplot_txt = '%s.txt' % splitext(out_biplot)[0]
    cmd = run_export(out_biplot, out_biplot_txt, 'biplot')
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_emperor_biplot(meta: str, biplot: str, out_plot: str,
                         cur_sh: TextIO, taxonomy: str) -> None:
    """
    Generates an interactive ordination plot where the user can visually
    integrate sample metadata.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param meta: The sample metadata.
    :param pcoa: The principal coordinates matrix to be plotted.
    :param out_plot: VISUALIZATION.
    :param cur_sh: writing file handle.
    """
    biplot_txt = '%s.txt' % splitext(biplot)[0]
    if isfile(biplot_txt):
        ordi = OrdinationResults.read(biplot_txt)
        ordi.features = ordi.features.iloc[:,:3]
        ordi.samples= ordi.samples.iloc[:,:3]
        ordi.eigvals = ordi.eigvals[:3]
        ordi.proportion_explained = ordi.proportion_explained[:3]
        ordi.write(biplot_txt)
    cmd = run_import(biplot_txt, biplot, "PCoAResults % Properties('biplot')")
    cur_sh.write('%s\n\n' % cmd)

    cmd = 'qiime emperor biplot \\\n'
    cmd += '--i-biplot %s \\\n' % biplot
    cmd += '--m-sample-metadata-file %s \\\n' % meta
    if taxonomy != 'missing':
        tax_tmp = '%s_taxonomy.tmp' % splitext(biplot)[0]
        tax_pd = pd.read_csv(taxonomy, header=0, sep='\t')
        if 'Taxon' in tax_pd:
            tax_split_pd = get_split_taxonomy(tax_pd.Taxon.tolist())
            tax_pd = pd.concat([tax_pd, tax_split_pd], axis=1, sort=False)
            tax_pd.to_csv(tax_tmp, index=False, sep='\t')
        else:
            tax_tmp = taxonomy
        cmd += '--m-feature-metadata-file %s \\\n' % tax_tmp
    cmd += '--p-number-of-features 20 \\\n'
    cmd += '--o-visualization %s\n' % out_plot
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_emperor(meta: str, pcoa: str, out_plot: str, cur_sh: TextIO) -> None:
    """
    Generates an interactive ordination plot where the user can visually
    integrate sample metadata.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param meta: The sample metadata.
    :param pcoa: The principal coordinates matrix to be plotted.
    :param out_plot: VISUALIZATION.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime emperor plot \\\n'
    cmd += '--i-pcoa %s \\\n' % pcoa
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s\n' % out_plot
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_seqs_fasta(out_fp_seqs_fasta: str, out_fp_seqs_qza: str,
                     tsv_pd: pd.DataFrame) -> str:
    """
    Write the fasta sequences.

    :param out_fp_seqs_fasta: output sequences fasta file name.
    :param out_fp_seqs_qza: output sequences qiime2 Artefact file name.
    :param tsv_pd: table which feature names are sequences.
    :param cur_sh: writing file handle.
    """
    with open(out_fp_seqs_fasta, 'w') as fas_o:
        for seq in tsv_pd.index:
            fas_o.write('>%s\n%s\n' % (seq.strip(), seq.strip()))
    cmd = run_import(out_fp_seqs_fasta, out_fp_seqs_qza, 'FeatureData[Sequence]')
    return cmd


def write_fragment_insertion(out_fp_seqs_qza: str, ref_tree_qza: str,
                             out_fp_sepp_tree: str, out_fp_sepp_plac: str,
                             qza: str, qza_in: str, qza_out: str,
                             cur_sh: TextIO) -> None:
    """
    Perform fragment insertion of sequences using the SEPP algorithm.
    https://docs.qiime2.org/2019.10/plugins/available/fragment-insertion/sepp/

    and

    Filters fragments not inserted into a phylogenetic tree from a feature-
    table. Some fragments computed by e.g. Deblur or DADA2 are too remote to
    get inserted by SEPP into a reference phylogeny. To be able to use the
    feature-table for downstream analyses like computing Faith's PD or
    UniFrac, the feature-table must be cleared of fragments that are not part
    of the phylogenetic tree, because their path length can otherwise not be
    determined. Typically, the number of rejected fragments is low (<= 10),
    but it might be worth to inspect the ratio of reads assigned to those
    rejected fragments.
    https://docs.qiime2.org/2019.10/plugins/available/fragment-insertion/filter-features/

    :param out_fp_seqs_qza: The sequences to insert into the reference tree.
    :param ref_tree_qza: The reference database to insert the representative sequences into.
    :param out_fp_sepp_tree: The tree with inserted feature data.
    :param out_fp_sepp_plac: Information about the feature placements within the reference tree.
    :param qza: A feature-table which needs to filtered down to those fragments
                       that are contained in the tree, e.g. result of a Deblur or DADA2 run.
    :param qza_in: The input table minus those fragments that were not
                       part of the tree. This feature-table can be used for
                       downstream analyses like phylogenetic alpha- or beta-
                       diversity computation.
    :param qza_out: Those fragments that got removed from the input table,
                       because they were not part of the tree. This table is
                       mainly used for quality control, e.g. to inspect the
                       ratio of removed reads per sample from the input table.
                       You can ignore this table for downstream analyses.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime fragment-insertion sepp \\\n'
    cmd += '--i-representative-sequences %s \\\n' % out_fp_seqs_qza
    cmd += '--i-reference-database %s \\\n' % ref_tree_qza
    cmd += '--o-tree %s \\\n' % out_fp_sepp_tree
    cmd += '--o-placements %s \\\n' % out_fp_sepp_plac
    cmd += '--p-threads 24\n'
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)
    cmd = 'qiime fragment-insertion filter-features \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--i-tree %s \\\n' % out_fp_sepp_tree
    cmd += '--o-filtered-table %s \\\n' % qza_in
    cmd += '--o-removed-table %s\n' % qza_out
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)
    qza_in_tsv = '%s.tsv' % splitext(qza_in)[0]
    cmd = run_export(qza_in, qza_in_tsv, 'FeatureTable')
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_doc(qza: str, fp: str, fa: str, new_meta: str, new_qza: str,
              new_tsv: str, cur_rad: str, new_tsv_token: str,
              cur_rad_token: str, n_nodes: str, n_procs: str,
              doc_params: dict, cur_sh: TextIO) -> None:

    if not isfile(new_qza):
        cmd = '\nqiime feature-table filter-samples \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--m-metadata-file %s \\\n' % new_meta
        cmd += '--o-filtered-table %s\n' % new_qza
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

        cmd += run_export(new_qza, new_tsv, 'FeatureTable')
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    cmd = 'mkdir -p %s\n' % cur_rad_token
    cmd += 'mkdir -p %s\n' % cur_rad
    cmd += 'cp %s %s\n' % (new_tsv, new_tsv_token)
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)

    cmd = '\nXDOC \\\n'
    cmd += '--i-otu %s \\\n' % new_tsv_token
    cmd += '--o-outdir %s \\\n' % cur_rad_token
    cmd += '-fp %s \\\n' % fp
    cmd += '-fa %s \\\n' % fa
    cmd += '--p-cpus %s \\\n' % (int(n_nodes) * int(n_procs))
    cmd += '--p-r %s \\\n' % doc_params['r']
    cmd += '--p-subr %s \\\n' % doc_params['subr']
    cmd += '--p-mov-avg %s \\\n' % doc_params['mov_avg']
    # for ci in doc_params['ci']:
    #     cmd += '--p-ci %s \\\n' % ci
    cmd += '--p-span %s \\\n' % doc_params['span']
    cmd += '--p-degree %s \\\n' % doc_params['degree']
    cmd += '--p-family %s \\\n' % doc_params['family']
    cmd += '--p-iterations %s \\\n' % doc_params['iterations']
    cmd += '--p-surface %s \\\n' % doc_params['surface']
    cmd += '--p-nulls %s \\\n' % doc_params['nulls']
    if doc_params['non_zero']:
        cmd += '--non-zero \\\n'
    if doc_params['null']:
        cmd += '--null \\\n'
    cmd += '--verbose\n\n'
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)

    cmd = 'rsync -am %s %s\n\n' % (cur_rad_token, cur_rad)
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)

    cmd = 'rm -rf %s\n\n' % cur_rad_token
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def write_deicode_biplot(qza: str, new_meta: str, new_qza: str, ordi_qza: str,
                         new_mat_qza: str, ordi_qzv: str, cur_sh: TextIO) -> None:
    """
    Performs robust center log-ratio transform robust PCA and
    ranks the features by the loadings of the resulting SVD.
    https://library.qiime2.org/plugins/deicode/19/

    :param qza: The feature table from which samples should be filtered.
    :param new_meta: Sample metadata containing formula terms.
    :param new_qza: The resulting feature table filtered by sample & Input table of counts.
    :param ordi_qza: A biplot of the (Robust Aitchison) RPCA feature loadings.
    :param new_mat_qza: The Aitchison distance ofthe sample loadings from RPCA.
    :param ordi_qzv: VISUALIZATION
    :param cur_sh: writing file handle.
    """
    cmd = '\n'
    if not isfile(new_qza):
        cmd += '\nqiime feature-table filter-samples \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--m-metadata-file %s \\\n' % new_meta
        cmd += '--o-filtered-table %s\n' % new_qza
    if not isfile(ordi_qza) or not isfile(new_mat_qza):
        cmd += 'qiime deicode rpca \\\n'
        cmd += '--i-table %s \\\n' % new_qza
        cmd += '--p-min-feature-count 10 \\\n'
        cmd += '--p-min-sample-count 500 \\\n'
        cmd += '--o-biplot %s \\\n' % ordi_qza
        cmd += '--o-distance-matrix %s\n' % new_mat_qza
    cmd += 'qiime emperor biplot \\\n'
    cmd += '--i-biplot %s \\\n' % ordi_qza
    cmd += '--m-sample-metadata-file %s \\\n' % new_meta
    cmd += '--o-visualization %s \\\n' % ordi_qzv
    cmd += '--p-number-of-features 20\n'
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def add_q2_types_to_meta(new_meta_pd: pd.DataFrame, new_meta: str) -> None:
    """
    Merge the q2-types to the metadata for PERMANOVA.

    :param new_meta_pd: metadata table.
    :param new_meta: metadata table file name.
    """
    q2types = pd.DataFrame(
        [(['#q2:types'] + ['categorical'] * new_meta_pd.shape[1])],
        columns=new_meta_pd.reset_index().columns.tolist(),
    )
    q2types.rename(columns={q2types.columns.tolist()[0]: new_meta_pd.index.name}, inplace=True)
    q2types.set_index(new_meta_pd.index.name, inplace=True)
    new_meta_pd = pd.concat([q2types, new_meta_pd])
    new_meta_pd.reset_index().to_csv(new_meta, index=False, sep='\t')


def check_absence_mat(mat_qzas: list, first_print: int, analysis: str) -> bool:
    """
    Check the absence of the beta diversity matrix.

    :param mat_qzas: beta diveristy matrices files names.
    :param first_print: whether to print the first message or not.
    :param analysis: Current analysis.
    :return: Whether the beta diversity matrix is absent or not.
    """
    presence_mat = [mat_qza for mat_qza in mat_qzas if isfile(mat_qza)]
    if not presence_mat:
        if not first_print:
            print('Beta diversity, distances matrices must be generated already to automatise %s\n'
                  '\t(re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)' % analysis)
            first_print += 1
        return True
    return False


def write_diversity_beta_group_significance(new_meta: str, mat_qza: str, new_mat_qza: str,
                                            testing_group: str, beta_type: str, new_qzv: str,
                                            new_html: str, cur_sh: TextIO) -> None:
    """
    Determine whether groups of samples are significantly different from one
    another using a permutation-based statistical test.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-group-significance/

    Includes calls to:
    filter-distance-matrix: Filter samples from a distance matrix
    https://docs.qiime2.org/2019.10/plugins/available/diversity/filter-distance-matrix/
    filter-samples: Filter samples from table
    https://docs.qiime2.org/2019.10/plugins/available/feature-table/filter-samples/

    :param new_meta: Sample metadata containing formula terms.
    :param mat_qza: Distance matrix to filter by sample.
    :param new_mat_qza: Matrix of distances between pairs of samples.
    :param qza: The feature table from which samples should be filtered.
    :param new_qza: The resulting feature table filtered by sample.
    :param testing_group: Categorical sample metadata column.
    :param new_qzv: VISUALIZATION.
    :param cur_sh: writing file handle.
    """
    if not isfile(new_qzv):
        cmd = 'qiime diversity filter-distance-matrix \\\n'
        cmd += '--m-metadata-file %s \\\n' % new_meta
        cmd += '--i-distance-matrix %s \\\n' % mat_qza
        cmd += '--o-filtered-distance-matrix %s\n' % new_mat_qza
        cmd += 'qiime diversity beta-group-significance \\\n'
        cmd += '--i-distance-matrix %s \\\n' % new_mat_qza
        cmd += '--p-method %s \\\n' % beta_type
        cmd += '--m-metadata-file %s \\\n' % new_meta
        cmd += '--m-metadata-column "%s" \\\n' % testing_group
        cmd += '--p-permutations 2999 \\\n'
        cmd += '--o-visualization %s\n' % new_qzv
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write(cmd)
    if not isfile(new_html):
        cmd = run_export(new_qzv, new_html, 'perms')
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write(cmd)


def write_diversity_adonis(new_meta: str, mat_qza: str, new_mat_qza: str,
                           formula: str, new_qzv: str, cur_sh: TextIO) -> None:
    """
    Determine whether groups of samples are significantly different from one
    another using the ADONIS permutation-based statistical test in vegan-R.
    The function partitions sums of squares of a multivariate data set, and is
    directly analogous to MANOVA (multivariate analysis of variance). This
    action differs from beta_group_significance in that it accepts R formulae
    to perform multi-way ADONIS tests; beta_group_signficance only performs
    one-way tests. For more details see
    http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html
    https://docs.qiime2.org/2019.10/plugins/available/diversity/adonis/

    Includes calls to:
    * filter-distance-matrix: Filter samples from a distance matrix
    https://docs.qiime2.org/2019.10/plugins/available/diversity/filter-distance-matrix/
    * filter-samples: Filter samples from table
    https://docs.qiime2.org/2019.10/plugins/available/feature-table/filter-samples/

    :param new_meta: Sample metadata containing formula terms.
    :param mat_qza: Distance matrix to filter by sample.
    :param new_mat_qza: Matrix of distances between pairs of samples.
    :param qza: The feature table from which samples should be filtered.
    :param new_qza: The resulting feature table filtered by sample.
    :param formula: Model formula containing only independent terms
                       contained in the sample metadata. These can be
                       continuous variables or factors, and they can have
                       interactions as in a typical R formula. E.g., the
                       formula "treatment+block" would test whether the input
                       distance matrix partitions based on "treatment" and
                       "block" sample metadata. The formula "treatment*block"
                       would test both of those effects as well as their
                       interaction. Enclose formulae in quotes to avoid
                       unpleasant surprises.
    :param new_qzv: VISUALIZATION.
    :param cur_sh: writing file handle.
    """
    cmd = '\nqiime diversity filter-distance-matrix \\\n'
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--i-distance-matrix %s \\\n' % mat_qza
    cmd += '--o-filtered-distance-matrix %s\n' % new_mat_qza
    cmd += 'qiime diversity adonis \\\n'
    cmd += '--i-distance-matrix %s \\\n' % new_mat_qza
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--p-formula "%s" \\\n' % formula
    cmd += '--p-permutations 2999 \\\n'
    cmd += '--p-n-jobs 6 \\\n'
    cmd += '--o-visualization %s\n' % new_qzv
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def write_procrustes(common_meta_fp: str, dm1: str, dm2: str,
                     dm_out1: str, dm_out2: str, biplot: str,
                     cur_sh: TextIO) -> None:

    dm_out_tsv1 = '%s.tsv' % splitext(dm_out1)[0]
    dm_out_tsv2 = '%s.tsv' % splitext(dm_out2)[0]
    pcoa_out1 = '%s_PCoA.qza' % splitext(dm_out1)[0]
    pcoa_out2 = '%s_PCoA.qza' % splitext(dm_out2)[0]
    ref_pcoa = '%s_ref.qza' % splitext(pcoa_out1)[0]
    oth_pcoa = '%s_oth.qza' % splitext(pcoa_out2)[0]

    cmd = ''
    if not isfile(dm_out1) or not isfile(dm_out_tsv1):
        cmd += '\nqiime diversity filter-distance-matrix \\\n'
        cmd += '--m-metadata-file %s \\\n' % common_meta_fp
        cmd += '--i-distance-matrix %s \\\n' % dm1
        cmd += '--o-filtered-distance-matrix %s\n' % dm_out1
        cmd += run_export(dm_out1, dm_out_tsv1, '')
    if not isfile(dm_out2) or not isfile(dm_out_tsv2):
        cmd += '\nqiime diversity filter-distance-matrix \\\n'
        cmd += '--m-metadata-file %s \\\n' % common_meta_fp
        cmd += '--i-distance-matrix %s \\\n' % dm2
        cmd += '--o-filtered-distance-matrix %s\n' % dm_out2
        cmd += run_export(dm_out2, dm_out_tsv2, '')
    if not isfile(pcoa_out1):
        cmd += '\nqiime diversity pcoa \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm_out1
        cmd += '--o-pcoa %s\n' % pcoa_out1
    if not isfile(pcoa_out2):
        cmd += '\nqiime diversity pcoa \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm_out2
        cmd += '--o-pcoa %s\n' % pcoa_out2
    if not isfile(ref_pcoa) or not isfile(oth_pcoa):
        cmd += '\nqiime diversity procrustes-analysis \\\n'
        cmd += '--i-reference %s \\\n' % pcoa_out1
        cmd += '--i-other %s \\\n' % pcoa_out2
        cmd += '--o-transformed-reference %s \\\n' % ref_pcoa
        cmd += '--o-transformed-other %s\n' % oth_pcoa
    if not isfile(biplot):
        cmd += '\nqiime emperor procrustes-plot \\\n'
        cmd += '--i-reference-pcoa %s \\\n' % ref_pcoa
        cmd += '--i-other-pcoa %s \\\n' % oth_pcoa
        cmd += '--m-metadata-file %s \\\n' % common_meta_fp
        cmd += '--o-visualization %s\n' % biplot
    if isfile(common_meta_fp):
        cmd += 'rm %s\n' % common_meta_fp
    if isfile(pcoa_out1):
        cmd += 'rm %s\n' % pcoa_out1
    if isfile(pcoa_out2):
        cmd += 'rm %s\n' % pcoa_out2
    if isfile(dm_out1):
        cmd += 'rm %s\n' % dm_out1
    if isfile(dm_out2):
        cmd += 'rm %s\n' % dm_out2
    if isfile(ref_pcoa):
        cmd += 'rm %s\n' % ref_pcoa
    if isfile(oth_pcoa):
        cmd += 'rm %s\n' % oth_pcoa
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n' % cmd)


def get_metric(metrics: list, file_name: str) -> str:
    """
    Get the current diversity from the file name.
    :param metrics: all the diversity metrics.
    :param file_name: file name.
    :return: either diversity metric or nothing.
    """
    for metric in metrics:
        if metric in file_name:
            return metric
    return ''


def get_case(case_vals: list, case_var: str, form: str = None) -> str:
    """
    Get the current case, which is the concatenation of:
     - diversity metric.
     - metadata variable.
     - metadata variable's values.
     - formula.
    :param case_vals: variable's values
    :param case_var: metadata variable.
    :param form: formula.
    :return: current case.
    """
    if len(case_vals):
        case = '%s_%s' % (case_var, '-'.join(
            [x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
    else:
        case = case_var
    if form:
        case = '%s_%s' % (case, form)
    case = case.replace('__', '_').replace(' ', '-').replace('/', '').replace('(', '').replace(')', '')
    return case


def write_longitudinal_volatility(out_fp: str, meta_alphas: str,
                                  time_point: str, cur_sh: TextIO) -> None:
    """
    Generate an interactive control chart depicting the longitudinal
    volatility of sample metadata and/or feature frequencies across time (as
    set using the "state_column" parameter). Any numeric metadata column (and
    metadata-transformable artifacts, e.g., alpha diversity results) can be
    plotted on the y-axis, and are selectable using the "metric_column"
    selector. Metric values are averaged to compare across any categorical
    metadata column using the "group_column" selector. Longitudinal volatility
    for individual subjects sampled over time is co-plotted as "spaghetti"
    plots if the "individual_id_column" parameter is used. state_column will
    typically be a measure of time, but any numeric metadata column can be
    used.
    https://docs.qiime2.org/2019.10/plugins/available/longitudinal/volatility/

    :param out_fp: VISUALIZATION.
    :param meta_alphas: Sample metadata file containing arguments will be individual-id-column.
    :param time_point: Metadata column containing state (time) variable information.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime longitudinal volatility \\\n'
    cmd += '--m-metadata-file %s \\\n' % meta_alphas
    cmd += '--p-state-column "%s" \\\n' % time_point
    cmd += '--p-individual-id-column "host_subject_id"'
    cmd += '--o-visualization %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_diversity_alpha_correlation(out_fp: str, qza: str, method: str,
                                      meta: str, cur_sh: TextIO) -> None:
    """
    Determine whether numeric sample metadata columns are correlated with alpha diversity.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-correlation/

    :param out_fp: Vector containing per-sample alpha diversities.
    :param qza: Vector of alpha diversity values by sample.
    :param method: The correlation test to be applied.
    :param meta: The sample metadata.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime diversity alpha-correlation \\\n'
    cmd += '--i-alpha-diversity %s \\\n' % qza
    cmd += '--p-method %s \\\n' % method
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_diversity_alpha(out_fp: str, datasets_phylo: dict, trees: dict, dat: str,
                          qza: str, metric: str, cur_sh: TextIO) -> bool:
    """
    Computes a user-specified alpha diversity metric for all samples in a feature table.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha/

    :param out_fp: Vector containing per-sample alpha diversities.
    :param datasets_phylo: phylogenetic decision.
    :param trees: The feature table containing the samples for which alpha diversity should be computed + tree.
    :param dat: dataset.
    :param qza: The feature table containing the samples for which alpha diversity should be computed.
    :param metric: The alpha diversity metric to be computed.
    :param cur_sh: writing file handle.
    :return: whether the command is to be skipped or not.
    """
    if metric in ['faith_pd']:
        if not datasets_phylo[dat][0] or dat not in trees:
            return True
        cmd = 'qiime diversity alpha-phylogenetic \\\n'
        if datasets_phylo[dat][1]:
            cmd += '--i-table %s \\\n' % trees[dat][0]
        else:
            cmd += '--i-table %s \\\n' % qza
        cmd += '--i-phylogeny %s \\\n' % trees[dat][1]
        cmd += '--p-metric %s \\\n' % metric
        cmd += '--o-alpha-diversity %s\n' % out_fp
    else:
        cmd = 'qiime diversity alpha \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--p-metric %s \\\n' % metric
        cmd += '--o-alpha-diversity %s\n' % out_fp
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)
    return False


def write_metadata_tabulate(out_fp: str, divs: list, meta: str, cur_sh: TextIO) -> None:
    """
    Generate a tabular view of Metadata. The output visualization supports
    interactive filtering, sorting, and exporting to common file formats.
    https://docs.qiime2.org/2019.10/plugins/available/metadata/tabulate/

    :param out_fp: VISUALIZATION.
    :param divs: The alpha diversity vectors to tabulate.
    :param meta: The metadata to tabulate.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime metadata tabulate \\\n'
    cmd += '--o-visualization %s \\\n' % out_fp
    for div in divs:
        cmd += '--m-input-file %s \\\n' % div
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_alpha_group_significance_cmd(alpha: str, metadata: str, visu: str, cur_sh: TextIO) -> None:
    """
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha-group-significance/

    :param alpha: Vector of alpha diversity values by sample.
    :param metadata: The sample metadata.
    :param visu: VISUALIZATION.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime diversity alpha-group-significance \\\n'
    cmd += '--i-alpha-diversity %s \\\n' % alpha
    cmd += '--m-metadata-file %s \\\n' % metadata
    cmd += '--o-visualization %s\n' % visu
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def get_new_meta_pd(meta_pd: pd.DataFrame, case: str,
                    case_var: str, case_vals: list) -> pd.DataFrame:
    """
    Perform subset.

    :param meta_pd: metadata table.
    :param case: concatenation of the current subset / test groups.
    :param case_var: current variable for subset.
    :param case_vals: current variable's values for subset.
    :return: Subsetted metadata table.
    """
    if 'ALL' in case:
        new_meta_pd = meta_pd.copy()
    elif len([x for x in case_vals if x[0] == '>' or x[0] == '<']):
        new_meta_pd = meta_pd.copy()
        for case_val in case_vals:
            if case_val[0] == '>':
                new_meta_pd = new_meta_pd[new_meta_pd[case_var].astype(float) >= float(case_val[1:])].copy()
            elif case_val[0] == '<':
                new_meta_pd = new_meta_pd[new_meta_pd[case_var].astype(float) <= float(case_val[1:])].copy()
    else:
        new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()
    return new_meta_pd


def get_new_alpha_div(case: str, div_qza: str, cur_rad: str,
                      new_meta_pd: pd.DataFrame, cur_sh: TextIO) -> str:
    """
    Subset the current alpha diversity vector and export.

    :param case: concatenation of the current subset / test groups.
    :param div_qza: alpha diversity qiime2 Artefact.
    :param cur_rad: radical of the file.
    :param new_meta_pd: metadata table.
    :param cur_sh: writing file handle.
    :return: subsetted alpha diversity qiime2 Artefact.
    """
    new_div = '%s.qza' % cur_rad
    if 'ALL' in case:
        cur_sh.write('echo "cp %s %s"\n' % (div_qza, new_div))
        cur_sh.write('cp %s %s\n' % (div_qza, new_div))
    else:
        new_tsv = '%s.tsv' % cur_rad
        new_tsv_pd = pd.read_csv(div_qza.replace('.qza', '.tsv'), header=0, sep='\t', dtype=str)
        new_tsv_pd.rename(columns={new_tsv_pd.columns.tolist()[0]: 'Feature ID'}, inplace=True)
        new_tsv_pd.set_index('Feature ID', inplace=True)
        new_tsv_pd = new_tsv_pd.loc[new_meta_pd.index.tolist(), :]
        new_tsv_pd.reset_index().to_csv(new_tsv, index=False, sep='\t')
        cmd = run_import(new_tsv, new_div, 'SampleData[AlphaDiversity]')
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)
    return new_div


def get_padded_new_rows_list(new_rows, max_new_rows):
    # create a new 'Not available' matrix of shape (n_features x n_fields)
    padded_new_rows_list = []
    for padded_row in new_rows:
        # update row if not of max length
        n_pad = max_new_rows - len(padded_row)
        if n_pad:
            to_pad = ['Not_available']*n_pad
            padded_row = padded_row + to_pad
        padded_new_rows_list.append(padded_row)
    return padded_new_rows_list


def add_alpha_level_label(taxa, padded_new_rows_list, max_new_rows):
    # add a alpha label for rank-identification prupose
    ALPHA = 'ABCDEFGHIJKL'
    cols = ['Node_level_%s' % (ALPHA[idx]) for idx in range(1, max_new_rows + 1)]
    padded_new_rows_pd = pd.DataFrame(
        padded_new_rows_list,
        index = taxa,
        columns = cols
    )
    return padded_new_rows_pd


def extend_split_taxonomy(split_taxa_pd: pd.DataFrame):
    to_concat = []
    for col in split_taxa_pd.columns.tolist():
        if split_taxa_pd[col].unique().size > 50:
            continue
        split_taxa_dummy = split_taxa_pd[col].str.get_dummies()
        split_taxa_dummy.columns = ['%s__%s' % (x, col) for x in split_taxa_dummy.columns]
        to_concat.append(split_taxa_dummy)
    if len(to_concat):
        return pd.concat(to_concat, axis=1)
    else:
        return pd.DataFrame()


def get_split_taxonomy(taxa, taxo_sep=';'):

    # get the taxon name split per "taxon" level
    split_chars = taxo_sep
    if len([1 for x in taxa if '|' in x]) > (0.5 * len(taxa)):
        split_chars += "|\|"
    if len([1 for x in taxa if '.' in x]) > (0.5 * len(taxa)):
        split_chars += "|\."

    split_taxa_pd = pd.DataFrame([
        pd.Series(
            [x.strip() for x in re.split(split_chars, str(taxon))
             if len(x.strip()) and x[0] != 'x']
        ) for taxon in taxa
    ])
    ALPHA = 'ABCDEFGHIJKL'
    split_taxa_pd_cols =['Taxolevel_%s' % (ALPHA[idx]) for idx in range(split_taxa_pd.shape[1])]
    split_taxa_pd.columns = split_taxa_pd_cols
    # get the max number of fields
    # max_new_rows = max([len(new_row) for new_row in new_rows])
    # padded_new_rows_list = get_padded_new_rows_list(new_rows, max_new_rows)
    # padded_new_rows_pd = add_alpha_level_label(taxa, padded_new_rows_list, max_new_rows)
    padded_extended_pd = extend_split_taxonomy(split_taxa_pd)
    if padded_extended_pd.shape[0]:
        split_taxa_pd = pd.concat([split_taxa_pd, padded_extended_pd], axis=1)
    # split_taxa_pd = split_taxa_pd.reset_index().rename(columns={'index': 'Taxon'})
    return split_taxa_pd

