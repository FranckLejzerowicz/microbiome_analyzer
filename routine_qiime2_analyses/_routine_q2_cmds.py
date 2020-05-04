# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from typing import TextIO
from os.path import isdir, isfile, splitext

from skbio.stats.ordination import OrdinationResults


def get_subset(tsv_pd: pd.DataFrame, subset: str, meta_subset: str, subset_regex: list) -> int:
    """
    Make a feature metadata from the regex
    to get the names of the features to keep.

    :param tsv_pd: Table containing the features to subset.
    :param meta_subset: Feature metadata to create.
    :param subset_regex: subsetting regex.
    """
    to_keep_feats = {}
    for regex in subset_regex:
        to_keep_feats[regex.lower()] = tsv_pd.index.str.lower().str.contains(regex.lower())
    to_keep_feats_pd = pd.DataFrame(to_keep_feats)
    to_keep_feats = to_keep_feats_pd.any(axis=1)

    feats_subset = tsv_pd.index[to_keep_feats].tolist()
    subset_pd = pd.DataFrame({'Feature ID': feats_subset, 'Subset': [subset]*len(feats_subset)})
    subset_pd.to_csv(meta_subset, index=False, sep='\t')

    return len(feats_subset)


def write_filter_features(qza: str, qza_subset: str, meta_subset: str, cur_sh: TextIO) -> None:
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
    cmd = 'qiime feature-table filter-features \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % meta_subset
    cmd += '--o-filtered-table %s\n' % qza_subset
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
                    conditionals_tsv: str, biplot_tsv: str,
                    batch: str, learn: str, epoch: str,
                    prior: str, thresh_feat: str, latent_dim: str,
                    train_column: str, n_example: str, gpu: bool,
                    standalone: bool, cur_sh: TextIO) -> None:
    """
    Performs bi-loglinear multinomial regression and calculates the
    conditional probability ranks of metabolite co-occurence given the microbe
    presence.

    :param meta_fp:
    :param qza1:
    :param qza2:
    :param res_dir:
    :param conditionals_tsv:
    :param biplot_tsv:
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
        cmd += '--checkpoint-interval 10 \\\n'
        cmd += '--summary-interval 10 \\\n'
        cmd += '--summary-dir %s \\\n' % res_dir
        cmd += '--ranks-file %s\n' % conditionals_tsv
    else:
        conditionals_qza = '%s.qza' % splitext(conditionals_tsv)[0]
        biplot_qza = '%s.qza' % splitext(biplot_tsv)[0]
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
        cmd += '--p-summary-interval 10 \\\n'
        cmd += '--o-conditionals %s \\\n' % conditionals_qza
        cmd += '--o-conditional-biplot %s\n' % biplot_qza
        if not isfile(conditionals_tsv):
            cmd += run_export(conditionals_qza, conditionals_tsv, '')

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
                       thresh_sample: str, thresh_feat: str, n_random: str,
                       diffs: str, diffs_qza: str, stats: str, plot: str,
                       base_diff_qza: str, base_stats: str, base_plot: str,
                       tensor: str, tensor_dir: str, cur_sh: TextIO) -> None:
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
    :param n_random:
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

    # could evolve
    base_formula = '1'

    if not isfile(new_qza):
        cmd = filter_feature_table(qza, new_qza, new_meta)
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if not isfile(diffs_qza):
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
        cmd += ' --p-num-random-test-examples %s \\\n' % n_random
        cmd += ' --p-summary-interval 2 \\\n'
        cmd += ' --o-differentials %s \\\n' % diffs_qza
        cmd += ' --o-regression-stats %s \\\n' % stats
        cmd += ' --o-regression-biplot %s\n' % plot
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if not isfile(diffs):
        cmd = run_export(diffs_qza, diffs, '')
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if not isfile(base_diff_qza):
        cmd = '\nqiime songbird multinomial \\\n'
        cmd += ' --i-table %s \\\n' % new_qza
        cmd += ' --m-metadata-file %s \\\n' % new_meta
        cmd += ' --p-formula "%s" \\\n' % base_formula
        cmd += ' --p-epochs %s \\\n' % epoch
        cmd += ' --p-batch-size %s \\\n' % batch
        cmd += ' --p-differential-prior %s \\\n' % diff_prior
        cmd += ' --p-learning-rate %s \\\n' % learn
        cmd += ' --p-min-sample-count %s \\\n' % thresh_sample
        cmd += ' --p-min-feature-count %s \\\n' % thresh_feat
        cmd += ' --p-num-random-test-examples %s \\\n' % n_random
        cmd += ' --p-summary-interval 2 \\\n'
        cmd += ' --o-differentials %s \\\n' % base_diff_qza
        cmd += ' --o-regression-stats %s \\\n' % base_stats
        cmd += ' --o-regression-biplot %s\n' % base_plot
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if not isfile(tensor):
        cmd = '\n\nqiime songbird summarize-paired \\\n'
        cmd += ' --i-regression-stats %s \\\n' % stats
        cmd += ' --i-baseline-stats %s \\\n' % base_stats
        cmd += ' --o-visualization %s\n' % tensor
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)

    if not isdir(tensor_dir):
        cmd = run_export(tensor, tensor_dir, '')
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n' % cmd)


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
        elif 'songbird' in typ:
            cmd += 'mv %s/index.html %s\n' % (splitext(output_path)[0], output_path)
        else:
            cmd += 'mv %s/*.tsv %s\n' % (splitext(output_path)[0], output_path)
        cmd += 'rm -rf %s\n' % splitext(output_path)[0]
    return cmd


def write_diversity_beta(out_fp: str, datasets_phylo: dict, trees: dict, dat: str,
                         qza: str, metric: str, cur_sh: TextIO) -> bool:
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
        if not datasets_phylo[dat][0] or dat not in trees:
            return True
        cmd = 'qiime diversity beta-phylogenetic \\\n'
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
    return False


def write_diversity_pcoa(DM: str, out_pcoa: str, cur_sh: TextIO) -> None:
    """
    Apply principal coordinate analysis.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/

    :param DM: The distance matrix on which PCoA should be computed.
    :param out_pcoa: The resulting PCoA matrix.
    :param cur_sh: writing file handle.
    """
    cmd = 'qiime diversity pcoa \\\n'
    cmd += '--i-distance-matrix %s \\\n' % DM
    cmd += '--o-pcoa %s\n' % out_pcoa
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_diversity_biplot(in_tab: str, out_pcoa: str, out_biplot: str, cur_sh: TextIO) -> None:
    """
    pcoa-biplot: Principal Coordinate Analysis Biplot.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa-biplot/

    :param DM: The distance matrix on which PCoA should be computed.
    :param out_pcoa: The resulting PCoA matrix.
    :param cur_sh: writing file handle.
    """
    rel = '%s_rel.qza' % splitext(in_tab)[0]
    cmd = 'qiime feature-table relative-frequency \\\n'
    cmd += '--i-table %s \\\n' % in_tab
    cmd += '--o-relative-frequency-table %s\n' % rel
    cmd += 'qiime diversity pcoa-biplot \\\n'
    cmd += '--i-pcoa %s \\\n' % out_pcoa
    cmd += '--i-features %s \\\n' % rel
    cmd += '--o-biplot %s\n' % out_biplot
    cmd += 'rm %s\n' % rel
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)

    out_biplot_txt = '%s.txt' % splitext(out_biplot)[0]
    cmd = run_export(out_biplot, out_biplot_txt, 'biplot')
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write('%s\n\n' % cmd)


def write_emperor(meta: str, pcoa_biplot: str, out_plot: str, cur_sh: TextIO, taxonomy: str) -> None:
    """
    Generates an interactive ordination plot where the user can visually
    integrate sample metadata.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    :param meta: The sample metadata.
    :param pcoa: The principal coordinates matrix to be plotted.
    :param out_plot: VISUALIZATION.
    :param cur_sh: writing file handle.
    """
    if taxonomy:
        pcoa_biplot_txt = '%s.txt' % splitext(pcoa_biplot)[0]
        if isfile(pcoa_biplot_txt):
            ordi = OrdinationResults.read(pcoa_biplot_txt)
            ordi.features = ordi.features.iloc[:,:3]
            ordi.samples= ordi.samples.iloc[:,:3]
            ordi.eigvals = ordi.eigvals[:3]
            ordi.proportion_explained = ordi.proportion_explained[:3]
            ordi.write(pcoa_biplot_txt)
        cmd = run_import(pcoa_biplot_txt, pcoa_biplot, "PCoAResults % Properties('biplot')")
        cur_sh.write('echo "%s"\n' % cmd)
        cur_sh.write('%s\n\n' % cmd)

        cmd = 'qiime emperor biplot \\\n'
        cmd += '--i-biplot %s \\\n' % pcoa_biplot
        cmd += '--m-sample-metadata-file %s \\\n' % meta
        cmd += '--m-feature-metadata-file %s \\\n' % taxonomy
        cmd += '--p-number-of-features 20 \\\n'
        cmd += '--o-visualization %s\n' % out_plot
    else:
        cmd = 'qiime emperor plot \\\n'
        cmd += '--i-pcoa %s \\\n' % pcoa_biplot
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
    cmd = '\nqiime feature-table filter-samples \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--o-filtered-table %s\n' % new_qza
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
                                            qza: str, new_qza: str, testing_group: str,
                                            new_qzv: str, cur_sh: TextIO) -> None:
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
    cmd = 'qiime diversity filter-distance-matrix \\\n'
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--i-distance-matrix %s \\\n' % mat_qza
    cmd += '--o-filtered-distance-matrix %s\n' % new_mat_qza
    cmd += 'qiime feature-table filter-samples \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--o-filtered-table %s\n' % new_qza
    cmd += 'qiime diversity beta-group-significance \\\n'
    cmd += '--i-distance-matrix %s \\\n' % new_mat_qza
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--m-metadata-column "%s" \\\n' % testing_group
    cmd += '--p-permutations 2999 \\\n'
    cmd += '--o-visualization %s\n' % new_qzv
    cur_sh.write('echo "%s"\n' % cmd)
    cur_sh.write(cmd)


def write_diversity_adonis(new_meta: str, mat_qza: str, new_mat_qza: str,
                           qza: str, new_qza: str, formula: str,
                           new_qzv: str, cur_sh: TextIO) -> None:
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
    cmd += 'qiime feature-table filter-samples \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--o-filtered-table %s\n' % new_qza
    cmd += 'qiime diversity adonis \\\n'
    cmd += '--i-distance-matrix %s \\\n' % new_mat_qza
    cmd += '--m-metadata-file %s \\\n' % new_meta
    cmd += '--p-formula "%s" \\\n' % formula
    cmd += '--p-permutations 2999 \\\n'
    cmd += '--p-n-jobs 6 \\\n'
    cmd += '--o-visualization %s\n' % new_qzv
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


def get_case(case_vals: list, metric: str, case_var: str, form: str = None) -> str:
    """
    Get the current case, which is the concatenation of:
     - diversity metric.
     - metadata variable.
     - metadata variable's values.
     - formula.
    :param case_vals: variable's values
    :param metric: diversity metric.
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
    case = case.replace('__', '_').replace(' ', '-')
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
