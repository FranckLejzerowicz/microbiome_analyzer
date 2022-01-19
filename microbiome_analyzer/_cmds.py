# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
from os.path import basename, dirname, isfile, splitext
from skbio.stats.ordination import OrdinationResults


def run_import(
        input_path: str,
        output_path: str,
        typ: str) -> str:
    """
    Return the import qiime2 command.
    
    Parameters
    ----------
    input_path : str
        Input file to import to qiime2 artefact.
    output_path : str
        Output qiime2 artefact.
    typ : str
        qiime2 type.

    Returns
    -------
    cmd : str
        Qiime2 command.
    """
    cmd = ''
    infile = input_path
    if typ.startswith("FeatureTable"):
        if not input_path.endswith('biom'):
            infile = '%s.biom' % splitext(input_path)[0]
            cmd += 'biom convert \\\n'
            cmd += '  -i %s \\\n' % input_path
            cmd += '  -o %s \\\n' % infile
            cmd += '  --table-type="OTU table" \\\n'
            cmd += '  --to-hdf5\n\n'
        cmd += 'qiime tools import \\\n'
        cmd += '  --input-path %s \\\n' % infile
        cmd += '  --output-path %s \\\n' % output_path
        cmd += '  --type "FeatureTable[Frequency]"\n'
    else:
        cmd += 'qiime tools import \\\n'
        cmd += '  --input-path %s \\\n' % infile
        cmd += '  --output-path %s \\\n' % output_path
        cmd += '  --type "%s"\n' % typ
    return cmd


def run_export(
        input_path: str,
        output_path: str,
        typ: str
) -> str:
    """
    Return the export qiime2 command.

    Parameters
    ----------
    input_path : str
        Input qiime2 artefact.
    output_path : str
        Output file to import to qiime2 artefact.
    typ : str
        qiime2 type.

    Returns
    -------
    cmd : str
        Qiime2 command.
    """
    inp = splitext(input_path)[0]
    out = splitext(output_path)[0]

    cmd = 'qiime tools export \\\n'
    cmd += '  --input-path %s \\\n' % input_path
    cmd += '  --output-path %s\n' % out
    if typ.startswith("FeatureTable"):
        if output_path.endswith('biom'):
            cmd += 'mv %s/*.biom %s\n' % (out, output_path)
            cmd += 'rm -rf %s\n' % inp
        else:
            cur_biom = '%s.biom' % out
            cmd += 'mv %s/*.biom %s\n' % (out, cur_biom)
            cmd += 'biom convert'
            cmd += '  -i %s \\\n' % cur_biom
            cmd += '  -o %s.tmp \\\n' % output_path
            cmd += '  --to-tsv\n\n'
            cmd += 'tail -n +2 %s.tmp > %s\n\n' % (output_path, output_path)
            cmd += 'rm -rf %s %s.tmp\n' % (out, output_path)
    else:
        if 'phylogeny' in typ:
            cmd += 'mv %s/*.nwk %s\n' % (out, output_path)
        elif typ in ['pcoa', 'biplot', 'mmvec']:
            cmd += 'mv %s/*.txt %s\n' % (out, output_path)
        elif typ in ['perms', 'mantel', 'songbird', 'mmvec_summary']:
            cmd += 'mv %s/index.html %s\n' % (out, output_path)
        else:
            cmd += 'mv %s/*.tsv %s\n' % (out, output_path)
        cmd += 'rm -rf %s\n' % out
    return cmd


def write_rarefy(
        qza: str,
        qza_out: str,
        depth: int
) -> str:
    """
    Subsample frequencies from all samples so that the sum of frequencies in
    each sample is equal to sampling-depth.
    https://docs.qiime2.org/2019.10/plugins/available/feature-table/rarefy/

    Parameters
    ----------
    qza
    qza_out
    depth

    Returns
    -------

    """
    cmd = 'qiime feature-table rarefy \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--p-sampling-depth %s \\\n' % str(depth)
    cmd += '--o-rarefied-table %s\n' % qza_out
    return cmd


def write_fasta(
        out_fp_seqs_fasta: str,
        out_fp_seqs_qza: str,
        biom: biom.Table,
        tsv_fp: str = ''
) -> str:
    """
    Write the fasta sequences.

    Parameters
    ----------
    out_fp_seqs_fasta
    out_fp_seqs_qza
    biom
    tsv_fp

    Returns
    -------

    """
    with open(out_fp_seqs_fasta, 'w') as fas_o:
        for seq in biom.ids(axis='observation'):
            fas_o.write('>%s\n%s\n' % (seq.strip(), seq.strip()))
    cmd = '# Write features as fasta file:\n'
    cmd += '#  - Features from: %s\n' % tsv_fp
    cmd += '# Snippet:\n'
    cmd += '# ```:\n'
    cmd += "# with open(fasta_out, 'w') as o:\n"
    cmd += "#     for seq in tsv_pd.index:\n"
    cmd += "#         o.write('>%s\\n%s\\n' % (seq.strip(), seq.strip()))\n"
    cmd += '# ```:\n'
    cmd += run_import(
        out_fp_seqs_fasta, out_fp_seqs_qza, 'FeatureData[Sequence]')
    return cmd


def write_taxonomy_sklearn(
        out_qza: str,
        out_fp_seqs_qza: str,
        ref_classifier_qza: str
) -> str:
    """
    Classify reads by taxon using a fitted classifier.
    https://docs.qiime2.org/2020.2/plugins/available/feature-classifier/classify-sklearn

    Parameters
    ----------
    out_qza
    out_fp_seqs_qza
    ref_classifier_qza

    Returns
    -------

    """
    cmd = 'qiime feature-classifier classify-sklearn \\\n'
    cmd += '--i-reads %s \\\n' % out_fp_seqs_qza
    cmd += '--i-classifier %s \\\n' % ref_classifier_qza
    cmd += '--p-n-jobs %s \\\n' % '4'
    cmd += '--o-classification %s\n' % out_qza
    return cmd


def write_collapse_taxo(
        tab_qza: str,
        tax_qza: str,
        collapsed_qza: str,
        collapsed_tsv: str,
        meta_fp: str,
        collapsed_meta: str,
        level: int,
        remove_empty: set
) -> str:
    """

    Parameters
    ----------
    tab_qza
    tax_qza
    collapsed_qza
    collapsed_tsv
    meta_fp
    collapsed_meta
    level
    remove_empty

    Returns
    -------

    """
    cmd = ''
    if not isfile(collapsed_qza):
        cmd += 'qiime taxa collapse \\\n'
        cmd += '--i-table %s \\\n' % tab_qza
        cmd += '--i-taxonomy %s \\\n' % tax_qza
        cmd += '--p-level %s \\\n' % level
        cmd += '--o-collapsed-table %s\n\n' % collapsed_qza
    if remove_empty:
        tax_tmp = '%s_filtempty%s.tsv' % (splitext(tax_qza)[0], level)
        cmd += '\n# Metadata to remove empty/catch-all features: %s\n' % (
            tax_tmp)
        with open(tax_tmp, 'w') as o:
            o.write('Feature ID\tTaxon\n')
            for tax in remove_empty:
                cmd += '#  - %s\n' % tax
                o.write('%s\tremove\n' % tax)
        cmd += '\n\nqiime feature-table filter-features \\\n'
        cmd += '--i-table %s \\\n' % collapsed_qza
        cmd += '--m-metadata-file %s \\\n' % tax_tmp
        cmd += '--o-filtered-table %s2.qza \\\n' % collapsed_qza
        cmd += '--p-exclude-ids\n'
        cmd += 'mv %s2.qza %s\n' % (collapsed_qza, collapsed_qza)
    if not isfile(collapsed_tsv):
        cmd += run_export(collapsed_qza, collapsed_tsv, 'FeatureTable')
    if not isfile(collapsed_meta):
        cmd += 'cp %s %s\n' % (meta_fp, collapsed_meta)
    return cmd


def write_collapse(
        tab_qza: str,
        tax_qza: str,
        collapsed_qza: str,
        collapsed_tsv: str,
        meta_fp: str,
        collapsed_meta: str,
        level: int,
        remove_empty: set
) -> str:
    """

    Parameters
    ----------
    tab_qza
    tax_qza
    collapsed_qza
    collapsed_tsv
    meta_fp
    collapsed_meta
    level
    remove_empty

    Returns
    -------

    """
    cmd = ''
    if not isfile(collapsed_qza):
        cmd += 'qiime taxa collapse \\\n'
        cmd += '--i-table %s \\\n' % tab_qza
        cmd += '--i-taxonomy %s \\\n' % tax_qza
        cmd += '--p-level %s \\\n' % level
        cmd += '--o-collapsed-table %s\n\n' % collapsed_qza
    if remove_empty:
        tax_tmp = '%s_filtempty%s.tsv' % (splitext(tax_qza)[0], level)
        cmd += '\n# Metadata to remove empty/catch-all features: %s\n' % (
            tax_tmp)
        with open(tax_tmp, 'w') as o:
            o.write('Feature ID\tTaxon\n')
            for tax in remove_empty:
                cmd += '#  - %s\n' % tax
                o.write('%s\tremove\n' % tax)
        cmd += '\n\nqiime feature-table filter-features \\\n'
        cmd += '--i-table %s \\\n' % collapsed_qza
        cmd += '--m-metadata-file %s \\\n' % tax_tmp
        cmd += '--o-filtered-table %s2.qza \\\n' % collapsed_qza
        cmd += '--p-exclude-ids\n'
        cmd += 'mv %s2.qza %s\n' % (collapsed_qza, collapsed_qza)
    if not isfile(collapsed_tsv):
        cmd += run_export(collapsed_qza, collapsed_tsv, 'FeatureTable')
    if not isfile(collapsed_meta):
        cmd += 'cp %s %s\n' % (meta_fp, collapsed_meta)
    return cmd


def write_sepp(
        out_fp_seqs_qza: str,
        ref_tree_qza: str,
        out_fp_sepp_tree: str,
        qza: str,
        qza_in: str,
        qza_in_tsv: str,
        qza_out: str
) -> str:
    """
    Perform fragment insertion of sequences using the SEPP algorithm.
    https://docs.qiime2.org/2019.10/plugins/available/fragment-insertion/sepp/

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

    Parameters
    ----------
    out_fp_seqs_qza
    ref_tree_qza
    out_fp_sepp_tree
    qza
    qza_in
    qza_in_tsv
    qza_out

    Returns
    -------

    """
    out_fp_sepp_plac = '%s/plac_%s' % (dirname(out_fp_sepp_tree),
                                       basename(out_fp_sepp_tree)[4:])
    cmd = ''
    if not isfile(out_fp_sepp_tree):
        cmd += 'qiime fragment-insertion sepp \\\n'
        cmd += '--i-representative-sequences %s \\\n' % out_fp_seqs_qza
        cmd += '--i-reference-database %s \\\n' % ref_tree_qza
        cmd += '--o-tree %s \\\n' % out_fp_sepp_tree
        cmd += '--o-placements %s \\\n' % out_fp_sepp_plac
        cmd += '--p-threads 24\n'
    if not isfile(qza_in_tsv):
        cmd += 'qiime fragment-insertion filter-features \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--i-tree %s \\\n' % out_fp_sepp_tree
        cmd += '--o-filtered-table %s \\\n' % qza_in
        cmd += '--o-removed-table %s\n' % qza_out
        cmd += run_export(qza_in, qza_in_tsv, 'FeatureTable')
    return cmd


def write_alpha(
        qza: str,
        qza_out: str,
        datasets_phylo: list,
        trees: list,
        metric: str
) -> str:
    """
    Computes a user-specified alpha diversity metric for all samples
    in a feature table.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/alpha/

    Parameters
    ----------
    qza
    qza_out
    datasets_phylo
    trees
    metric

    Returns
    -------

    """
    if metric in ['faith_pd']:
        if not datasets_phylo[0] or not trees or not isfile(trees[1]):
            return ''
        cmd = 'qiime diversity alpha-phylogenetic \\\n'
        if datasets_phylo[1]:
            cmd += '--i-table %s \\\n' % trees[0]
        else:
            cmd += '--i-table %s \\\n' % qza
        cmd += '--i-phylogeny %s \\\n' % trees[1]
        cmd += '--p-metric %s \\\n' % metric
        cmd += '--o-alpha-diversity %s\n' % qza_out
    else:
        cmd = 'qiime diversity alpha \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--p-metric %s \\\n' % metric
        cmd += '--o-alpha-diversity %s\n\n' % qza_out
    return cmd


def write_filter(
        qza: str,
        qza_subset: str,
        meta_subset: str
) -> str:
    """
    filter-features: Filter features from table¶
    https://docs.qiime2.org/2020.2/plugins/available/feature-table/filter-features/

    Filter features from table based on frequency and/or metadata. Any samples
    with a frequency of zero after feature filtering will also be removed. See
    the filtering tutorial on https://docs.qiime2.org for additional details.

    Parameters
    ----------
    qza
    qza_subset
    meta_subset

    Returns
    -------

    """

    cmd = 'qiime feature-table filter-features \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % meta_subset
    cmd += '--o-filtered-table %s\n\n' % qza_subset
    return cmd


def write_barplots(
        qza: str,
        qzv: str,
        meta: str,
        tax_qza: str
) -> str:
    """
    Visualize taxonomy with an interactive bar plot¶
    https://docs.qiime2.org/2020.2/plugins/available/taxa/barplot/

    This visualizer produces an interactive barplot visualization of
    taxonomies. Interactive features include multi-level sorting, plot
    recoloring, sample relabeling, and SVG figure export.

    Parameters
    ----------
    qza
    qzv
    meta
    tax_qza

    Returns
    -------

    """
    cmd = 'qiime taxa barplot \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--i-taxonomy %s \\\n' % tax_qza
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s \\\n' % qzv
    return cmd


def write_tabulate(
        out_qza: str,
        alphas: list
) -> str:
    """
    Generate a tabular view of Metadata. The output visualization supports
    interactive filtering, sorting, and exporting to common file formats.
    https://docs.qiime2.org/2019.10/plugins/available/metadata/tabulate/

    Parameters
    ----------
    out_qza
    alphas

    Returns
    -------

    """
    cmd = 'qiime metadata tabulate \\\n'
    cmd += '--o-visualization %s \\\n' % out_qza
    for alpha in alphas:
        cmd += '--m-input-file %s \\\n' % alpha[0]
    cmd += '\n\n'
    return cmd


def write_alpha_correlation(
        qza: str,
        qzv: str,
        method: str,
        meta: str
) -> str:
    """
    Determine whether numeric sample metadata columns
    are correlated with alpha diversity.
    https://docs.qiime2.org/2019.10/plugins/available/
    diversity/alpha-correlation/

    Parameters
    ----------
    qza
    qzv
    method
    meta

    Returns
    -------

    """
    cmd = 'qiime diversity alpha-correlation \\\n'
    cmd += '--i-alpha-diversity %s \\\n' % qza
    cmd += '--p-method %s \\\n' % method
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s\n\n' % qzv
    return cmd


def write_alpha_rarefaction(
        qza: str,
        qzv: str,
        metric: str,
        phylo: tuple,
        tree: tuple,
        meta: str,
        raref: str
) -> str:
    """

    Parameters
    ----------
    qza
    qzv
    metric
    phylo
    tree
    meta
    raref

    Returns
    -------

    """
    cmd = 'qiime diversity alpha-rarefaction \\\n'
    if metric in ['faith_pd']:
        if not phylo[0] or not tree[-1]:
            return ''
        if phylo[1]:
            cmd += '--i-table %s \\\n' % tree[0]
        else:
            cmd += '--i-table %s \\\n' % qza
        cmd += '--i-phylogeny %s \\\n' % tree[1]
    else:
        cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--p-metrics %s \\\n' % metric
    if raref:
        cmd += '--p-max-depth %s \\\n' % raref.split('raref')[-1]
    else:
        cmd += '--p-max-depth 5000 \\\n'
    cmd += '--o-visualization %s\n\n' % qzv
    return cmd


def write_volatility(
        meta: str,
        qzv: str,
        timepoint: str,
        host: str
) -> str:
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

    Parameters
    ----------
    meta
    qzv
    timepoint
    host

    Returns
    -------

    """
    cmd = 'qiime longitudinal volatility \\\n'
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--p-state-column "%s" \\\n' % timepoint
    cmd += '--p-individual-id-column "%s" \\\n' % host
    cmd += '--o-visualization %s\n' % qzv
    return cmd


def write_beta(
        qza: str,
        dm_qza: str,
        phylo: tuple,
        metric: str,
        tree: tuple,
        config
) -> str:
    """
    Computes a user-specified beta diversity metric for all pairs of samples
    in a feature table.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta/
    and
    Computes a user-specified phylogenetic beta diversity metric for all pairs
    of samples in a feature table.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-phylogenetic/

    Parameters
    ----------
    qza
    dm_qza
    phylo
    metric
    tree
    config

    Returns
    -------

    """
    cmd = ''
    if 'unifrac' in metric:
        if phylo[0] and tree[1] and isfile(tree[1]):
            cmd += 'qiime diversity beta-phylogenetic \\\n'
            if tree[0]:
                cmd += '--i-table %s \\\n' % tree[0]
            else:
                cmd += '--i-table %s \\\n' % qza
            cmd += '--i-phylogeny %s \\\n' % tree[1]
    else:
        cmd = 'qiime diversity beta \\\n'
        cmd += '--i-table %s \\\n' % qza
    if cmd:
        cmd += '--p-metric %s \\\n' % metric
        n_nodes = config.run_params['beta']['n_nodes']
        n_procs = config.run_params['beta']['n_procs']
        if float(config.qiime_env.split('-')[1]) >= 2020.8 and 'phylo' in cmd:
            cmd += '--p-threads %s \\\n' % (int(n_nodes) * int(n_procs))
        else:
            cmd += '--p-n-jobs %s \\\n' % (int(n_nodes) * int(n_procs))
        cmd += '--o-distance-matrix %s\n' % dm_qza
    return cmd


def write_deicode(
        qza: str,
        meta: str,
        new_qza: str,
        ordi: str,
        dm_qza: str,
        qzv: str
) -> str:
    """
    Performs robust center log-ratio transform robust PCA and
    ranks the features by the loadings of the resulting SVD.
    https://library.qiime2.org/plugins/deicode/19/

    Parameters
    ----------
    qza
    meta
    new_qza
    ordi
    dm_qza
    qzv

    Returns
    -------

    """
    cmd = ''
    if not isfile(new_qza):
        cmd += 'qiime feature-table filter-samples \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--m-metadata-file %s \\\n' % meta
        cmd += '--o-filtered-table %s\n\n' % new_qza
    if not isfile(ordi) or not isfile(dm_qza):
        cmd += 'qiime deicode rpca \\\n'
        cmd += '--i-table %s \\\n' % new_qza
        # cmd += '--p-min-feature-count 10 \\\n'
        # cmd += '--p-min-sample-count 500 \\\n'
        cmd += '--p-n-components 2 \\\n'
        cmd += '--o-biplot %s \\\n' % ordi
        cmd += '--o-distance-matrix %s\n\n' % dm_qza
    cmd += 'qiime emperor biplot \\\n'
    cmd += '--i-biplot %s \\\n' % ordi
    cmd += '--m-sample-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s \\\n' % qzv
    cmd += '--p-number-of-features 10\n\n'
    cmd += 'rm %s %s\n\n' % (meta, new_qza)
    return cmd


def write_pcoa(
        dm: str,
        dm_filt: str,
        meta: str,
        meta_met: str,
        group: str,
        pcoa: str
) -> str:
    """
    Apply principal coordinate analysis.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/

    Parameters
    ----------
    dm
    dm_filt
    meta
    meta_met
    group
    pcoa

    Returns
    -------

    """
    cmd = ''
    if group != 'ALL':
        cmd += 'qiime diversity filter-distance-matrix \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm
        cmd += '--m-metadata-file %s \\\n' % meta_met
        cmd += '--o-filtered-distance-matrix %s\n' % dm_filt
    else:
        dm_filt = dm
    cmd += 'qiime diversity pcoa \\\n'
    cmd += '--i-distance-matrix %s \\\n' % dm_filt
    cmd += '--o-pcoa %s\n\n' % pcoa
    if group != 'ALL':
        cmd += 'rm %s\n\n' % dm_filt
    cmd += 'mv %s %s.tsv\n\n' % (meta_met, meta)
    return cmd


def write_biplot(
        tsv: str,
        qza: str,
        meta: str,
        pcoa: str,
        tax: str,
        biplot: str,
        biplot_tax: str
) -> str:
    """
    pcoa-biplot: Principal Coordinate Analysis Biplot.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa-biplot/

    Parameters
    ----------
    tsv
    qza
    meta
    pcoa
    tax
    biplot
    biplot_tax

    Returns
    -------

    """
    cmd = ''
    if isfile(tax):
        tax_dict = {}
        with open(tax) as f, open(biplot_tax, 'w') as o_tax:
            o_tax.write('Feature ID\tTaxon\tPrevious ID\n')
            n = 0
            for ldx, line in enumerate(f):
                line_split = line.strip().split('\t')
                if len(line_split) == 1:
                    line_split = [line_split[0], line_split[0]]
                if ldx and not line.startswith('#q2:types'):
                    new = 'x__%s;%s' % (n, line_split[1])
                    tax_dict[line_split[0]] = new
                    o_tax.write('%s\t%s\t%s\n' % (new, new, line_split[0]))
                    n += 1
        tab_tsv = '%s_table.tsv' % splitext(biplot)[0]
        with open(tsv) as f, open(tab_tsv, 'w') as o_tab:
            for ldx, line in enumerate(f):
                t = line.strip().split('\t')
                if t[0] in tax_dict:
                    o_tab.write('%s\t%s\n' % (tax_dict[t[0]], '\t'.join(t[1:])))
                else:
                    o_tab.write(line)
        tab_qza = '%s.qza' % splitext(tab_tsv)[0]
        cmd += run_import(tab_tsv, tab_qza, 'FeatureTable[Frequency]')
    else:
        tab_qza = qza

    tab_rel_qza_tmp = '%s_rel_tmp.qza' % splitext(tab_qza)[0]
    tab_rel_qza = '%s_rel.qza' % splitext(tab_qza)[0]

    cmd += '\nqiime feature-table relative-frequency \\\n'
    cmd += '--i-table %s \\\n' % tab_qza
    cmd += '--o-relative-frequency-table %s\n\n' % tab_rel_qza_tmp
    cmd += '\nqiime feature-table filter-samples \\\n'
    cmd += '--i-table %s \\\n' % tab_rel_qza_tmp
    cmd += '--m-metadata-file %s.tsv \\\n' % meta
    cmd += '--o-filtered-table %s\n\n' % tab_rel_qza

    cmd += 'qiime diversity pcoa-biplot \\\n'
    cmd += '--i-pcoa %s \\\n' % pcoa
    cmd += '--i-features %s \\\n' % tab_rel_qza
    cmd += '--o-biplot %s\n\n' % biplot

    cmd += 'rm %s %s\n' % (tab_rel_qza_tmp, tab_rel_qza)

    out_biplot_txt = '%s.txt' % splitext(biplot)[0]
    cmd += run_export(biplot, out_biplot_txt, 'biplot')
    return cmd


def write_emperor(
        pcoa: str,
        qzv: str,
        meta: str
) -> str:
    """
    Generates an interactive ordination plot where the user can visually
    integrate sample metadata.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    Parameters
    ----------
    pcoa
    qzv
    meta

    Returns
    -------

    """
    cmd = 'qiime emperor plot \\\n'
    cmd += '--i-pcoa %s \\\n' % pcoa
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-visualization %s\n\n' % qzv
    return cmd


def write_emperor_biplot(
        biplot: str,
        biplot_tax: str,
        meta: str,
        qzv: str,
        tax_pd: pd.DataFrame
) -> str:
    """
    Generates an interactive ordination plot where the user can visually
    integrate sample metadata.
    https://docs.qiime2.org/2019.10/plugins/available/emperor/

    Parameters
    ----------
    biplot
    biplot_tax
    meta
    qzv
    tax_pd

    Returns
    -------

    """
    cmd = ''
    biplot_txt = '%s.txt' % splitext(biplot)[0]
    if isfile(biplot_txt):
        ordi = OrdinationResults.read(biplot_txt)
        ordi.features = ordi.features.iloc[:, :3]
        ordi.samples = ordi.samples.iloc[:, :3]
        ordi.eigvals = ordi.eigvals[:3]
        ordi.proportion_explained = ordi.proportion_explained[:3]
        ordi.write(biplot_txt)
    cmd += run_import(biplot_txt, biplot, "PCoAResults % Properties('biplot')")
    cmd += 'qiime emperor biplot \\\n'
    cmd += '--i-biplot %s \\\n' % biplot
    cmd += '--m-sample-metadata-file %s \\\n' % meta
    if isfile(biplot_tax):
        cmd += '--m-feature-metadata-file %s \\\n' % biplot_tax
    cmd += '--p-number-of-features 10 \\\n'
    cmd += '--o-visualization %s\n' % qzv
    # if isfile(biplot_tax):
    #     cmd += 'rm %s\n' % tax_tmp
    return cmd


def write_empress(
        qza: str,
        pcoa_biplot: str,
        qzv: str,
        meta: str,
        feat_metas: list,
        tree
) -> str:
    """

    Parameters
    ----------
    qza
    pcoa_biplot
    qzv
    meta
    feat_metas
    tree

    Returns
    -------

    """
    cmd = ''
    if '/biplot' in pcoa_biplot:
        biplot_txt = '%s.txt' % splitext(pcoa_biplot)[0]
        if isfile(biplot_txt):
            ordi = OrdinationResults.read(biplot_txt)
            ordi.features = ordi.features.iloc[:, :3]
            ordi.samples = ordi.samples.iloc[:, :3]
            ordi.eigvals = ordi.eigvals[:3]
            ordi.proportion_explained = ordi.proportion_explained[:3]
            ordi.write(biplot_txt)
        cmd += run_import(
            biplot_txt, pcoa_biplot, "PCoAResults % Properties('biplot')")

    cmd += 'qiime empress community-plot \\\n'
    cmd += '--i-tree %s \\\n' % tree[1]
    cmd += '--i-pcoa %s \\\n' % pcoa_biplot
    if tree[0]:
        cmd += '--i-feature-table %s \\\n' % tree[0]
    else:
        cmd += '--i-feature-table %s \\\n' % qza
    cmd += '--m-sample-metadata-file %s \\\n' % meta
    for feat_meta in feat_metas:
        cmd += '--m-feature-metadata-file %s \\\n' % feat_meta
    cmd += '--p-number-of-features 15 \\\n'
    cmd += '--p-filter-extra-samples \\\n'
    cmd += '--o-visualization %s\n' % qzv
    return cmd


def write_permanova_permdisp(
        meta: str,
        test: str,
        typ: str,
        dm: str,
        dm_filt: str,
        qzv: str,
        html: str
) -> str:
    """
    Determine whether groups of samples are significantly different from one
    another using a permutation-based statistical test.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/beta-group-significance/

    Parameters
    ----------
    meta
    test
    typ
    dm
    dm_filt
    qzv
    html

    Returns
    -------

    """
    cmd = ''
    if not isfile(dm_filt):
        cmd += 'qiime diversity filter-distance-matrix \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm
        cmd += '--m-metadata-file %s \\\n' % meta
        cmd += '--o-filtered-distance-matrix %s\n\n' % dm_filt
    if not isfile(qzv):
        cmd += 'qiime diversity beta-group-significance \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm_filt
        cmd += '--p-method %s \\\n' % typ
        cmd += '--m-metadata-file %s \\\n' % meta
        cmd += '--m-metadata-column "%s" \\\n' % test
        cmd += '--p-permutations 499 \\\n'
        cmd += '--o-visualization %s\n\n' % qzv
    if not isfile(html):
        cmd += run_export(qzv, html, 'perms')
    return cmd


def write_adonis(
        meta: str,
        formula: str,
        variables: list,
        stratas: list,
        dms_metrics: list,
        out: str,
        template: list
) -> list:
    """

    Parameters
    ----------
    meta
    formula
    variables
    stratas
    dms_metrics
    out
    template

    Returns
    -------

    """
    n_perm = '499'
    r_script = []
    for line in template:
        l = line
        if 'NAME_METRIC' in line:
            for (dm, _, me) in dms_metrics:
                if not isfile(dm):
                    continue
                if stratas:
                    for strata in stratas:
                        if 'PERMUTATIONS' in line:
                            l = block_line.replace(
                                'BLOCK', strata).replace(
                                'NAME', strata)
                        l = l.replace('FORMULA', formula)
                        l = l.replace(
                            'NAME', strata).replace(
                            'METRIC', me).replace(
                            'DM', me).replace(
                            'STRATA', strata).replace(
                            'VARS',
                            '", "'.join((variables + [strata]))).replace(
                            'PERMUTATIONS', ', permutations = perm')
                else:
                    l = l.replace(
                        'FORMULA', formula).replace(
                        'METRIC', me).replace(
                        'NAME', 'no_strata').replace(
                        'STRATA', 'no_strata').replace(
                        'VARS', '", "'.join(variables)).replace(
                        'PERMUTATIONS', '')
        elif 'METRIC' in line:
            for (dm, _, me) in dms_metrics:
                if not isfile(dm):
                    continue
                l = l.replace(
                    'METRIC', me).replace(
                    'DM_FP', '%s.tsv' % splitext(dm)[0])
        elif 'META_FP' in line:
            l = l.replace('META_FP', meta)
        elif 'NPERM' in line:
            l = l.replace('NPERM', n_perm)
        elif 'BLOCK' in line:
            block_line = line
        elif 'OUT' in line:
            l = l.replace('OUT', out)
        r_script.append(l)
    return r_script


def write_tsne(
        dm: str,
        dm_filt: str,
        meta: str,
        meta_met: str,
        group: str,
        tsne: str
) -> str:
    """

    Parameters
    ----------
    dm
    dm_filt
    meta
    meta_met
    group
    tsne

    Returns
    -------

    """
    cmd = ''
    if group != 'ALL':
        cmd += 'qiime diversity filter-distance-matrix \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm
        cmd += '--m-metadata-file %s \\\n' % meta_met
        cmd += '--o-filtered-distance-matrix %s\n' % dm_filt
    else:
        dm_filt = dm
    perplexity = '25'
    learning_rate = '200'
    early_exaggeration = '10'
    cmd += 'qiime diversity tsne \\\n'
    cmd += '--i-distance-matrix %s \\\n' % dm_filt
    cmd += '--p-perplexity %s \\\n' % perplexity
    cmd += '--p-learning-rate %s \\\n' % learning_rate
    cmd += '--p-early-exaggeration %s \\\n' % early_exaggeration
    cmd += '--o-tsne %s\n\n' % tsne
    if group != 'ALL':
        cmd += 'rm %s\n\n' % dm_filt
    cmd += 'mv %s %s.tsv\n\n' % (meta_met, meta)
    return cmd


def write_umap(
        dm: str,
        dm_filt: str,
        meta: str,
        meta_met: str,
        group: str,
        umap: str
) -> str:
    """

    Parameters
    ----------
    dm
    dm_filt
    meta
    meta_met
    group
    umap

    Returns
    -------

    """
    cmd = ''
    if group != 'ALL':
        cmd += 'qiime diversity filter-distance-matrix \\\n'
        cmd += '--i-distance-matrix %s \\\n' % dm
        cmd += '--m-metadata-file %s \\\n' % meta_met
        cmd += '--o-filtered-distance-matrix %s\n' % dm_filt
    else:
        dm_filt = dm
    n_neighbors = '15'
    min_dist = '0.4'
    cmd += 'qiime diversity umap \\\n'
    cmd += '--i-distance-matrix %s \\\n' % dm_filt
    cmd += '--p-n-neighbors %s \\\n' % n_neighbors
    cmd += '--p-min-dist %s \\\n' % min_dist
    cmd += '--o-umap %s\n\n' % umap
    if group != 'ALL':
        cmd += 'rm %s\n\n' % dm_filt
    cmd += 'mv %s %s.tsv\n\n' % (meta_met, meta)
    return cmd


def write_procrustes(
        meta_fp: str,
        meta_me: str,
        d1: str,
        d2: str,
        d1f: str,
        d2f: str,
        qzv: str,
        dis: str,
        tsv: str
) -> str:
    """

    Parameters
    ----------
    meta_fp
    meta_me
    d1
    d2
    d1f
    d2f
    qzv
    dis
    tsv

    Returns
    -------

    """
    pcoa_out1 = '%s_pcoa.qza' % splitext(d1f)[0]
    pcoa_out2 = '%s_pcoa.qza' % splitext(d2f)[0]
    ref_pcoa = '%s_ref.qza' % splitext(pcoa_out1)[0]
    oth_pcoa = '%s_oth.qza' % splitext(pcoa_out2)[0]
    cmd = '\nqiime diversity filter-distance-matrix \\\n'
    cmd += '--m-metadata-file %s \\\n' % meta_me
    cmd += '--i-distance-matrix %s \\\n' % d1
    cmd += '--o-filtered-distance-matrix %s\n' % d1f
    cmd += '\nqiime diversity filter-distance-matrix \\\n'
    cmd += '--m-metadata-file %s \\\n' % meta_me
    cmd += '--i-distance-matrix %s \\\n' % d2
    cmd += '--o-filtered-distance-matrix %s\n' % d2f
    cmd += '\nqiime diversity pcoa \\\n'
    cmd += '--i-distance-matrix %s \\\n' % d1f
    cmd += '--o-pcoa %s\n' % pcoa_out1
    cmd += '\nqiime diversity pcoa \\\n'
    cmd += '--i-distance-matrix %s \\\n' % d2f
    cmd += '--o-pcoa %s\n' % pcoa_out2
    cmd += '\nqiime diversity procrustes-analysis \\\n'
    cmd += '--i-reference %s \\\n' % pcoa_out1
    cmd += '--i-other %s \\\n' % pcoa_out2
    cmd += '--o-disparity-results %s \\\n' % dis
    cmd += '--o-transformed-reference %s \\\n' % ref_pcoa
    cmd += '--o-transformed-other %s\n' % oth_pcoa
    cmd += '\nqiime emperor procrustes-plot \\\n'
    cmd += '--i-reference-pcoa %s \\\n' % ref_pcoa
    cmd += '--i-other-pcoa %s \\\n' % oth_pcoa
    cmd += '--m-metadata-file %s \\\n' % meta_me
    cmd += '--o-visualization %s\n' % qzv
    cmd += 'mv %s %s.tsv\n' % (meta_me, meta_fp)
    cmd += run_export(dis, tsv, '')
    cmd += 'rm %s %s %s %s %s %s %s\n' % (
        pcoa_out1, pcoa_out2, ref_pcoa, oth_pcoa, d1f, d2f, dis)
    return cmd


def write_mantel(
        dat1: str,
        dat2: str,
        meta_fp: str,
        meta_met: str,
        d1: str,
        d2: str,
        d1f: str,
        d2f: str,
        qzv: str,
        html: str
) -> str:
    """

    Parameters
    ----------
    dat1
    dat2
    meta_fp
    meta_met
    d1
    d2
    d1f
    d2f
    qzv
    html

    Returns
    -------

    """
    cmd = '\nqiime diversity filter-distance-matrix \\\n'
    cmd += '--m-metadata-file %s \\\n' % meta_met
    cmd += '--i-distance-matrix %s \\\n' % d1
    cmd += '--o-filtered-distance-matrix %s\n' % d1f
    cmd += '\nqiime diversity filter-distance-matrix \\\n'
    cmd += '--m-metadata-file %s \\\n' % meta_met
    cmd += '--i-distance-matrix %s \\\n' % d2
    cmd += '--o-filtered-distance-matrix %s\n' % d2f
    cmd += '\nqiime diversity mantel \\\n'
    cmd += '--i-dm1 %s \\\n' % d1f
    cmd += '--i-dm2 %s \\\n' % d2f
    cmd += '--p-label1 %s \\\n' % dat1
    cmd += '--p-label2 %s \\\n' % dat2
    cmd += '--o-visualization %s\n' % qzv
    cmd += run_export(qzv, html, 'mantel')
    cmd += 'rm %s %s %s\n' % (qzv, d1f, d2f)
    cmd += 'mv %s %s.tsv\n' % (meta_met, meta_fp)
    return cmd


def write_phate(
        config,
        qza: str,
        new_qza: str,
        new_tsv: str,
        meta: str,
        html: str
) -> str:
    cmd = ''
    if not isfile(new_tsv):
        cmd += '\nqiime feature-table filter-samples \\\n'
        cmd += '--i-table %s \\\n' % qza
        cmd += '--m-metadata-file %s \\\n' % meta
        cmd += '--o-filtered-table %s\n' % new_qza
        cmd += run_export(new_qza, new_tsv, 'FeatureTable')
        cmd += 'rm %s %s\n' % (new_qza, new_qza.replace('.qza', '.biom'))

    cmd += '\nXphate \\\n'
    cmd += '--i-table %s \\\n' % new_tsv
    cmd += '--m-metadata %s \\\n' % meta
    cmd += '--o-html %s \\\n' % html
    if 'labels' in config.phate:
        for label in config.phate['labels']:
            cmd += '--p-labels %s \\\n' % label
    # cmd += '-fp %s \\\n' % fp
    # cmd += '-fa %s \\\n' % fa
    n_nodes = config.run_params['beta']['n_nodes']
    n_procs = config.run_params['beta']['n_procs']
    cmd += '--p-cpus %s \\\n' % (int(n_nodes) * int(n_procs))

    params = {'t': (None,), 'decay': (15,), 'knn': (5,)}
    for k, vs_ in params.items():
        if 'params' in config.phate and k in config.phate['params']:
            vs = config.phate['params'][k]
        else:
            vs = vs_
        for v in vs:
            if v:
                cmd += '--p-%ss %s \\\n' % (k, v)
    cmd += '--no-clusters \\\n'
    cmd += '--separate \\\n'
    cmd += '--verbose\n\n'
    return cmd


def write_sourcetracking(
        qza: str,
        new_qza: str,
        new_tsv: str,
        new_meta: str,
        meth: str,
        fp: str,
        fa: str,
        cur_rad: str,
        column: str,
        sink: str,
        sources: list,
        sourcetracking_params: dict,
        loo: bool,
        n_nodes: str,
        n_procs: str,
        imports: set
) -> str:
    """

    Parameters
    ----------
    qza
    new_qza
    new_tsv
    new_meta
    meth
    fp
    fa
    cur_rad
    column
    sink
    sources
    sourcetracking_params
    loo
    n_nodes
    n_procs
    imports

    Returns
    -------

    """
    cmd = ''
    if not isfile(new_tsv) and new_tsv not in imports:
        cmd = '\nqiime feature-table filter-samples'
        cmd += ' --i-table %s' % qza
        cmd += ' --m-metadata-file %s' % new_meta
        cmd += ' --o-filtered-table %s\n' % new_qza
        cmd = run_export(new_qza, new_tsv, 'FeatureTable')
        imports.add(new_tsv)

    cmd = '\nXsourcetracking'
    cmd += ' -i %s' % new_tsv
    cmd += ' -m %s' % new_meta
    cmd += ' -o %s' % cur_rad
    cmd += ' -c %s' % column
    cmd += ' -si %s' % sink
    for source in sources:
        if source:
            cmd += ' -so %s' % source
    cmd += ' -fp %s' % fp
    cmd += ' -fa %s' % fa
    cmd += ' -meth %s' % meth
    cmd += ' --p-cpus %s' % (int(n_nodes) * int(n_procs))
    if sourcetracking_params['rarefaction']:
        cmd += ' --p-rarefaction %s' % sourcetracking_params['rarefaction']
    if sourcetracking_params['iterations']:
        cmd += ' --p-iterations-burnins %s' % sourcetracking_params['iterations']
    if meth == 'sourcetracker' and loo:
        cmd += ' --loo \n'
    cmd += ' --verbose \n'
    return cmd


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


def write_mmvec(
        meta_fp: str,
        qza1: str,
        qza2: str,
        res_dir: str,
        model_odir: str,
        null_odir: str,
        ranks_tsv: str,
        ordination_tsv: str,
        stats: str,
        ranks_null_tsv: str,
        ordination_null_tsv: str,
        stats_null: str,
        summary: str,
        batch: str,
        learn: str,
        epoch: str,
        input_prior: str,
        output_prior: str,
        thresh_feat: str,
        latent_dim: str,
        train_column: str,
        n_example: str,
        summary_interval: str,
        gpu: bool,
        qiime_env: str
) -> str:
    """
    Performs bi-loglinear multinomial regression and calculates the
    conditional probability ranks of metabolite co-occurence given the microbe
    presence.

    Parameters
    ----------
    meta_fp
    qza1
    qza2
    res_dir
    model_odir
    null_odir
    ranks_tsv
    ordination_tsv
    stats
    ranks_null_tsv
    ordination_null_tsv
    stats_null
    summary
    batch
    learn
    epoch
    input_prior
    output_prior
    thresh_feat
    latent_dim
    train_column
    n_example
    summary_interval
    gpu
    qiime_env

    Returns
    -------

    """
    cmd = ''
    if gpu:
        biom1 = '%s.biom' % splitext(qza1)[0]
        biom2 = '%s.biom' % splitext(qza2)[0]
        cmd += '\nmmvec paired-omics \\\n'
        cmd += '--arm-the-gpu \\\n'
        cmd += '--microbe-file %s \\\n' % biom1
        cmd += '--metabolite-file %s \\\n' % biom2
        cmd += '--min-feature-count %s \\\n' % thresh_feat
        cmd += '--epochs %s \\\n' % epoch
        cmd += '--batch-size %s \\\n' % batch
        cmd += '--latent-dim %s \\\n' % latent_dim
        cmd += '--input-prior %s \\\n' % input_prior
        cmd += '--output-prior %s \\\n' % output_prior
        cmd += '--learning-rate %s \\\n' % learn
        cmd += '--beta1 0.85 \\\n'
        cmd += '--beta2 0.90 \\\n'
        cmd += '--checkpoint-interval %s \\\n' % summary_interval
        cmd += '--summary-interval %s \\\n' % summary_interval
        cmd += '--summary-dir %s \\\n' % res_dir
        cmd += '--ranks-file %s\n' % ranks_tsv
    else:
        ranks_qza = '%s.qza' % splitext(ranks_tsv)[0]
        ranks_null_qza = '%s.qza' % splitext(ranks_null_tsv)[0]
        ordination_qza = '%s.qza' % splitext(ordination_tsv)[0]
        ordination_null_qza = '%s.qza' % splitext(ordination_null_tsv)[0]
        summary_html = '%s.html' % splitext(summary)[0]
        if not isfile(ranks_qza) or not isfile(ordination_qza) or not isfile(stats):
            cmd += '\ncd %s\n' % model_odir
            cmd_mmvec = '\nqiime mmvec paired-omics \\\n'
            cmd_mmvec += '--i-microbes %s \\\n' % qza1
            cmd_mmvec += '--i-metabolites %s \\\n' % qza2
            cmd_mmvec += '--m-metadata-file %s \\\n' % meta_fp
            if str(train_column) != 'None':
                cmd_mmvec += '--p-training-column %s \\\n' % train_column
            else:
                cmd_mmvec += '--p-num-testing-examples %s \\\n' % n_example
            cmd_mmvec += '--p-min-feature-count %s \\\n' % thresh_feat
            cmd_mmvec += '--p-epochs %s \\\n' % epoch
            cmd_mmvec += '--p-batch-size %s \\\n' % batch
            cmd_mmvec += '--p-latent-dim %s \\\n' % latent_dim
            cmd_mmvec += '--p-input-prior %s \\\n' % input_prior
            cmd_mmvec += '--p-output-prior %s \\\n' % output_prior
            cmd_mmvec += '--p-learning-rate %s \\\n' % learn
            cmd_mmvec += '--p-summary-interval %s \\\n' % summary_interval
            if 'qiime2-2020' in qiime_env:
                cmd_mmvec += '--p-equalize-biplot \\\n'
            cmd_mmvec += '--o-conditionals %s \\\n' % ranks_qza
            cmd_mmvec += '--o-conditional-biplot %s \\\n' % ordination_qza
            cmd_mmvec += '--o-model-stats %s \\\n' % stats
            cmd_mmvec += '--output-dir %s/logdir\n' % model_odir
            cmd_mmvec += '\nrm -rf %s/logdir\n' % model_odir
            cmd += cmd_mmvec

            cmd += '\ncd %s\n' % null_odir
            cmd += cmd_mmvec.replace(
                '--p-latent-dim %s' % latent_dim,
                '--p-latent-dim 0'
            ).replace(
                '--o-conditionals %s' % ranks_qza,
                '--o-conditionals %s' % ranks_null_qza
            ).replace(
                '--o-conditional-biplot %s' % ordination_qza,
                '--o-conditional-biplot %s' % ordination_null_qza
            ).replace(
                '--o-model-stats %s' % stats,
                '--o-model-stats %s' % stats_null
            ).replace(
                '%s/logdir' % model_odir,
                '%s/logdir' % null_odir
            )

            cmd += '\nqiime mmvec summarize-paired \\\n'
            cmd += '--i-model-stats %s \\\n' % stats
            cmd += '--i-baseline-stats %s \\\n' % stats_null
            cmd += '--o-visualization %s\n' % summary
            cmd += run_export(summary, summary_html, 'mmvec_summary')

        if not isfile(ranks_tsv):
            cmd += run_export(ranks_qza, ranks_tsv, '')
        if not isfile(ordination_tsv):
            cmd += run_export(ordination_qza, ordination_tsv, 'mmvec')
    return cmd


def run_add_metadata(
        input_path: str,
        output_path: str,
        meta: str
) -> str:
    """

    Parameters
    ----------
    input_path
    output_path
    meta

    Returns
    -------

    """
    cmd = 'biom add-metadata \\\n'
    cmd += '  -i %s \\\n' % input_path
    cmd += '  -o %s \\\n' % output_path
    cmd += '  --sample-metadata-fp %s\n' % meta
    return cmd


def write_songbird(
        qza,
        new_qza,
        new_meta,
        nsams,
        params,
        formula,
        bformula,
        out_paths,
        force
) -> tuple:
    """

    Parameters
    ----------
    qza
    new_qza
    new_meta
    nsams
    params
    formula
    bformula
    out_paths
    force

    Returns
    -------

    """
    fcmd = ''
    if not isfile(new_qza):
        fcmd += filter_feature_table(qza, new_qza, new_meta)

    batches = int(params['batches'])

    cmd = ''
    if force or not isfile(out_paths['diff_qza']):
        cmd = '# model\n'
        cmd += '\nqiime songbird multinomial \\\n'
        cmd += ' --i-table %s \\\n' % new_qza
        cmd += ' --m-metadata-file %s \\\n' % new_meta
        cmd += ' --p-formula "%s" \\\n' % formula
        cmd += ' --p-epochs %s \\\n' % params['epochs']
        if batches > 0.8 * nsams:
            cmd += ' --p-batch-size %s \\\n' % str(int(nsams * 0.8))
        else:
            cmd += ' --p-batch-size %s \\\n' % params['batches']
        cmd += ' --p-differential-prior %s \\\n' % params['diff_priors']
        cmd += ' --p-learning-rate %s \\\n' % params['learns']
        cmd += ' --p-min-sample-count %s \\\n' % params['thresh_samples']
        cmd += ' --p-min-feature-count %s \\\n' % params['thresh_feats']
        if 'examples' in params:
            cmd += ' --p-num-random-test-examples %s \\\n' % params['examples']
        else:
            cmd += ' --p-training-column %s \\\n' % params['train']
        cmd += ' --p-summary-interval %s \\\n' % params['summary_interval']
        cmd += ' --o-differentials %s \\\n' % out_paths['diff_qza']
        cmd += ' --o-regression-stats %s \\\n' % out_paths['stat']
        cmd += ' --o-regression-biplot %s\n' % out_paths['plot']

    if force or not isfile(out_paths['diff']):
        cmd += run_export(out_paths['diff_qza'], out_paths['diff'], '')

    stat_tsv = '%s.txt' % splitext(out_paths['stat'])[0]
    if force or not isfile(stat_tsv):
        cmd += run_export(out_paths['stat'], stat_tsv, '')

    bcmd = ''
    if force or len(out_paths['bdiff_qza']) and not isfile(out_paths['bstat']):
        bcmd += '\nqiime songbird multinomial \\\n'
        bcmd += ' --i-table %s \\\n' % new_qza
        bcmd += ' --m-metadata-file %s \\\n' % new_meta
        bcmd += ' --p-formula "%s" \\\n' % bformula
        bcmd += ' --p-epochs %s \\\n' % params['epochs']
        if batches > 0.8 * nsams:
            bcmd += ' --p-batch-size %s \\\n' % str(int(nsams * 0.8))
        else:
            bcmd += ' --p-batch-size %s \\\n' % params['batches']
        bcmd += ' --p-differential-prior %s \\\n' % params['diff_priors']
        bcmd += ' --p-learning-rate %s \\\n' % params['learns']
        bcmd += ' --p-min-sample-count %s \\\n' % params['thresh_samples']
        bcmd += ' --p-min-feature-count %s \\\n' % params['thresh_feats']
        if 'examples' in params:
            bcmd += ' --p-num-random-test-examples %s \\\n' % params['examples']
        else:
            bcmd += ' --p-training-column %s \\\n' % params['train']
        bcmd += ' --p-summary-interval %s \\\n' % params['summary_interval']
        bcmd += ' --o-differentials %s \\\n' % out_paths['bdiff_qza']
        bcmd += ' --o-regression-stats %s \\\n' % out_paths['bstat']
        bcmd += ' --o-regression-biplot %s\n' % out_paths['bplot']

    if force or not isfile(out_paths['tens']):
        cmd += '\nqiime songbird summarize-paired \\\n'
        cmd += ' --i-regression-stats %s \\\n' % out_paths['stat']
        cmd += ' --i-baseline-stats %s \\\n' % out_paths['bstat']
        cmd += ' --o-visualization %s\n' % out_paths['tens']

    if force or not isfile(out_paths['html']):
        cmd += run_export(out_paths['tens'], out_paths['html'], 'songbird')

    return cmd, fcmd, bcmd


def get_biplot_commands(
        ordi_edit_fp,
        qza,
        qzv,
        omic_feature,
        omic_sample,
        meta1_fp,
        meta2_fp,
        n_edit,
        max_r
):

    cmd = '\n'
    if max_r >= 2:
        if not isfile(qza):
            cmd += '\nqiime tools import'
            cmd += ' --input-path %s' % ordi_edit_fp
            cmd += ' --output-path %s' % qza
            cmd += ' --type "PCoAResults %s Properties([\'biplot\'])"\nsleep 3' % '%'
        cmd += '\nqiime emperor biplot'
        cmd += ' --i-biplot %s' % qza
        cmd += ' --m-%s-metadata-file %s' % (omic_feature, meta1_fp)
        cmd += ' --m-%s-metadata-file %s' % (omic_sample, meta2_fp)
        cmd += ' --p-number-of-features %s' % n_edit
        cmd += ' --o-visualization %s\n' % qzv
    return cmd


def get_xmmvec_commands(
        ordi_edit_fp,
        omic1,
        omic2,
        meta1_fp,
        meta2_fp,
        xmmvecs,
        pair
):
    cmd = '\n'
    ranks_fp = ordi_edit_fp.replace('ordination.txt', 'ranks.tsv')
    nrows = None
    ncols = None
    if isfile(ranks_fp):
        nrows = pd.read_table(ranks_fp, usecols=[0, 1, 2]).shape[0]
        ncols = pd.read_table(ranks_fp, nrows=3, index_col=0).shape[1]

    ranks_html = ordi_edit_fp.replace('ordination.txt', 'ranks.html')
    # if not isfile(ranks_html):
    if 1:
        cmd += '\nXmmvec'
        cmd += ' --i-ranks-path %s' % ranks_fp
        cmd += ' --o-ranks-explored %s' % ranks_html
        cmd += ' --p-omic1-name %s' % omic1
        cmd += ' --p-omic2-name %s' % omic2
        if nrows > 50:
            cmd += ' --p-omic1-max 50'
        if ncols > 50:
            cmd += ' --p-omic2-max 50'
        if xmmvecs and pair in xmmvecs:
            cmd += ' --p-omic1-metadata %s' % meta1_fp
            cmd += ' --p-omic2-metadata %s' % meta2_fp
            if omic1 in xmmvecs[pair]:
                if 'color_variable' in xmmvecs[pair][omic1]:
                    cmd += ' --p-omic1-column %s' % xmmvecs[
                        pair][omic1]['color_variable']
            if omic2 in xmmvecs[pair]:
                if 'color_variable' in xmmvecs[pair][omic2]:
                    cmd += ' --p-omic2-column %s' % xmmvecs[
                        pair][omic2]['color_variable']
        cmd += '\n'
    return cmd


def get_paired_heatmaps_command(
        ranks_fp: str,
        omic1_common_fp: str,
        omic2_common_fp: str,
        taxonomy_tsv: str,
        features_names: list,
        topn: int,
        paired_heatmap_qzv: str,
        pre_paired_fp: str
):

    cmd = ''
    # if not isfile(paired_heatmap_qzv):
    if 1:

        omic1_tmp = '%s_tmp.tsv' % splitext(omic1_common_fp)[0]
        omic2_tmp = '%s_tmp.tsv' % splitext(omic2_common_fp)[0]
        omic1_qza_tmp = '%s_tmp.qza' % splitext(omic1_common_fp)[0]
        omic2_qza_tmp = '%s_tmp.qza' % splitext(omic2_common_fp)[0]
        taxonomy_tsv_tmp = '%s_tmp.tsv' % splitext(taxonomy_tsv)[0]
        ranks_fp_tmp = '%s_tmp.tsv' % splitext(ranks_fp)[0]
        ranks_qza_tmp = '%s_tmp.qza' % splitext(ranks_fp)[0]

        py = '%s.py' % splitext(paired_heatmap_qzv)[0]
        with open(py, 'w') as o, open(pre_paired_fp) as f:
            for line in f:
                if "'OMIC1_COMMON_FP_TMP'" in line:
                    o.write(line.replace('OMIC1_COMMON_FP_TMP', omic1_tmp))
                elif "'OMIC2_COMMON_FP_TMP'" in line:
                    o.write(line.replace('OMIC2_COMMON_FP_TMP', omic2_tmp))
                elif "'OMIC1_COMMON_QZA_TMP'" in line:
                    o.write(line.replace('OMIC1_COMMON_QZA_TMP', omic1_qza_tmp))
                elif "'OMIC2_COMMON_QZA_TMP'" in line:
                    o.write(line.replace('OMIC2_COMMON_QZA_TMP', omic2_qza_tmp))
                elif "'OMIC1_COMMON_FP'" in line:
                    o.write(line.replace('OMIC1_COMMON_FP', omic1_common_fp))
                elif "'OMIC2_COMMON_FP'" in line:
                    o.write(line.replace('OMIC2_COMMON_FP', omic2_common_fp))
                elif "'TAXONOMY_TSV_TMP'" in line:
                    o.write(line.replace('TAXONOMY_TSV_TMP', taxonomy_tsv_tmp))
                elif "'TAXONOMY_TSV'" in line:
                    o.write(line.replace('TAXONOMY_TSV', taxonomy_tsv))
                elif "'RANKS_FP_TMP'" in line:
                    o.write(line.replace('RANKS_FP_TMP', ranks_fp_tmp))
                elif "'RANKS_QZA_TMP'" in line:
                    o.write(line.replace('RANKS_QZA_TMP', ranks_qza_tmp))
                elif "'RANKS_FP'" in line:
                    o.write(line.replace('RANKS_FP', ranks_fp))
                else:
                    o.write(line)

        cmd += '\npython3 %s\n' % py

        cmd += '\nqiime mmvec paired-heatmap'
        cmd += ' --i-ranks %s' % ranks_qza_tmp
        cmd += ' --i-microbes-table %s' % omic1_qza_tmp
        cmd += ' --i-metabolites-table %s' % omic2_qza_tmp
        cmd += ' --m-microbe-metadata-file %s' % taxonomy_tsv_tmp
        cmd += ' --m-microbe-metadata-column Taxon'
        if features_names:
            cmd += ' --p-top-k-microbes 0'
            for features_name in features_names:
                cmd += ' --p-features %s' % features_name
        else:
            cmd += ' --p-top-k-microbes %s' % topn
        cmd += ' --p-normalize rel_row'
        cmd += ' --p-top-k-metabolites 100'
        cmd += ' --p-level 6'
        cmd += ' --o-visualization %s\n' % paired_heatmap_qzv

        cmd += '\nrm %s %s %s %s\n' % (
            ranks_qza_tmp, omic1_qza_tmp, omic2_qza_tmp, taxonomy_tsv_tmp)
    return cmd
