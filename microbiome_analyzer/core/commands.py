# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
from os.path import basename, dirname, isfile, splitext
from skbio.stats.ordination import OrdinationResults

from microbiome_analyzer._scratch import io_update, to_do, rep


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
            cmd += 'biom convert'
            cmd += ' -i %s' % input_path
            cmd += ' -o %s' % infile
            cmd += ' --table-type="OTU table"'
            cmd += ' --to-hdf5\n'
        output_path_tmp = '%s.tmp.qza' % output_path
        cmd += 'qiime tools import'
        cmd += ' --input-path %s' % infile
        cmd += ' --output-path %s' % output_path_tmp
        cmd += ' --type "FeatureTable[Frequency]"\n'
        cmd += 'qiime feature-table filter-samples'
        cmd += ' --i-table %s' % output_path_tmp
        cmd += ' --p-min-frequency 1'
        cmd += ' --p-filter-empty-features'
        cmd += ' --o-filtered-table %s\n' % output_path
        cmd += 'rm %s\n' % output_path_tmp
    else:
        cmd += 'qiime tools import'
        cmd += ' --input-path %s' % infile
        cmd += ' --output-path %s' % output_path
        cmd += ' --type "%s"\n' % typ
    return cmd


def run_export(
        input_path: str,
        output_path: str,
        typ: str = ''
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

    cmd = 'qiime tools export'
    cmd += '  --input-path %s' % input_path
    cmd += '  --output-path %s\n' % out
    if typ.startswith("FeatureTable"):
        if output_path.endswith('biom'):
            cmd += 'mv %s/*.biom %s\n' % (out, output_path)
            cmd += 'rm -rf %s\n' % inp
        else:
            cur_biom = '%s.biom' % out
            cmd += 'mv %s/*.biom %s\n' % (out, cur_biom)
            cmd += 'biom convert'
            cmd += '  -i %s' % cur_biom
            cmd += '  -o %s.tmp' % output_path
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
    cmd = 'qiime feature-table rarefy'
    cmd += ' --i-table %s' % qza
    cmd += ' --p-sampling-depth %s' % str(depth)
    cmd += ' --o-rarefied-table %s\n' % qza_out
    return cmd


def write_fasta(
        seqs_fasta: str,
        seqs_qza: str,
        biom_table: biom.Table,
        tsv_fp: str = ''
) -> str:
    """
    Write the fasta sequences.

    Parameters
    ----------
    seqs_fasta
    seqs_qza
    biom_table
    tsv_fp

    Returns
    -------

    """
    with open(rep(seqs_fasta), 'w') as fas_o:
        for seq in biom_table.ids(axis='observation'):
            fas_o.write('>%s\n%s\n' % (seq.strip(), seq.strip()))
    cmd = '# Write features as fasta file:\n'
    cmd += '#  - Features from: %s\n' % tsv_fp
    cmd += '# Snippet:\n'
    cmd += '# ```:\n'
    cmd += "# with open(fasta_out, 'w') as o:\n"
    cmd += "#     for seq in tsv_pd.index:\n"
    cmd += "#         o.write('>%s\\n%s\\n' % (seq.strip(), seq.strip()))\n"
    cmd += '# ```:\n'
    cmd += run_import(seqs_fasta, seqs_qza, 'FeatureData[Sequence]')
    return cmd


def write_taxonomy_sklearn(
        classifier_qza: str,
        seqs_qza: str,
        out_qza: str
) -> str:
    """
    Classify reads by taxon using a fitted classifier.
    https://docs.qiime2.org/2020.2/plugins/available/feature-classifier/classify-sklearn

    Parameters
    ----------
    classifier_qza
    seqs_qza
    out_qza

    Returns
    -------

    """
    cmd = 'qiime feature-classifier classify-sklearn'
    cmd += ' --i-classifier %s' % classifier_qza
    cmd += ' --i-reads %s' % seqs_qza
    cmd += ' --p-n-jobs %s' % '4'
    cmd += ' --o-classification %s\n' % out_qza
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
        cmd += 'qiime taxa collapse'
        cmd += ' --i-table %s' % tab_qza
        cmd += ' --i-taxonomy %s' % tax_qza
        cmd += ' --p-level %s' % level
        cmd += ' --o-collapsed-table %s\n\n' % collapsed_qza
    if remove_empty:
        tax_tmp = '%s_filtempty%s.tsv' % (splitext(tax_qza)[0], level)
        cmd += '\n# Metadata to remove empty/catch-all features: %s\n' % (
            tax_tmp)
        with open(tax_tmp, 'w') as o:
            o.write('Feature ID\tTaxon\n')
            for tax in remove_empty:
                cmd += '#  - %s\n' % tax
                o.write('%s\tremove\n' % tax)
        cmd += '\n\nqiime feature-table filter-features'
        cmd += ' --i-table %s' % collapsed_qza
        cmd += ' --m-metadata-file %s' % tax_tmp
        cmd += ' --o-filtered-table %s2.qza' % collapsed_qza
        cmd += ' --p-exclude-ids\n'
        cmd += 'mv %s2.qza %s\n' % (collapsed_qza, collapsed_qza)
    if not isfile(collapsed_tsv):
        cmd += run_export(collapsed_qza, collapsed_tsv, 'FeatureTable')
    if not isfile(collapsed_meta):
        cmd += 'cp %s %s\n' % (meta_fp, collapsed_meta)
    return cmd


def write_collapse(
        self,
        dat: str,
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
    self
    dat : str
    tab_qza : str
    tax_qza : str
    collapsed_qza : str
    collapsed_tsv : str
    meta_fp : str
    collapsed_meta : str
    level : int
    remove_empty : set

    Returns
    -------

    """
    cmd = ''
    if to_do(collapsed_qza):
        cmd += 'qiime taxa collapse'
        cmd += ' --i-table %s' % tab_qza
        cmd += ' --i-taxonomy %s' % tax_qza
        cmd += ' --p-level %s' % level
        cmd += ' --o-collapsed-table %s\n\n' % collapsed_qza
        io_update(self, i_f=[tab_qza, tax_qza], o_f=collapsed_qza, key=dat)

    if remove_empty:
        tax_tmp = '%s_filtempty%s.tsv' % (splitext(tax_qza)[0], level)
        cmd += '\n# Metadata to remove empty/catch-all features: %s\n' % (
            tax_tmp)
        with open(rep(tax_tmp), 'w') as o:
            o.write('Feature ID\tTaxon\n')
            for tax in remove_empty:
                cmd += '#  - %s\n' % tax
                o.write('%s\tremove\n' % tax)
        cmd += '\n\nqiime feature-table filter-features'
        cmd += ' --i-table %s' % collapsed_qza
        cmd += ' --m-metadata-file %s' % tax_tmp
        cmd += ' --o-filtered-table %s2.qza' % collapsed_qza
        cmd += ' --p-exclude-ids\n'
        cmd += 'mv %s2.qza %s\n' % (collapsed_qza, collapsed_qza)
        io_update(self, i_f=tax_tmp, o_f=collapsed_qza, key=dat)

    if to_do(collapsed_tsv):
        cmd += run_export(collapsed_qza, collapsed_tsv, 'FeatureTable')
        if isfile(collapsed_qza):
            io_update(self, i_f=collapsed_qza, o_f=collapsed_tsv, key=dat)
        else:
            io_update(self, o_f=collapsed_tsv, key=dat)

    if to_do(collapsed_meta):
        cmd += 'cp %s %s\n' % (meta_fp, collapsed_meta)
        io_update(self, i_f=meta_fp, o_f=collapsed_meta, key=dat)

    return cmd


def write_sepp(
        self,
        dat: str,
        seqs_qza: str,
        tree_qza: str,
        sepp_tree: str,
        qza: str,
        qza_in: str,
        tsv_in: str,
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
    self
    dat
    seqs_qza
    tree_qza
    o_tree
    qza
    qza_in
    tsv_in
    qza_out

    Returns
    -------

    """
    plac = '%s/plac_%s' % (dirname(sepp_tree), basename(sepp_tree))
    cmd = ''
    if to_do(o_tree):
        cmd += 'qiime fragment-insertion sepp'
        cmd += ' --i-representative-sequences %s' % seqs_qza
        cmd += ' --i-reference-database %s' % tree_qza
        cmd += ' --o-tree %s' % o_tree
        cmd += ' --o-placements %s' % plac
        cmd += ' --p-threads %s\n' % self.run_params['sepp']['cpus']
        io_update(self, i_f=[seqs_qza, tree_qza], o_f=[o_tree, plac], key=dat)
    if to_do(tsv_in):
        cmd += 'qiime fragment-insertion filter-features'
        cmd += ' --i-table %s' % qza
        cmd += ' --i-tree %s' % o_tree
        cmd += ' --o-filtered-table %s' % qza_in
        cmd += ' --o-removed-table %s\n' % qza_out
        cmd += run_export(qza_in, tsv_in, 'FeatureTable')
        i_f = [qza]
        if not to_do(o_tree):
            i_f.append(o_tree)
        io_update(self, i_f=i_f, o_f=[qza_in, qza_out, tsv_in], key=dat)
    return cmd


def write_alpha(
        self,
        dat: str,
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
    self
    dat
    qza
    qza_out
    datasets_phylo
    trees
    metric

    Returns
    -------

    """
    i_f = []
    if metric in ['faith_pd']:
        if not datasets_phylo[0] or not trees or to_do(trees[1]):
            return ''
        cmd = 'qiime diversity alpha-phylogenetic'
        if datasets_phylo[1]:
            cmd += ' --i-table %s' % trees[0]
            i_f.append(trees[0])
        else:
            cmd += ' --i-table %s' % qza
            i_f.append(qza)
        cmd += ' --i-phylogeny %s' % trees[1]
        i_f.append(trees[1])
        cmd += ' --p-metric %s' % metric
        cmd += ' --o-alpha-diversity %s\n' % qza_out
    else:
        i_f.append(qza)
        cmd = 'qiime diversity alpha'
        cmd += ' --i-table %s' % qza
        cmd += ' --p-metric %s' % metric
        cmd += ' --o-alpha-diversity %s\n\n' % qza_out
    io_update(self, i_f=i_f, o_f=qza_out, key=dat)
    return cmd


def write_sample_filter(
        qza: str,
        qza_subset: str,
        meta_subset: str
) -> str:
    """
    filter-samples: Filter samples from table
    https://docs.qiime2.org/2023.2/plugins/available/feature-table/filter-samples/

    Filter samples from table based on frequency and/or metadata. Any features
    with a frequency of zero after sample filtering will also be removed. See
    the filtering tutorial on https://docs.qiime2.org for additional details.


    Parameters
    ----------
    qza : str
        Path to the raw data table
    qza_subset : str
        Path to the filtered data table
    meta_subset : str
        Metadata table used for sample filering

    Returns
    -------
    cmd : str
        Sample filering command line
    """
    cmd = 'qiime feature-table filter-samples'
    cmd += ' --i-table %s' % qza
    cmd += ' --m-metadata-file %s' % meta_subset
    cmd += ' --o-filtered-table %s\n\n' % qza_subset
    return cmd


def write_feat_filter(
        qza: str,
        qza_subset: str,
        meta_subset: str
) -> str:
    """
    filter-features: Filter features from table
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

    cmd = 'qiime feature-table filter-features'
    cmd += ' --i-table %s' % qza
    cmd += ' --m-metadata-file %s' % meta_subset
    cmd += ' --o-filtered-table %s\n\n' % qza_subset
    return cmd


def write_barplots(
        qza: str,
        tax_qza: str,
        meta: str,
        qzv: str
) -> str:
    """
    Visualize taxonomy with an interactive bar plot
    https://docs.qiime2.org/2020.2/plugins/available/taxa/barplot/

    This visualizer produces an interactive barplot visualization of
    taxonomies. Interactive features include multi-level sorting, plot
    recoloring, sample relabeling, and SVG figure export.

    Parameters
    ----------
    qza
    tax_qza
    meta
    qzv

    Returns
    -------

    """
    cmd = 'qiime taxa barplot'
    cmd += ' --i-table %s' % qza
    cmd += ' --i-taxonomy %s' % tax_qza
    cmd += ' --m-metadata-file %s' % meta
    cmd += ' --o-visualization %s' % qzv
    return cmd


def write_krona(
        qza: str,
        tax_qza: str,
        meta: str,
        qzv: str
) -> str:
    """
    Plugin for creating Krona plots
    https://library.qiime2.org/plugins/q2-krona/39/

    q2-krona is developed to make is easy to generate Krona plots, because
    the tool needs some rearrangement on FeatureTable[Frequency] to be able
    to work.

    Parameters
    ----------
    qza : str
        Path to the input data qza file
    tax_qza : str
        Path to the input taxonomy qza file
    meta : str
        Path to the metadata
    qzv : str
        Path to the krona plot qzv

    Returns
    -------
    cmd : str
        Command line to generate the krona plot
    """
    cmd = 'qiime krona collapse-and-plot'
    cmd += ' --i-table %s' % qza
    cmd += ' --i-taxonomy %s' % tax_qza
    cmd += ' --o-krona-plot %s' % qzv
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
    cmd = 'qiime metadata tabulate'
    cmd += ' --o-visualization %s' % out_qza
    for alpha in alphas:
        cmd += ' --m-input-file %s' % alpha[0]
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
    cmd = 'qiime diversity alpha-correlation'
    cmd += ' --i-alpha-diversity %s' % qza
    cmd += ' --p-method %s' % method
    cmd += ' --m-metadata-file %s' % meta
    cmd += ' --o-visualization %s\n\n' % qzv
    return cmd


def write_alpha_rarefaction(
        self,
        dat: str,
        qza: str,
        qzv: str,
        metric: str,
        raref: str,
        data
) -> str:
    """

    Parameters
    ----------
    self
    dat
    qza
    qzv
    metric
    reref
    data

    Returns
    -------

    """
    phylo = data.phylo
    tree = data.tree

    i_f = [data.meta]
    cmd = 'qiime diversity alpha-rarefaction'
    if metric in ['faith_pd']:
        if not phylo[0] or not tree[-1]:
            return ''
        if phylo[1]:
            cmd += ' --i-table %s' % tree[0]
            i_f.append(tree[0])
        else:
            cmd += ' --i-table %s' % qza
            i_f.append(qza)
        cmd += ' --i-phylogeny %s' % tree[1]
        i_f.append(tree[1])
    else:
        cmd += ' --i-table %s' % qza
        i_f.append(qza)
    cmd += ' --m-metadata-file %s' % data.meta
    cmd += ' --p-metrics %s' % metric
    if raref:
        cmd += ' --p-max-depth %s' % raref.split('raref')[-1]
    else:
        cmd += ' --p-max-depth 5000'
    cmd += ' --o-visualization %s\n\n' % qzv
    io_update(self, i_f=i_f, o_f=qzv, key=dat)
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
    cmd = 'qiime longitudinal volatility'
    cmd += ' --m-metadata-file %s' % meta
    cmd += ' --p-state-column "%s"' % timepoint
    cmd += ' --p-individual-id-column "%s"' % host
    cmd += ' --o-visualization %s\n' % qzv
    return cmd


def write_beta(
        self,
        dat: str,
        qza: str,
        dm_qza: str,
        metric: str,
        data
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
    self
    dat : str
    qza
    dm_qza
    metric
    data

    Returns
    -------

    """
    phylo = data.phylo
    tree = data.tree

    i_f = []
    cmd = ''
    if 'unifrac' in metric:
        if phylo[0] and tree[1] and not to_do(tree[1]):
            cmd += 'qiime diversity beta-phylogenetic'
            if tree[0]:
                cmd += ' --i-table %s' % tree[0]
                i_f.append(tree[0])
            else:
                cmd += ' --i-table %s' % qza
                i_f.append(qza)
            cmd += ' --i-phylogeny %s' % tree[1]
            i_f.append(tree[1])
    else:
        cmd = 'qiime diversity beta'
        cmd += ' --i-table %s' % qza
        i_f.append(qza)
    if cmd:
        cmd += ' --p-metric %s' % metric
        nodes = self.config.run_params['beta']['nodes']
        cpus = self.config.run_params['beta']['cpus']
        qiime_env = self.config.qiime_env
        if float(qiime_env.split('-')[1]) >= 2020.8 and 'phylo' in cmd:
            cmd += ' --p-threads %s' % (int(nodes) * int(cpus))
        else:
            cmd += ' --p-n-jobs %s' % (int(nodes) * int(cpus))
        cmd += ' --o-distance-matrix %s\n' % dm_qza
        io_update(self, i_f=i_f, o_f=dm_qza, key=dat)
    return cmd


def write_rpca(
        self,
        dat: str,
        qza: str,
        meta: str,
        new_qza: str,
        ordi: str,
        dist: str,
        tree_qza: str,
        table: str,
        taxon: str,
        qzv: str,
        qurro: str,
        data
) -> str:
    """
    Performs Phylogenetic Robust Aitchison PCA.
    https://github.com/biocore/gemelli

    Parameters
    ----------
    self
    dat : str
    qza : str
    meta : str
    new_qza : str
    ordi : str
    dist : str
    tree_qza : str
    table : str
    taxon : str
    qzv : str
    qurro : str
    data

    Returns
    -------
    cmd : str
    """
    phylo = data.phylo
    tree = data.tree
    tax = data.tax
    taxa = data.taxa
    i_f, o_f = [], []
    cmd, cmd_exp = '', ''
    is_phylo = False
    qza_in = qza
    if phylo[0] and tree[1] and not to_do(tree[1]):
        is_phylo = True
        cmd += 'qiime gemelli phylogenetic-rpca-with-taxonomy'
        if tree[0]:
            qza_in = tree[0]
            i_f.append(tree[0])
        else:
            i_f.append(qza)
        cmd += ' --i-table %s' % new_qza
        cmd += ' --i-phylogeny %s' % tree[1]
        cmd += ' --m-taxonomy-file %s' % tax[1]
        cmd += ' --o-counts-by-node-tree %s' % tree_qza
        cmd += ' --o-counts-by-node %s' % table
        cmd += ' --o-t2t-taxonomy %s' % taxon
        tree_nwk = '%s.nwk' % splitext(tree_qza)[0]
        cmd_exp += run_export(tree_qza, tree_nwk, 'phylogeny')
        i_f.extend([tree[1], tax[1]])
        o_f.extend([tree_qza, tree_nwk, table, taxon])
    else:
        cmd += 'qiime gemelli auto-rpca'
        cmd += ' --i-table %s' % new_qza
        i_f.append(qza)
    if to_do(ordi):
        cmd += ' --p-min-feature-count %s' % self.config.run_params['rpca'][
            'min_feature_count']
        cmd += ' --p-min-sample-count %s' % self.config.run_params['rpca'][
            'min_sample_count']
        cmd += ' --o-biplot %s' % ordi
        cmd += ' --o-distance-matrix %s\n' % dist
        o_f.extend([ordi, dist])
        io_update(self, i_f=i_f, o_f=o_f, key=dat)
    else:
        io_update(self, i_f=i_f, o_f=o_f, key=dat)
        cmd = ''

    if self.config.jobs:
        cmd += 'rm %s %s\n\n' % (meta, new_qza)

    # export distance matrix
    dist_tsv = dist.replace('.qza', '.tsv')
    cmd += run_export(dist, dist_tsv)
    # export ordination
    ordi_tsv = ordi.replace('.qza', '.txt')
    cmd += run_export(ordi, ordi_tsv, 'pcoa')
    if cmd_exp:
        cmd += cmd_exp

    rm_cmd, feat_cmd = '', ''
    if data.feat_meta:
        nid_tsv = ordi_tsv + '.nID.txt'
        cmd += 'grep "^n[0-9]*\\t" %s | cut -f1 > %s\n' % (ordi_tsv, nid_tsv)
        for feat_meta in data.feat_meta:
            cmd += 'cat %s %s > %s.nID.tsv\n' % (feat_meta, nid_tsv, feat_meta)
            rm_cmd += 'rm %s.nID.tsv\n' % feat_meta
            feat_cmd += ' --m-feature-metadata-file %s.nID.tsv' % feat_meta
    if is_phylo:
        cmd += 'qiime empress community-plot'
        cmd += ' --i-tree %s' % tree_qza
        cmd += ' --i-feature-table %s' % table
        cmd += ' --i-pcoa %s' % ordi
        cmd += ' --m-feature-metadata-file %s' % taxon
        cmd += ' --p-filter-missing-features'
        cmd += ' --p-ignore-missing-samples'
    else:
        cmd += 'qiime emperor biplot'
        cmd += ' --i-biplot %s' % ordi
        cmd += ' --m-feature-metadata-file %s' % tax[1]
        cmd += ' --m-feature-metadata-file %s' % taxa[-1]

    if data.feat_meta:
        cmd += feat_cmd
    cmd += ' --p-number-of-features %s' % self.config.run_params['rpca'][
        'number_of_features']
    cmd += ' --m-sample-metadata-file %s' % data.meta
    cmd += ' --o-visualization %s\n\n' % qzv

    if to_do(ordi):
        io_update(self, i_f=meta, o_f=[qzv, ordi_tsv], key=dat)
    else:
        io_update(self, i_f=[ordi, meta], o_f=[qzv, ordi_tsv], key=dat)

    cmd_final = ''
    if self.config.jobs:
        cmd_final += 'rm -rf $TMPDIR\n'
    cmd_final += 'qiime feature-table filter-samples'
    cmd_final += ' --i-table %s' % qza_in
    cmd_final += ' --p-min-frequency 1'
    cmd_final += ' --p-filter-empty-features'
    cmd_final += ' --m-metadata-file %s' % meta
    cmd_final += ' --o-filtered-table %s\n\n' % new_qza
    cmd_final += cmd + rm_cmd

    cmd_final += 'qiime qurro loading-plot'
    cmd_final += ' --i-ranks %s' % ordi
    if phylo[0] and tree[1] and not to_do(tree[1]):
        cmd_final += ' --i-table %s' % table
        cmd_final += ' --m-feature-metadata-file %s' % taxon
    else:
        cmd_final += ' --i-table %s' % new_qza
        cmd_final += ' --m-feature-metadata-file %s' % tax[1]
    cmd_final += ' --m-sample-metadata-file %s' % data.meta
    cmd_final += ' --o-visualization %s\n' % qurro

    io_update(self, i_f=[qza_in, meta], o_f=qurro, key=dat)

    return cmd_final


def write_deicode(
        self,
        dat: str,
        qza: str,
        meta: str,
        new_qza: str,
        ordi: str,
        ordi_tsv: str,
        dm_qza: str,
        qzv: str
) -> str:
    """
    Performs robust center log-ratio transform robust PCA and
    ranks the features by the loadings of the resulting SVD.
    https://library.qiime2.org/plugins/deicode/19/

    Parameters
    ----------
    self
    dat
    qza
    meta
    new_qza
    ordi
    ordi_tsv
    dm_qza
    qzv

    Returns
    -------

    """
    cmd = ''
    if to_do(new_qza):
        cmd += 'qiime feature-table filter-samples'
        cmd += ' --i-table %s' % qza
        cmd += ' --p-min-frequency 1'
        cmd += ' --p-filter-empty-features'
        cmd += ' --m-metadata-file %s' % meta
        cmd += ' --o-filtered-table %s\n\n' % new_qza
        io_update(self, i_f=[qza, meta], o_f=new_qza, key=dat)

    if to_do(ordi) or to_do(dm_qza):
        cmd += 'qiime deicode rpca'
        cmd += ' --i-table %s' % new_qza
        # cmd += ' --p-min-feature-count 10'
        # cmd += ' --p-min-sample-count 500'
        cmd += ' --p-n-components 2'
        cmd += ' --o-biplot %s' % ordi
        cmd += ' --o-distance-matrix %s\n\n' % dm_qza
        if not to_do(new_qza):
            io_update(self, i_f=new_qza, o_f=[ordi, dm_qza], key=dat)
        else:
            io_update(self, o_f=[ordi, dm_qza], key=dat)

    cmd += 'qiime emperor biplot'
    cmd += ' --i-biplot %s' % ordi
    cmd += ' --m-sample-metadata-file %s' % meta
    cmd += ' --o-visualization %s' % qzv
    cmd += ' --p-number-of-features 10\n\n'
    cmd += 'rm %s %s\n\n' % (meta, new_qza)
    cmd += run_export(ordi, ordi_tsv, 'pcoa')
    if to_do(ordi):
        io_update(self, i_f=meta, o_f=[qzv, ordi_tsv], key=dat)
    else:
        io_update(self, i_f=[ordi, meta], o_f=[qzv, ordi_tsv], key=dat)

    return cmd


def write_pcoa(
        self,
        dat: str,
        dm: str,
        dm_filt: str,
        meta: str,
        meta_met: str,
        group: str,
        pcoa: str,
        pcoa_tsv: str
) -> str:
    """
    Apply principal coordinate analysis.
    https://docs.qiime2.org/2019.10/plugins/available/diversity/pcoa/

    Parameters
    ----------
    self
    dat : str
    dm
    dm_filt
    meta
    meta_met
    group
    pcoa
    pcoa_tsv

    Returns
    -------

    """
    cmd = ''
    if group != 'ALL':
        cmd += 'qiime diversity filter-distance-matrix'
        cmd += ' --i-distance-matrix %s' % dm
        cmd += ' --m-metadata-file %s' % meta_met
        cmd += ' --o-filtered-distance-matrix %s\n' % dm_filt
    else:
        dm_filt = dm

    cmd += 'qiime diversity pcoa'
    cmd += ' --i-distance-matrix %s' % dm_filt
    cmd += ' --o-pcoa %s\n\n' % pcoa
    if group != 'ALL':
        cmd += 'rm %s\n\n' % dm_filt
    cmd += 'mv %s %s.tsv\n\n' % (meta_met, meta)
    cmd += run_export(pcoa, pcoa_tsv, 'pcoa')

    io_update(self, i_f=[dm, meta_met], o_f=[pcoa, pcoa_tsv, meta], key=dat)

    return cmd


def write_biplot(
        self,
        dat:str,
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
    self
    dat
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
    # if not to_do(tax):
    #     tax_dict = {}
    #     with open(rep(tax)) as f, open(rep(biplot_tax), 'w') as o_tax:
    #         o_tax.write('Feature ID\tTaxon\tPrevious ID\n')
    #         n = 0
    #         for ldx, line in enumerate(f):
    #             line_split = line.strip().split('\t')
    #             if len(line_split) == 1:
    #                 line_split = [line_split[0], line_split[0]]
    #             if ldx and not line.startswith('#q2:types'):
    #                 new = 'x__%s;%s' % (n, line_split[1])
    #                 tax_dict[line_split[0]] = new
    #                 o_tax.write('%s\t%s\t%s\n' % (new, new, line_split[0]))
    #                 n += 1
    #
    #     tab_tsv = '%s_table.tsv' % splitext(biplot)[0]
    #     with open(rep(tsv)) as f, open(rep(tab_tsv), 'w') as o_tab:
    #         for ldx, line in enumerate(f):
    #             t = line.strip().split('\t')
    #             if t[0] in tax_dict:
    #                 o_tab.write('%s\t%s\n' % (tax_dict[t[0]], '\t'.join(t[1:])))
    #             else:
    #                 o_tab.write(line)
    #     tab_qza = '%s.qza' % splitext(tab_tsv)[0]
    #     cmd += run_import(tab_tsv, tab_qza, 'FeatureTable[Frequency]')
    #     io_update(self, i_f=tab_tsv, key=dat)
    # else:
    #     tab_qza = qza
    #     io_update(self, i_f=qza, key=dat)

    tab_rel_qza_tmp = '%s_rel_tmp.qza' % splitext(qza)[0]
    tab_rel_qza = '%s_rel.qza' % splitext(qza)[0]
    # tab_rel_qza_tmp = '%s_rel_tmp.qza' % splitext(tab_qza)[0]
    # tab_rel_qza = '%s_rel.qza' % splitext(tab_qza)[0]

    cmd += '\nqiime feature-table relative-frequency'
    # cmd += ' --i-table %s' % tab_qza
    cmd += ' --i-table %s' % qza
    cmd += ' --o-relative-frequency-table %s\n\n' % tab_rel_qza_tmp
    cmd += '\nqiime feature-table filter-samples'
    cmd += ' --i-table %s' % tab_rel_qza_tmp
    cmd += ' --p-filter-empty-features'
    cmd += ' --m-metadata-file %s.tsv' % meta
    cmd += ' --o-filtered-table %s\n\n' % tab_rel_qza

    cmd += 'qiime diversity pcoa-biplot'
    cmd += ' --i-pcoa %s' % pcoa
    cmd += ' --i-features %s' % tab_rel_qza
    cmd += ' --o-biplot %s\n\n' % biplot

    cmd += 'rm %s %s\n' % (tab_rel_qza_tmp, tab_rel_qza)

    out_biplot_txt = '%s.txt' % splitext(biplot)[0]
    cmd += run_export(biplot, out_biplot_txt, 'biplot')

    io_update(
        self, i_f=[meta, qza, pcoa], o_f=[biplot, out_biplot_txt], key=dat)

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
    cmd = 'qiime emperor plot'
    cmd += ' --i-pcoa %s' % pcoa
    cmd += ' --m-metadata-file %s' % meta
    cmd += ' --o-visualization %s\n\n' % qzv
    return cmd


def write_emperor_biplot(
        self,
        dat: str,
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
    self
    dat
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
    if not to_do(biplot_txt):
        ordi = OrdinationResults.read(rep(biplot_txt))
        ordi.features = ordi.features.iloc[:, :3]
        ordi.samples = ordi.samples.iloc[:, :3]
        ordi.eigvals = ordi.eigvals[:3]
        ordi.proportion_explained = ordi.proportion_explained[:3]
        ordi.write(rep(biplot_txt))
    cmd += run_import(biplot_txt, biplot, "PCoAResults % Properties('biplot')")
    cmd += 'qiime emperor biplot'
    cmd += ' --i-biplot %s' % biplot
    cmd += ' --m-sample-metadata-file %s' % meta
    i_f = [biplot_txt, meta]
    if not to_do(biplot_tax):
        cmd += ' --m-feature-metadata-file %s' % biplot_tax
        i_f.append(biplot_tax)
    cmd += ' --p-number-of-features 10'
    cmd += ' --o-visualization %s\n' % qzv
    # if isfile(biplot_tax):
    #     cmd += 'rm %s\n' % tax_tmp
    io_update(self, i_f=i_f, o_f=qzv, key=dat)

    return cmd


def write_empress(
        self,
        dat: str,
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
    self
    dat
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

    i_f = [tree[1], meta]
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
        i_f.append(biplot_txt)
    else:
        i_f.append(pcoa_biplot)

    cmd += 'qiime empress community-plot'
    cmd += ' --i-tree %s' % tree[1]
    cmd += ' --i-pcoa %s' % pcoa_biplot
    if tree[0]:
        cmd += ' --i-feature-table %s' % tree[0]
        i_f.append(tree[0])
    else:
        cmd += ' --i-feature-table %s' % qza
        i_f.append(qza)
    cmd += ' --m-sample-metadata-file %s' % meta
    for feat_meta in feat_metas:
        cmd += ' --m-feature-metadata-file %s' % feat_meta
        i_f.append(feat_meta)
    cmd += ' --p-number-of-features 15'
    cmd += ' --p-filter-extra-samples'
    cmd += ' --o-visualization %s\n' % qzv
    io_update(self, i_f=i_f, o_f=qzv, key=dat)
    return cmd


def write_permanova_permdisp(
        self,
        dat: str,
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
    self
    dat
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
    if to_do(dm_filt):
        cmd += 'qiime diversity filter-distance-matrix'
        cmd += ' --i-distance-matrix %s' % dm
        cmd += ' --m-metadata-file %s' % meta
        cmd += ' --o-filtered-distance-matrix %s\n\n' % dm_filt
        io_update(self, i_f=[dm, meta], o_f=dm_filt, key=dat)
    else:
        io_update(self, i_f=dm_filt, key=dat)

    if to_do(qzv) or to_do(html):
        cmd += 'qiime diversity beta-group-significance'
        cmd += ' --i-distance-matrix %s' % dm_filt
        cmd += ' --p-method %s' % typ
        cmd += ' --m-metadata-file %s' % meta
        cmd += ' --m-metadata-column "%s"' % test
        cmd += ' --p-permutations 499'
        cmd += ' --o-visualization %s\n\n' % qzv
        cmd += run_export(qzv, html, 'perms')
        i_f = [meta]
        if isfile(dm_filt):
            i_f.append(dm_filt)
        io_update(self, i_f=i_f, o_f=[qzv, html], key=dat)

    return cmd


def write_adonis(
        self,
        dat: str,
        meta_fp: str,
        formula: str,
        variables: list,
        stratas: list,
        dms_metrics: list,
        out: str,
) -> list:
    """

    Parameters
    ----------
    self
    dat
    meta_fp
    formula
    variables
    stratas
    dms_metrics
    out

    Returns
    -------
    r : list
        R scripts line
    """
    def adonis_collector(n, s, me, r):
        r.append("res[[%s]]$strata <- '%s'" % (n, s))
        r.append("res[[%s]]$metric <- '%s'" % (n, me))
        r.append("res[[%s]]$factors <- row.names(res[[%s]])" % (n, n))
        r.append("row.names(res[[%s]]) <- NULL" % n)

    i_f = [meta_fp]
    if not self.config.jobs:
        dms_metrics = [(rep(x), y, z) for x, y, z in dms_metrics]
        meta_fp = rep(meta_fp)
        out = rep(out)

    n_perm = '499'
    r = ["library(vegan)", "library(readr)", "library(data.table)"]
    r.append("res <- list()")
    r.append("meta_fp <- '%s'" % meta_fp)
    r.append("meta <- read.table(meta_fp, header=TRUE, check.names=FALSE, "
             "sep='\\t', colClasses=c('sample_name'='character'))")
    r.append("row.names(meta) <- meta[,'sample_name']")
    # r.append("meta <- meta[,-1]")
    meta = "meta[, c('%s'" % "', '".join(variables)
    if stratas:
        meta += ", '%s')" % "', '".join(stratas)
    else:
        meta += ")"
    r.append("meta <- %s, drop=FALSE]" % meta)
    r.append("meta <- na.omit(meta)")
    r.append("perm <- how(nperm = %s)" % n_perm)
    for (dm, _, me) in dms_metrics:
        if to_do(dm):
            continue
        r.append("dm_fp <- '%s.tsv'" % splitext(dm)[0])
        r.append("dm <- read.table(dm_fp, check.names = FALSE)")
        i_f.append('%s.tsv' % splitext(dm)[0])
        r.append("sdm <- dm[row.names(meta), row.names(meta)]")
        a2 = "as.data.frame(adonis2(sdm ~ %s, data=meta" % (formula)
        if stratas:
            for s in stratas:
                n = "'%s_%s'" % (s, me)
                r.append("setBlocks(perm) <- with(meta, meta[,'%s'])" % s)
                r.append("res[[%s]] <- %s, permutations = perm))" % (n, a2))
                adonis_collector(n, s, me, r)
        else:
            n = "'no_strata_%s'" % me
            r.append("res[[%s]] <- %s))" % (n, a2))
            adonis_collector(n, 'no_strata', me, r)
    r.append("out_fp <- '%s'" % out)
    r.append("write.table(x=rbindlist(res), file=out_fp, quote=FALSE, "
             "sep='\\t', row.names=FALSE)")
    io_update(self, i_f=i_f, o_f=out, key=dat)
    return r


def write_tsne(
        self,
        dat: str,
        dm: str,
        dm_filt: str,
        meta: str,
        meta_met: str,
        group: str,
        tsne: str,
        tsne_tsv: str
) -> str:
    """

    Parameters
    ----------
    self
    dat : str
    dm
    dm_filt
    meta
    meta_met
    group
    tsne
    tsne_tsv

    Returns
    -------

    """
    cmd = ''
    if group != 'ALL':
        cmd += 'qiime diversity filter-distance-matrix'
        cmd += ' --i-distance-matrix %s' % dm
        cmd += ' --m-metadata-file %s' % meta_met
        cmd += ' --o-filtered-distance-matrix %s\n' % dm_filt
    else:
        dm_filt = dm

    perplexity = '25'
    learning_rate = '200'
    early_exaggeration = '10'
    cmd += 'qiime diversity tsne'
    cmd += ' --i-distance-matrix %s' % dm_filt
    cmd += ' --p-perplexity %s' % perplexity
    cmd += ' --p-learning-rate %s' % learning_rate
    cmd += ' --p-early-exaggeration %s' % early_exaggeration
    cmd += ' --o-tsne %s\n\n' % tsne
    if group != 'ALL':
        cmd += 'rm %s\n\n' % dm_filt
    cmd += 'mv %s %s.tsv\n\n' % (meta_met, meta)
    cmd += run_export(tsne, tsne_tsv, 'pcoa')

    io_update(self, i_f=[dm, meta_met], o_f=[tsne, tsne_tsv, meta], key=dat)

    return cmd


def write_umap(
        self,
        dat: str,
        dm: str,
        dm_filt: str,
        meta: str,
        meta_met: str,
        group: str,
        umap: str,
        umap_tsv: str
) -> str:
    """

    Parameters
    ----------
    self
    dat : str
    dm
    dm_filt
    meta
    meta_met
    group
    umap
    umap_tsv

    Returns
    -------

    """
    cmd = ''
    if group != 'ALL':
        cmd += 'qiime diversity filter-distance-matrix'
        cmd += ' --i-distance-matrix %s' % dm
        cmd += ' --m-metadata-file %s' % meta_met
        cmd += ' --o-filtered-distance-matrix %s\n' % dm_filt
    else:
        dm_filt = dm

    n_neighbors = '15'
    min_dist = '0.4'
    cmd += 'qiime diversity umap'
    cmd += ' --i-distance-matrix %s' % dm_filt
    cmd += ' --p-n-neighbors %s' % n_neighbors
    cmd += ' --p-min-dist %s' % min_dist
    cmd += ' --o-umap %s\n\n' % umap
    if group != 'ALL':
        cmd += 'rm %s\n\n' % dm_filt
    cmd += 'mv %s %s.tsv\n\n' % (meta_met, meta)
    cmd += run_export(umap, umap_tsv, 'pcoa')

    io_update(self, i_f=[dm, meta_met], o_f=[umap, umap_tsv, meta], key=dat)

    return cmd


def write_procrustes(
        self,
        dat1: str,
        dat2: str,
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
    self
    dat1 : str
    dat2 : str
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
    cmd = '\nqiime diversity filter-distance-matrix'
    cmd += ' --m-metadata-file %s' % meta_me
    cmd += ' --i-distance-matrix %s' % d1
    cmd += ' --o-filtered-distance-matrix %s\n' % d1f
    cmd += '\nqiime diversity filter-distance-matrix'
    cmd += ' --m-metadata-file %s' % meta_me
    cmd += ' --i-distance-matrix %s' % d2
    cmd += ' --o-filtered-distance-matrix %s\n' % d2f
    cmd += '\nqiime diversity pcoa'
    cmd += ' --i-distance-matrix %s' % d1f
    cmd += ' --o-pcoa %s\n' % pcoa_out1
    cmd += '\nqiime diversity pcoa'
    cmd += ' --i-distance-matrix %s' % d2f
    cmd += ' --o-pcoa %s\n' % pcoa_out2
    cmd += '\nqiime diversity procrustes-analysis'
    cmd += ' --i-reference %s' % pcoa_out1
    cmd += ' --i-other %s' % pcoa_out2
    cmd += ' --o-disparity-results %s' % dis
    cmd += ' --o-transformed-reference %s' % ref_pcoa
    cmd += ' --o-transformed-other %s\n' % oth_pcoa
    cmd += '\nqiime emperor procrustes-plot'
    cmd += ' --i-reference-pcoa %s' % ref_pcoa
    cmd += ' --i-other-pcoa %s' % oth_pcoa
    cmd += ' --m-metadata-file %s' % meta_me
    cmd += ' --o-visualization %s\n' % qzv
    meta_out = '%s.tsv' % meta_fp
    cmd += 'mv %s %s\n' % (meta_me, meta_out)
    cmd += run_export(dis, tsv, '')
    cmd += 'rm %s\n' % pcoa_out1
    cmd += 'rm %s\n' % pcoa_out2
    cmd += 'rm %s\n' % ref_pcoa
    cmd += 'rm %s\n' % oth_pcoa
    cmd += 'rm %s\n' % d1f
    cmd += 'rm %s\n' % d2f
    cmd += 'rm %s\n' % dis
    dat = '%s__%s' % (dat1, dat2)
    io_update(self, i_f=[meta_me, d1, d2], o_f=[meta_out, qzv, tsv], key=dat)

    return cmd


def write_mantel(
        self,
        dat1: str,
        r1: str,
        dat2: str,
        r2: str,
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
    self
    dat1 : str
    r1 : str
    dat2 : str
    r2 : str
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
    cmd = '\nqiime diversity filter-distance-matrix'
    cmd += ' --m-metadata-file %s' % meta_met
    cmd += ' --i-distance-matrix %s' % d1
    cmd += ' --o-filtered-distance-matrix %s\n' % d1f
    cmd += '\nqiime diversity filter-distance-matrix'
    cmd += ' --m-metadata-file %s' % meta_met
    cmd += ' --i-distance-matrix %s' % d2
    cmd += ' --o-filtered-distance-matrix %s\n' % d2f
    cmd += '\nqiime diversity mantel'
    cmd += ' --i-dm1 %s' % d1f
    cmd += ' --i-dm2 %s' % d2f
    cmd += ' --p-intersect-ids'
    cmd += ' --p-label1 %s' % (dat1 + r1)
    cmd += ' --p-label2 %s' % (dat2 + r2)
    cmd += ' --o-visualization %s\n' % qzv
    cmd += run_export(qzv, html, 'mantel')
    cmd += 'rm %s %s %s\n' % (qzv, d1f, d2f)
    meta_out = '%s.tsv' % meta_fp
    cmd += 'mv %s %s' % (meta_met, meta_out)

    dat = '%s__%s' % (dat1, dat2)
    io_update(self, i_f=[meta_met, d1, d2], o_f=[meta_out, qzv, html], key=dat)

    return cmd


def write_phate(
        self,
        dat: str,
        qza: str,
        new_qza: str,
        new_tsv: str,
        meta: str,
        html: str
) -> str:
    cmd = ''
    if to_do(new_tsv):
        cmd += '\nqiime feature-table filter-samples'
        cmd += ' --i-table %s' % qza
        cmd += ' --m-metadata-file %s' % meta
        cmd += ' --o-filtered-table %s\n' % new_qza
        cmd += run_export(new_qza, new_tsv, 'FeatureTable')
        cmd += 'rm %s %s\n' % (new_qza, new_qza.replace('.qza', '.biom'))
        io_update(self, i_f=qza, key=dat)

    cmd += '\nXphate'
    cmd += ' --i-table %s' % new_tsv
    cmd += ' --m-metadata %s' % meta
    cmd += ' --o-html %s' % html
    if 'labels' in self.config.phate:
        for label in self.config.phate['labels']:
            cmd += ' --p-labels %s' % label
    # cmd += '-fp %s' % fp
    # cmd += '-fa %s' % fa
    nodes = self.config.run_params['beta']['nodes']
    cpus = self.config.run_params['beta']['cpus']
    cmd += ' --p-cpus %s' % (int(nodes) * int(cpus))

    params = {'t': (None,), 'decay': (15,), 'knn': (5,)}
    for k, vs_ in params.items():
        if 'params' in self.config.phate and k in self.config.phate['params']:
            vs = self.config.phate['params'][k]
        else:
            vs = vs_
        for v in vs:
            if v:
                cmd += ' --p-%ss %s' % (k, v)
    cmd += ' --no-clusters'
    cmd += ' --separate'
    cmd += ' --verbose\n\n'
    i_f = [meta]
    if not to_do(new_tsv):
        i_f.append(new_tsv)
    io_update(self, i_f=i_f, o_d=dirname(html), key=dat)
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
        nodes: str,
        cpus: str,
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
    nodes
    cpus
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
    cmd += ' --p-cpus %s' % (int(nodes) * int(cpus))
    if sourcetracking_params['rarefaction']:
        cmd += ' --p-rarefaction %s' % sourcetracking_params['rarefaction']
    if sourcetracking_params['iterations']:
        cmd += ' --p-iterations-burnins %s' % sourcetracking_params['iterations']
    if meth == 'sourcetracker' and loo:
        cmd += ' --loo \n'
    cmd += ' --verbose \n'
    cmd += run_export(pcoa, pcoa_tsv, 'pcoa')
    return cmd


def filter_feature_table(qza: str, new_qza: str, meta: str) -> str:
    """
    :param qza:
    :param new_qza:
    :param meta:
    :return:
    """
    cmd = '\nqiime feature-table filter-samples'
    cmd += ' --i-table %s' % qza
    cmd += ' --m-metadata-file %s' % meta
    cmd += ' --o-filtered-table %s\n' % new_qza
    return cmd


def write_mmvec(
        self,
        pair: str,
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
    self
    pair
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
        cmd += '\nmmvec paired-omics'
        cmd += ' --arm-the-gpu'
        cmd += ' --microbe-file %s' % biom1
        cmd += ' --metabolite-file %s' % biom2
        cmd += ' --min-feature-count %s' % thresh_feat
        cmd += ' --epochs %s' % epoch
        cmd += ' --batch-size %s' % batch
        cmd += ' --latent-dim %s' % latent_dim
        cmd += ' --input-prior %s' % input_prior
        cmd += ' --output-prior %s' % output_prior
        cmd += ' --learning-rate %s' % learn
        cmd += ' --beta1 0.85'
        cmd += ' --beta2 0.90'
        cmd += ' --checkpoint-interval %s' % summary_interval
        cmd += ' --summary-interval %s' % summary_interval
        cmd += ' --summary-dir %s' % res_dir
        cmd += ' --ranks-file %s\n' % ranks_tsv
        io_update(
            self, i_f=[biom1, biom2], o_d=res_dir, o_f=ranks_tsv, key=pair)
    else:
        ranks_qza = '%s.qza' % splitext(ranks_tsv)[0]
        ranks_null_qza = '%s.qza' % splitext(ranks_null_tsv)[0]
        ordination_qza = '%s.qza' % splitext(ordination_tsv)[0]
        ordination_null_qza = '%s.qza' % splitext(ordination_null_tsv)[0]
        summary_html = '%s.html' % splitext(summary)[0]
        if to_do(ranks_qza) or to_do(ordination_qza) or to_do(stats):
            cmd += '\ncd %s\n' % model_odir
            cmd_mmvec = '\nqiime mmvec paired-omics'
            cmd_mmvec += ' --i-microbes %s' % qza1
            cmd_mmvec += ' --i-metabolites %s' % qza2
            cmd_mmvec += ' --m-metadata-file %s' % meta_fp
            if str(train_column) != 'None':
                cmd_mmvec += ' --p-training-column %s' % train_column
            else:
                cmd_mmvec += ' --p-num-testing-examples %s' % n_example
            cmd_mmvec += ' --p-min-feature-count %s' % thresh_feat
            cmd_mmvec += ' --p-epochs %s' % epoch
            cmd_mmvec += ' --p-batch-size %s' % batch
            cmd_mmvec += ' --p-latent-dim %s' % latent_dim
            cmd_mmvec += ' --p-input-prior %s' % input_prior
            cmd_mmvec += ' --p-output-prior %s' % output_prior
            cmd_mmvec += ' --p-learning-rate %s' % learn
            cmd_mmvec += ' --p-summary-interval %s' % summary_interval
            if 'qiime2-2020' in qiime_env:
                cmd_mmvec += ' --p-equalize-biplot'
            cmd_mmvec += ' --o-conditionals %s' % ranks_qza
            cmd_mmvec += ' --o-conditional-biplot %s' % ordination_qza
            cmd_mmvec += ' --o-model-stats %s' % stats
            cmd_mmvec += ' --output-dir %s/logdir\n' % model_odir
            cmd_mmvec += '\nrm -rf %s/logdir\n' % model_odir
            cmd_mmvec += run_export(ranks_qza, ranks_tsv, '')
            cmd_mmvec += run_export(ordination_qza, ordination_tsv, 'mmvec')
            cmd += cmd_mmvec

            cmd += '\ncd %s\n' % null_odir
            cmd += cmd_mmvec.replace(
                ' --p-latent-dim %s' % latent_dim,
                ' --p-latent-dim 0'
            ).replace(
                ' --o-conditionals %s' % ranks_qza,
                ' --o-conditionals %s' % ranks_null_qza
            ).replace(
                ' --o-conditional-biplot %s' % ordination_qza,
                ' --o-conditional-biplot %s' % ordination_null_qza
            ).replace(
                ' --o-model-stats %s' % stats,
                ' --o-model-stats %s' % stats_null
            ).replace(
                '%s/logdir' % model_odir,
                '%s/logdir' % null_odir
            )

            cmd += '\nqiime mmvec summarize-paired'
            cmd += ' --i-model-stats %s' % stats
            cmd += ' --i-baseline-stats %s' % stats_null
            cmd += ' --o-visualization %s\n' % summary
            cmd += run_export(summary, summary_html, 'mmvec_summary')

            io_update(
                self, i_f=[qza1, qza2, meta_fp], key=pair,
                o_f=[ranks_qza, ranks_tsv, ordination_qza, ordination_tsv,
                     stats, ranks_null_qza,  ordination_null_qza, stats_null,
                     summary, summary_html])
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
    cmd = 'biom add-metadata'
    cmd += '  -i %s' % input_path
    cmd += '  -o %s' % output_path
    cmd += ' --sample-metadata-fp %s\n' % meta
    return cmd


def write_songbird(
        self,
        dat,
        qza,
        new_qza,
        new_meta,
        nsams,
        params,
        formula,
        bformula,
        out_paths
) -> tuple:
    """

    Parameters
    ----------
    self
    dat
    qza
    new_qza
    new_meta
    nsams
    params
    formula
    bformula
    out_paths

    Returns
    -------

    """
    fcmd = ''
    if to_do(new_qza):
        fcmd += filter_feature_table(qza, new_qza, new_meta)
        io_update(self, i_f=[qza, new_meta], o_f=new_qza, key=(dat, 'f'))

    cmd = ''
    force = self.config.force
    stat = out_paths['stat']
    diff_qza = out_paths['diff_qza']
    if force or to_do(diff_qza):

        diff_tsv = out_paths['diff']
        stat_tsv = '%s.txt' % splitext(stat)[0]
        plot = out_paths['plot']

        cmd = '# model\n'
        cmd += '\nqiime songbird multinomial'
        cmd += ' --i-table %s' % new_qza
        cmd += ' --m-metadata-file %s' % new_meta
        cmd += ' --p-formula "%s"' % formula
        cmd += ' --p-epochs %s' % params['epochs']
        if int(params['batches']) > 0.8 * nsams:
            cmd += ' --p-batch-size %s' % str(int(nsams * 0.8))
        else:
            cmd += ' --p-batch-size %s' % params['batches']
        cmd += ' --p-differential-prior %s' % params['diff_priors']
        cmd += ' --p-learning-rate %s' % params['learns']
        cmd += ' --p-min-sample-count %s' % params['thresh_samples']
        cmd += ' --p-min-feature-count %s' % params['thresh_feats']
        if 'examples' in params:
            cmd += ' --p-num-random-test-examples %s' % params['examples']
        else:
            cmd += ' --p-training-column %s' % params['train']
        cmd += ' --p-summary-interval %s' % params['summary_interval']
        cmd += ' --o-differentials %s' % diff_qza
        cmd += ' --o-regression-stats %s' % stat
        cmd += ' --o-regression-biplot %s\n' % plot
        cmd += run_export(diff_qza, diff_tsv, '')
        cmd += run_export(stat, stat_tsv, '')
        io_update(self, i_f=[new_qza, new_meta], key=(dat, ''),
                  o_f=[diff_qza, diff_tsv, stat, stat_tsv, plot])

    bcmd = ''
    bdiff_qza = out_paths['bdiff_qza']
    bstat = out_paths['bstat']
    if force or len(bdiff_qza) and to_do(bstat):
        bplot = out_paths['bplot']
        bcmd += '\nqiime songbird multinomial'
        bcmd += ' --i-table %s' % new_qza
        bcmd += ' --m-metadata-file %s' % new_meta
        bcmd += ' --p-formula "%s"' % bformula
        bcmd += ' --p-epochs %s' % params['epochs']
        if int(params['batches']) > 0.8 * nsams:
            bcmd += ' --p-batch-size %s' % str(int(nsams * 0.8))
        else:
            bcmd += ' --p-batch-size %s' % params['batches']
        bcmd += ' --p-differential-prior %s' % params['diff_priors']
        bcmd += ' --p-learning-rate %s' % params['learns']
        bcmd += ' --p-min-sample-count %s' % params['thresh_samples']
        bcmd += ' --p-min-feature-count %s' % params['thresh_feats']
        if 'examples' in params:
            bcmd += ' --p-num-random-test-examples %s' % params['examples']
        else:
            bcmd += ' --p-training-column %s' % params['train']
        bcmd += ' --p-summary-interval %s' % params['summary_interval']
        bcmd += ' --o-differentials %s' % bdiff_qza
        bcmd += ' --o-regression-stats %s' % bstat
        bcmd += ' --o-regression-biplot %s\n' % bplot
        io_update(self, i_f=[new_qza, new_meta], o_f=[bdiff_qza, bstat, bplot],
                  key=(dat, 'b'))

    tens = out_paths['tens']
    html = out_paths['html']
    if force or to_do(tens):
        cmd += '\nqiime songbird summarize-paired'
        cmd += ' --i-regression-stats %s' % stat
        cmd += ' --i-baseline-stats %s' % bstat
        cmd += ' --o-visualization %s\n' % tens
        cmd += run_export(tens, html, 'songbird')
        io_update(self, i_f=[stat, bstat], o_f=[tens, html], key=(dat, ''))

    return cmd, fcmd, bcmd


def get_biplot_commands(
        self,
        pair,
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
        if to_do(qza):
            cmd += '\nqiime tools import'
            cmd += ' --input-path %s' % ordi_edit_fp
            cmd += ' --output-path %s' % qza
            cmd += ' --type "PCoAResults %s Properties([\'biplot\'])"\nsleep 3' % '%'
            io_update(self, i_f=ordi_edit_fp, o_f=qza, key=pair)
        else:
            io_update(self, i_f=qza, key=pair)
        cmd += '\nqiime emperor biplot'
        cmd += ' --i-biplot %s' % qza
        cmd += ' --m-%s-metadata-file %s' % (omic_feature, meta1_fp)
        cmd += ' --m-%s-metadata-file %s' % (omic_sample, meta2_fp)
        cmd += ' --p-number-of-features %s' % n_edit
        cmd += ' --o-visualization %s\n' % qzv
        io_update(self, i_f=[meta1_fp, meta2_fp], o_f=qzv, key=pair)
    return cmd


def get_xmmvec_commands(
        self,
        pair,
        ordi_edit_fp,
        omic1,
        omic2,
        meta1_fp,
        meta2_fp,
):
    cmd = '\n'
    nrows = None
    ncols = None
    ranks_fp = ordi_edit_fp.replace('ordination.txt', 'ranks.tsv')
    if not to_do(ranks_fp):
        nrows = pd.read_table(rep(ranks_fp), usecols=[0, 1, 2]).shape[0]
        ncols = pd.read_table(rep(ranks_fp), nrows=3, index_col=0).shape[1]

    ranks_html = ordi_edit_fp.replace('ordination.txt', 'ranks.html')
    # if self.config.force or not isfile(ranks_html):
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
        if self.xmmvecs and pair in self.xmmvecs:
            cmd += ' --p-omic1-metadata %s' % meta1_fp
            cmd += ' --p-omic2-metadata %s' % meta2_fp
            if omic1 in self.xmmvecs[pair]:
                if 'color_variable' in self.xmmvecs[pair][omic1]:
                    cmd += ' --p-omic1-column %s' % self.xmmvecs[
                        pair][omic1]['color_variable']
            if omic2 in self.xmmvecs[pair]:
                if 'color_variable' in self.xmmvecs[pair][omic2]:
                    cmd += ' --p-omic2-column %s' % self.xmmvecs[
                        pair][omic2]['color_variable']
        cmd += '\n'
        io_update(self, i_f=[ranks_fp, meta1_fp, meta2_fp], o_f=ranks_html,
                  key=pair)
    return cmd
