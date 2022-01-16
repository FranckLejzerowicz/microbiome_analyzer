# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from os.path import isfile, splitext

import pandas as pd
from biom import load_table

from microbiome_analyzer._io import (
    read_meta_pd, convert_to_biom, check_features, get_output, get_cohort)
from microbiome_analyzer._rarefy import get_dat_depths, get_digit_depth
from microbiome_analyzer._taxonomy import (
    get_tax_tables, parse_split_taxonomy, edit_split_taxonomy)


class Data(object):
    """
    """
    wols = []
    dnas = []
    names = []

    def __init__(self, dataset):
        self.dat = dataset
        Data.names.append(dataset)
        self.tsv = {}
        self.biom = {}
        self.qza = {}
        self.meta = ''
        self.data = {}
        self.metadata = pd.DataFrame()
        self.feat_meta = []
        self.sb = None
        self.phylo = ('', 0)
        self.features = {}
        self.rarefs = ['']
        self.raref_depths = []
        self.tax = None
        self.taxa = None
        self.tree = ('', '', '')
        self.seqs = None
        self.collapsed = {}
        self.alpha = {}
        self.alphas = []
        self.subsets = {}
        self.sample_subsets = {}
        self.rpca = {}
        self.beta = {}
        self.perms = {}
        self.adonis = {}
        self.r2 = {}
        self.pcoa = {}
        self.tsne = {}
        self.umap = {}
        self.phate = {}
        self.sourcetracking = {}
        self.biplot = {}
        self.doc = {}
        self.filts = None
        self.filt = None
        self.taxon = None
        self.source = dataset
        self.subset = ''
        self.path = None

    def read_biom(self, index=''):
        biom = load_table(self.biom[index])
        if biom.shape[0] >= 10:
            self.data[index] = biom

    def read_tsv(self):
        data_pd = read_meta_pd(self.tsv[''], '#OTU ID').set_index('#OTU ID')
        self.data[''] = data_pd

    def read_meta_pd(self):
        meta_pd = read_meta_pd(self.meta, 'sample_name')
        self.metadata = meta_pd

    def check_gid_or_dna(self):
        gids, dna, correct = check_features(self.data[''])
        if len(gids) == self.data[''].shape[0]:
            Data.wols.append(self.dat)
            self.features = gids
            if correct:
                id_map = dict((x, x.replace(r'[; ]+', '|'))
                              for x in self.data[''].ids(axis='observation'))
                self.data[''].update_ids(id_map, axis='observation')
                self.phylo = ('wol', 1)
            else:
                self.phylo = ('wol', 0)
        if dna:
            Data.dnas.append(self.dat)
            self.phylo = ('amplicon', 0)

    def get_feature_metadata(self):
        feat = self.meta.split('metadata/meta_')[0]
        feat_fps = glob.glob('%s/qiime/feature_metadata/feat_*' % feat)
        if feat_fps:
            for feat_fp in feat_fps:
                # feature metadata can be large: no formatting check, yet!
                # idx = 'Feature ID'
                # feat_pd = read_meta_pd(feat_fp, idx).set_index(idx)
                self.feat_meta.append(feat_fp)


class Datasets(object):
    """Collect the data associated with each dataset passed but the user
    """
    raw_filt = {}
    filt_raw = {}
    raw_coll = {}
    coll_raw = {}

    def __init__(self, config) -> None:
        """Initialize the class instance with the dataset name"""
        self.config = config
        self.datasets = {}
        self.collect_datasets()
        self.edits = {}

    def collect_datasets(self):
        for dat in self.config.datasets:
            data = Data(dat)
            tsv = '%s/data/tab_%s.tsv' % (self.config.folder, dat)
            biom = '%s.biom' % splitext(tsv)[0]
            if not isfile(biom) and not isfile(tsv):
                print('[skipping] Not tsv/biom table for %s' % dat)
                continue
            meta = '%s/metadata/meta_%s.tsv' % (self.config.folder, dat)
            if not isfile(meta):
                print('[skipping] Not metadata table for %s' % dat)
                continue
            data.meta = meta
            data.tsv[''] = tsv
            if isfile(biom):
                data.biom[''] = biom
                data.read_biom()
            else:
                data.read_tsv()
                data.data[''] = convert_to_biom(data.data[''])
            data.read_meta_pd()
            data.check_gid_or_dna()
            data.get_feature_metadata()
            self.datasets[dat] = data

    def set_rarefaction_paths(self):
        for dataset, data in self.datasets.items():
            if dataset in Datasets.filt_raw:
                data.raref_depths = self.datasets[
                    Datasets.filt_raw[dataset]
                ].raref_depths
            if not data.raref_depths:
                continue
            for depth_ in data.raref_depths[1]:
                depth = '_raref%s' % get_digit_depth(
                    depth_, data.data[''].sum())
                data.rarefs.append(depth)

    def set_rarefaction_depths(self):
        for dataset, data in self.datasets.items():
            sam_sum = pd.Series(data.data[''].sum(axis='sample'))
            skip, depths = get_dat_depths(
                dataset, self.config.folder, self.config.rarefs, sam_sum)
            if skip:
                continue
            data.raref_depths = depths

    def set_rarefaction(self):
        self.set_rarefaction_depths()
        self.set_rarefaction_paths()

    def set_taxonomy_paths(self, method):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            odir = get_output(self.config.folder, 'taxonomy/%s' % dataset)
            if data.phylo and data.phylo[0] == 'amplicon':
                tax_tsv = '%s/tax_%s_%s.tsv' % (odir, dataset, method)
                meth = method
            else:
                tax_tsv = '%s/tax_%s.tsv' % (odir, dataset)
                if data.phylo and data.phylo[0] == 'wol':
                    meth = 'wol'
                else:
                    meth = 'feat'
            tax_qza = '%s.qza' % splitext(tax_tsv)[0]
            data.tax = [meth, tax_qza, tax_tsv]

    def set_tree_paths(self):
        for dataset, data in self.datasets.items():
            if dataset in Datasets.filt_raw:
                continue
            if data.phylo:
                odir = get_output(self.config.folder, 'phylo/%s' % dataset)
                tree_nwk = '%s/tree_%s.nwk' % (odir, dataset)
                tree_qza = '%s.qza' % splitext(tree_nwk)[0]
                if data.phylo[0] == 'amplicon':
                    intree_qza = '%s_inTree.qza' % splitext(tree_nwk)[0]
                    data.tree = (intree_qza, tree_qza, tree_nwk)
                else:
                    data.tree = ('', tree_qza, tree_nwk)

    def set_seqs_paths(self):
        for dataset, data in self.datasets.items():
            if data.phylo and data.phylo[0] == 'amplicon':
                odir = get_output(self.config.folder, 'seqs/%s' % dataset)
                seqs_fas = '%s/seq_%s.fasta' % (odir, dataset)
                seqs_qza = '%s.qza' % splitext(seqs_fas)[0]
                data.seqs = (seqs_qza, seqs_fas)

    @staticmethod
    def _get_filt_raw(dataset_):
        if dataset_ in Datasets.filt_raw:
            dataset = Datasets.filt_raw[dataset_]
        else:
            dataset = dataset_
        return dataset

    def get_precomputed_taxonomy(self, method='sklearn'):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            folder = get_output(self.config.folder, 'taxonomy/%s' % dataset)
            tax_qza = '%s/tax_%s_%s.qza' % (folder, dataset, method)
            tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            if isfile(tax_tsv) and isfile(tax_qza):
                data.tax = ['', tax_qza, tax_tsv]
            tax_qza = '%s/tax_%s.qza' % (folder, dataset)
            tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            if isfile(tax_tsv) and isfile(tax_qza):
                data.tax = ['', tax_qza, tax_tsv]

    def get_precomputed_trees(self):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            folder = get_output(self.config.folder, 'phylo/%s' % dataset)
            tree_qza = '%s/tree_%s.qza' % (folder, dataset)
            tree_nwk = '%s.nwk' % splitext(tree_qza)[0]
            if isfile(tree_nwk) and isfile(tree_qza):
                data.tree = ('', tree_qza, tree_nwk)
                data.phylo = ('precpu', 0)

    def delete_non_filtered(self):
        to_delete = [
            dat for dat in sorted(self.datasets) if dat not in self.filt_raw]
        for to_del in to_delete:
            self.datasets.pop(to_del)

    def get_taxo_levels(self):
        for dat, data in self.datasets.items():
            if not data.tax:
                continue
            if not isfile(data.tax[2]):
                print("Can't split taxonomy: %s is not present" % data.tax[2])
                continue
            if dat in Datasets.filt_raw:
                data.taxa = tuple(self.datasets[Datasets.filt_raw[dat]].taxa)
                continue
            split_taxa_fp = '%s_splitaxa.txt' % splitext(data.tax[2])[0]
            tax_pd, split_taxa_pd = get_tax_tables(data.tax[2])
            if split_taxa_pd.shape[1] > 1:
                ranks = parse_split_taxonomy(split_taxa_pd)
                split_taxa_pd = edit_split_taxonomy(ranks, split_taxa_pd)
            split_taxa_pd.to_csv(split_taxa_fp, index=True, sep='\t')
            data.taxa = (tax_pd, split_taxa_pd, split_taxa_fp)

    def set_sample_subsets(self):
        sample_subsets_issues = {}
        for dat, data in self.datasets.items():
            if data.source != dat:
                continue
            for var, vals in self.config.subsets.items():
                if var != 'ALL' and var not in set(data.metadata.columns):
                    sample_subsets_issues.setdefault(var, []).append(dat)
                    continue
                data.sample_subsets[var] = list((map(list, vals)))

        if sample_subsets_issues:
            print('Issues with sample groups set using option `-g`')
            for subset, dats in sample_subsets_issues.items():
                print(' - %s:\n\t* %s' % (subset, '\n\t* '.join(dats)))

        # give groups allocated to source dataset to all sub datasets
        for dat, data in self.datasets.items():
            if not data.sample_subsets:
                data.sample_subsets = dict(
                    self.datasets[data.source].sample_subsets)

    def meta_subset(self, meta_pd: pd.DataFrame, group: str, vals: list) -> set:
        new_meta_pd = meta_pd.copy()
        if 'ALL' in group:
            pass
        elif len([x for x in vals if x[0] == '>' or x[0] == '<']):
            for case_val in vals:
                if case_val[0] == '>':
                    new_meta_pd = new_meta_pd[
                        new_meta_pd[group].astype(float) >= float(case_val[1:])
                        ].copy()
                elif case_val[0] == '<':
                    new_meta_pd = new_meta_pd[
                        new_meta_pd[group].astype(float) <= float(case_val[1:])
                        ].copy()
        else:
            new_meta_pd = meta_pd[meta_pd[group].isin(vals)].copy()
        return set(new_meta_pd.sample_name)

    def get_meta_subsets(self):
        self.set_sample_subsets()
        for dat, data in self.datasets.items():
            for raref, tab in data.biom.items():
                subsets = {}
                if raref not in data.data:
                    data.subsets[raref] = {}
                else:
                    tab = data.data[raref]
                    sams = set(pd.Series(tab.ids(axis='sample')))
                    for group, group_vals in data.sample_subsets.items():
                        for vals in group_vals:
                            cohort = get_cohort('subsetting', group, vals, data)
                            if not cohort:
                                continue
                            s = sams & self.meta_subset(
                                data.metadata, group, vals)
                            if len(s) < 10:
                                continue
                            subsets[cohort] = (list(s), group)
                    data.subsets[raref] = subsets
