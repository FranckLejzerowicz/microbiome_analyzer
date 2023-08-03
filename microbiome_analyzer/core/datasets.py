# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
from os.path import isdir, splitext

import pandas as pd
from biom import load_table

from microbiome_analyzer.core.metadata import get_sample_subset
from microbiome_analyzer._inputs import read_meta_pd
from microbiome_analyzer._scratch import to_do, rep
from microbiome_analyzer._io_utils import (
    convert_to_biom,
    check_features,
)
from microbiome_analyzer.analyses.rarefy import get_dat_depths, get_digit_depth
from microbiome_analyzer.analyses.taxonomy import (
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
        self.adonis = {}
        self.alpha = {}
        self.alphas = []
        self.beta = {}
        self.biom = {}
        self.biplot = {}
        self.collapsed = {}
        self.data = {}
        self.doc = {}
        self.feat_meta = []
        self.features = {}
        self.filt = None
        self.filts = None
        self.meta = ''
        self.metadata = pd.DataFrame()
        self.pcoa = {}
        self.path = None
        self.perms = {}
        self.phate = {}
        self.phylo = ('', 0)
        self.qza = {}
        self.r2 = {}
        self.rarefs = ['']
        self.raref_depths = []
        self.rpca = {}
        self.sample_subsets = {}
        self.sourcetracking = {}
        self.sb = None
        self.samples = set()
        self.seqs = None
        self.source = dataset
        self.subset = ''
        self.subsets = {}
        self.tax = None
        self.taxa = None
        self.taxon = None
        self.tree = ('', '', '')
        self.tsne = {}
        self.tsv = {}
        self.umap = {}

    def read_biom(self, index=''):
        biom = load_table(rep(self.biom[index]))
        biom.remove_empty(axis='whole', inplace=True)
        self.samples = set(biom.ids(axis='sample'))
        if biom.shape[0] >= 10:
            self.data[index] = biom

    def read_tsv(self):
        data_pd = read_meta_pd(rep(self.tsv['']), 'Feature ID')
        data_pd = data_pd.set_index('Feature ID')
        empty = (data_pd.sum() == 0)
        if sum(empty):
            print('Warning: empty samples in "%s": %s' % (
                self.dat, set(empty.loc[empty].index)))
            data_pd = data_pd.loc[:, (data_pd.sum() > 0)]
            data_pd = data_pd.loc[(data_pd.sum(1) > 0), :]
        self.samples = set(data_pd.columns)
        self.data[''] = data_pd

    def read_meta_pd(self):
        meta_pd = read_meta_pd(rep(self.meta), 'sample_name')
        meta_pd = meta_pd.loc[meta_pd['sample_name'].isin(self.samples)]
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


class Datasets(object):
    """Collect the data associated with each dataset passed but the user
    """
    raw_filt = {}
    filt_raw = {}
    raw_coll = {}
    coll_raw = {}

    def __init__(self, config) -> None:
        """Initialize the class instance with the dataset name"""
        self.dir = config.dir
        self.config = config
        self.datasets = {}
        self.edits = {}
        self.key_dir = None
        self.dirs = set()

    def collect_datasets(self):
        for dat in self.config.datasets:
            data = Data(dat)
            # tsv = '%s/data/%s/data.tsv' % (self.dir, dat)
            tsv = '%s/data/%s.tsv' % (self.dir, dat)
            biom = '%s.biom' % splitext(tsv)[0]
            if to_do(biom) and to_do(tsv):
                print('[skipping] Not tsv/biom table for %s' % dat)
                continue
            # meta = '%s/metadata/%s/metadata.tsv' % (self.dir, dat)
            meta = '%s/metadata/%s.tsv' % (self.dir, dat)
            if to_do(meta):
                print('[skipping] Not metadata table for %s' % dat)
                continue
            data.meta = meta
            data.tsv[''] = tsv
            if not to_do(biom):
                data.biom[''] = biom
                data.read_biom()
            else:
                data.read_tsv()
                data.data[''] = convert_to_biom(data.data[''])
            data.read_meta_pd()
            data.check_gid_or_dna()
            data.feat_meta = self.get_feat_meta(dat)
            self.datasets[dat] = data

    def get_feat_meta(self, dat):
        paths = glob.glob('%s/feature_metadata/%s/*' % (rep(self.dir), dat))
        paths = ['${SCRATCH_FOLDER}%s' % x for x in paths]
        return paths

    def set_key_dir(self, analysis, dataset):
        """
        Get the output folder name and collect it for later folder creations.
        """
        self.key_dir = '/'.join([self.dir, analysis, dataset])
        if not isdir(rep(self.key_dir)):
            os.makedirs(rep(self.key_dir))

    def set_rarefaction_paths(self):
        for dataset, data in self.datasets.items():
            if dataset in Datasets.filt_raw:
                data.raref_depths = self.datasets[
                    Datasets.filt_raw[dataset]].raref_depths
            if not data.raref_depths:
                continue
            for depth_ in data.raref_depths[1]:
                # depth = '/raref%s' % get_digit_depth(
                depth = '_raref%s' % get_digit_depth(
                    depth_, data.data[''].sum(axis='sample'))
                data.rarefs.append(depth)

    def set_rarefaction_depths(self):
        for dataset, data in self.datasets.items():
            sam_sum = pd.Series(data.data[''].sum(axis='sample'))
            skip, depths = get_dat_depths(dataset, self.config.rarefs,
                                          rep(self.dir), sam_sum)
            if skip:
                continue
            data.raref_depths = depths

    def set_rarefaction(self):
        self.set_rarefaction_depths()
        self.set_rarefaction_paths()

    def set_taxonomy_paths(self, method):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            self.set_key_dir('taxonomy', dataset)
            # tax_tsv = '%s/taxonomy.tsv' % self.key_dir
            tax_tsv = '%s/%s.tsv' % (self.key_dir, dataset)
            if not data.phylo or data.phylo[0] != 'amplicon':
                # tax_tsv = '%s/%s_%s.tsv' % (self.key_dir, dataset, method)
            # else:
                # tax_tsv = '%s/%s.tsv' % (self.key_dir, dataset)
                if data.phylo and data.phylo[0] == 'wol':
                    method = 'wol'
                else:
                    method = 'feat'
            tax_qza = '%s.qza' % splitext(tax_tsv)[0]
            data.tax = [method, tax_qza, tax_tsv]

    def set_tree_paths(self):
        for dataset, data in self.datasets.items():
            if dataset in Datasets.filt_raw:
                continue
            if data.phylo:
                self.set_key_dir('phylogeny', dataset)
                # tree_nwk = '%s/tree.nwk' % self.key_dir
                tree_nwk = '%s/%s.nwk' % (self.key_dir, dataset)
                tree_qza = '%s.qza' % splitext(tree_nwk)[0]
                if data.phylo[0] == 'amplicon':
                    intree_qza = '%s_inTree.qza' % splitext(tree_nwk)[0]
                    data.tree = (intree_qza, tree_qza, tree_nwk)
                else:
                    data.tree = ('', tree_qza, tree_nwk)

    def set_seqs_paths(self):
        for dataset, data in self.datasets.items():
            if data.phylo and data.phylo[0] == 'amplicon':
                self.set_key_dir('sequences', dataset)
                seqs_fas = '%s/%s.fasta' % (self.key_dir, dataset)
                # seqs_fas = '%s/sequences.fasta' % self.key_dir
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
            self.set_key_dir('taxonomy', dataset)
            # tax_qza = '%s/%s_%s.qza' % (self.key_dir, dataset, method)
            # # tax_qza = '%s/taxonomy.qza' % self.key_dir
            # tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            # if not to_do(tax_tsv) and not to_do(tax_qza):
            #     data.tax = ['', tax_qza, tax_tsv]
            tax_qza = '%s/%s.qza' % (self.key_dir, dataset)
            tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            if not to_do(tax_tsv) and not to_do(tax_qza):
                data.tax = ['', tax_qza, tax_tsv]

    def get_precomputed_trees(self):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            self.set_key_dir('phylogeny', dataset)
            tree_qza = '%s/%s.qza' % (self.key_dir, dataset)
            # tree_qza = '%s/tree.qza' % self.key_dir
            tree_nwk = '%s.nwk' % splitext(tree_qza)[0]
            if not to_do(tree_nwk):
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
            tax_tsv = data.tax[2]
            if to_do(tax_tsv):
                print("Can't split taxonomy: %s is not present" % rep(tax_tsv))
                continue
            if dat in Datasets.filt_raw:
                data.taxa = tuple(self.datasets[Datasets.filt_raw[dat]].taxa)
                continue
            split_taxa_fp = '%s_splitaxa.txt' % splitext(tax_tsv)[0]
            tax_pd, split_taxa_pd = get_tax_tables(rep(tax_tsv))
            if split_taxa_pd.shape[1] > 1:
                ranks = parse_split_taxonomy(split_taxa_pd)
                split_taxa_pd = edit_split_taxonomy(ranks, split_taxa_pd)
            split_taxa_pd.to_csv(rep(split_taxa_fp), index=True, sep='\t')
            data.taxa = (tax_pd, split_taxa_pd, split_taxa_fp)

    def set_sample_subsets(self):
        sample_subsets_issues = {}
        for dat, data in self.datasets.items():
            if data.source != dat:
                continue
            meta_cols = set(data.metadata.columns)
            for name, vars_vals in self.config.subsets.items():
                data.sample_subsets[name] = {}
                for var, vals in vars_vals.items():
                    if var not in meta_cols:
                        sample_subsets_issues.setdefault(var, []).append(dat)
                        continue
                    data.sample_subsets[name][var] = vals

        if sample_subsets_issues:
            print('[sample_subsets] Issues with sample groups set using option')
            for subset, dats in sample_subsets_issues.items():
                print(' - %s:\n\t* %s' % (subset, '\n\t* '.join(dats)))

        # give groups allocated to source dataset to all sub datasets
        for dat, data in self.datasets.items():
            if not data.sample_subsets:
                data.sample_subsets = dict(
                    self.datasets[data.source].sample_subsets)

    def get_meta_subsets(self):
        self.set_sample_subsets()
        for dat, data in self.datasets.items():
            for raref, tab in data.biom.items():
                subsets = {}
                tab = data.data[raref]
                tab_sams = set(pd.Series(tab.ids(axis='sample')))
                for name, vars_vals in data.sample_subsets.items():
                    if name == 'ALL':
                        sams = tab_sams
                    else:
                        vars_sams = get_sample_subset(
                            data.metadata, dat, vars_vals)
                        sams = tab_sams & set.intersection(*vars_sams)
                    if len(sams) < 4:
                        continue
                    subsets[name] = (list(sams), list(vars_vals))
                data.subsets[raref] = subsets
