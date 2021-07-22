# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import isfile, splitext
from routine_qiime2_analyses._routine_q2_io_utils import (
    check_input,
    read_yaml_file,
    read_meta_pd,
    check_features,
    get_analysis_folder
)
from routine_qiime2_analyses._routine_q2_rarefy import (
    get_dat_depths, get_digit_depth
)
from routine_qiime2_analyses._routine_q2_taxonomy import (
    get_tax_tables, parse_split_taxonomy, edit_split_taxonomy
)


class Data(object):
    """
    """
    wols = []
    dnas = []
    names = []

    def __init__(self, dataset, datasets_folder):
        self.dataset = dataset
        Data.names.append(dataset)
        self.tsv = ['%s/data/tab_%s.tsv' % (datasets_folder, dataset)]
        self.qza = ['%s/data/tab_%s.qza' % (datasets_folder, dataset)]
        self.meta = ['%s/metadata/meta_%s.tsv' % (datasets_folder, dataset)]
        self.data = []
        self.metadata = []
        self.phylo = None
        self.features = {}
        self.rarefs = ['']
        self.raref_depths = []
        self.tax = None
        self.taxonomy = None
        self.tax_split = None
        self.tree = None
        self.seqs = None
        self.alpha = []
        self.beta = []

    def read_data_pd(self):
        data_pd = read_meta_pd(self.tsv[0], '#OTU ID').set_index('#OTU ID')
        self.data.append(data_pd)

    def read_meta_pd(self):
        meta_pd = read_meta_pd(self.meta[0], 'sample_name')
        self.metadata.append(meta_pd)

    def check_gid_or_dna(self):
        gids, dna, correct = check_features(self.data[0])
        if len(gids) == self.data[0].shape[0]:
            Data.wols.append(self.dataset)
            self.features = gids
            if correct:
                new_index = self.data[0].index.str.replace(r'[; ]+', '|')
                self.data[0].index = new_index
                self.phylo = ('wol', 1)
            else:
                self.phylo = ('wol', 0)
        if dna:
            Data.dnas.append(self.dataset)
            self.phylo = ('amplicon', 0)


class Datasets(object):
    """Collect the data associated with each dataset passed but the user
    """
    raw_filt = {}
    filt_raw = {}
    raw_coll = {}
    coll_raw = {}

    def __init__(self, i_datasets, i_datasets_folder) -> None:
        """Initialize the class instance with the dataset name"""
        self.passed_datasets = i_datasets
        self.datasets_folder = check_input(i_datasets_folder)
        self.datasets = {}
        self.collect_datasets()
        self.read_datasets()

    def collect_datasets(self):
        for dataset in self.passed_datasets:
            data = Data(dataset, self.datasets_folder)
            if not isfile(data.tsv[0]):
                print(data.tsv[0], 'does not exist\n Skipping', data)
                continue
            if not isfile(data.meta[0]):
                print(data.meta, 'does not exist\n Skipping', data)
                continue
            self.datasets[dataset] = data

    def read_datasets(self):
        for dataset, data in self.datasets.items():
            data.read_data_pd()
            data.read_meta_pd()
            data.check_gid_or_dna()

    def set_rarefaction_paths(self, config):
        for dataset, data in self.datasets.items():
            if dataset in Datasets.filt_raw:
                data.raref_depths = self.datasets[
                    Datasets.filt_raw[dataset]
                ].raref_depths
            if not data.raref_depths:
                continue
            odir = get_analysis_folder(
                config.i_datasets_folder, 'rarefy/%s' % dataset)
            for depth_ in data.raref_depths[1]:
                depth = '_raref%s' % get_digit_depth(depth_, data.data[0].sum())
                data.tsv.append('%s/tab_%s%s.tsv' % (odir, dataset, depth))
                data.qza.append('%s/tab_%s%s.qza' % (odir, dataset, depth))
                data.meta.append('%s/meta_%s%s.tsv' % (odir, dataset, depth))
                data.rarefs.append(depth)

    def set_rarefaction_depths(self, config, depths_yml):
        for dataset, data in self.datasets.items():
            sam_sum = data.data[0].sum()
            skip, depths = get_dat_depths(
                dataset, config.i_datasets_folder, depths_yml, sam_sum)
            if skip:
                continue
            data.raref_depths = depths

    def set_rarefaction(self, config):
        depths_yml = read_yaml_file(config.raref_depths)
        self.set_rarefaction_depths(config, depths_yml)
        self.set_rarefaction_paths(config)

    def set_taxonomy_paths(self, config, method):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            odir = get_analysis_folder(config.i_datasets_folder,
                                       'taxonomy/%s' % dataset)
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

    def set_tree_paths(self, config):
        for dataset, data in self.datasets.items():
            if dataset in Datasets.filt_raw:
                continue
            if data.phylo:
                odir = get_analysis_folder(config.i_datasets_folder,
                                           'phylo/%s' % dataset)
                tree_nwk = '%s/tree_%s.nwk' % (odir, dataset)
                tree_qza = '%s.qza' % splitext(tree_nwk)[0]
                if data.phylo[0] == 'amplicon':
                    intree_qza = '%s_inTree.qza' % splitext(tree_nwk)[0]
                    data.tree = (intree_qza, tree_qza, tree_nwk)
                else:
                    data.tree = ('', tree_qza, tree_nwk)

    def set_seqs_paths(self, config):
        for dataset, data in self.datasets.items():
            if data.phylo and data.phylo[0] == 'amplicon':
                odir = get_analysis_folder(config.i_datasets_folder,
                                           'seqs/%s' % dataset)
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

    def get_precomputed_taxonomy(self, config, method='sklearn'):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            analysis_folder = get_analysis_folder(
                config.i_datasets_folder, 'taxonomy/%s' % dataset)
            tax_qza = '%s/tax_%s_%s.qza' % (analysis_folder, dataset, method)
            tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            if isfile(tax_tsv) and isfile(tax_qza):
                data.tax = ['', tax_qza, tax_tsv]

            tax_qza = '%s/tax_%s.qza' % (analysis_folder, dataset)
            tax_tsv = '%s.tsv' % splitext(tax_qza)[0]
            if isfile(tax_tsv) and isfile(tax_qza):
                data.tax = ['', tax_qza, tax_tsv]

    def get_precomputed_trees(self, config):
        for dataset_, data in self.datasets.items():
            dataset = self._get_filt_raw(dataset_)
            analysis_folder = get_analysis_folder(
                config.i_datasets_folder, 'phylo/%s' % dataset)
            tree_qza = '%s/tree_%s.qza' % (analysis_folder, dataset)
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
                data.tax_split = tuple(self.datasets[
                    Datasets.filt_raw[dat]].tax_split)
                continue
            split_taxa_fp = '%s_splitaxa.txt' % splitext(data.tax[2])[0]
            tax_pd, split_taxa_pd = get_tax_tables(data.tax[2])
            if split_taxa_pd.shape[1] > 1:
                ranks = parse_split_taxonomy(split_taxa_pd)
                split_taxa_pd = edit_split_taxonomy(ranks, split_taxa_pd)
            split_taxa_pd.to_csv(split_taxa_fp, index=True, sep='\t')
            data.tax_split = (split_taxa_pd, split_taxa_fp)
