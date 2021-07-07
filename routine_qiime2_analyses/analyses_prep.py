# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from os.path import isfile, splitext
from skbio.tree import TreeNode
from routine_qiime2_analyses.dataset_collection import Datasets, Data
from routine_qiime2_analyses._routine_q2_cmds import (
    run_import, run_export, write_rarefy, write_seqs_fasta,
    write_fragment_insertion
)
from routine_qiime2_analyses._routine_q2_filter import (
    get_thresholds, no_filtering, get_dat_filt, filtering_thresholds,
    harsh_filtering
)
from routine_qiime2_analyses._routine_q2_rarefy import (
    get_digit_depth,
)
from routine_qiime2_analyses._routine_q2_io_utils import (
    read_yaml_file
)
from routine_qiime2_analyses._routine_q2_taxonomy import (
    get_taxonomy_command, get_edit_taxonomy_command, get_gids,
    get_split_levels, fix_collapsed_data, write_collapse_taxo
)
from routine_qiime2_analyses._routine_q2_phylo import (
    get_wol_tree, get_sepp_tree
)


class AnalysisPrep(object):
    """
    """
    analyses_commands = {}

    # def __init__(self, config, analysis) -> None:
    def __init__(self, config, project) -> None:
        self.config = config
        self.project = project
        self.cmds = {}

    # def import_datasets(self, config, project):
    def import_datasets(self):
        for dataset, data in self.project.datasets.items():
            if self.config.force or not isfile(data.qza[0]):
                cmd = run_import(
                    data.tsv[0], data.qza[0], 'FeatureTable[Frequency]')
                self.cmds.setdefault(dataset, []).append(cmd)
        self.register_command('import')

    # def filter_rare_samples(self, config, project):
    def filter_rare_samples(self):
        thresholds = read_yaml_file(self.config.filt_threshs)
        project_filt = {}
        for dat, data in self.project.datasets.items():
            if dat not in thresholds:
                continue
            names, thresh_sam, thresh_feat = get_thresholds(thresholds[dat])
            if no_filtering(dat, thresh_sam, thresh_feat):
                continue
            dat_filt = get_dat_filt(dat, names, thresh_sam, thresh_feat)
            Datasets.filt_raw[dat_filt] = dat
            Datasets.raw_filt[dat] = dat_filt
            # register the filtered dataset as an additional dataset
            data_filt = Data(dat_filt, self.config.i_datasets_folder)
            qza_exists = isfile(data_filt.qza[0])
            meta_exists = isfile(data_filt.meta[0])
            if not self.config.force and qza_exists and meta_exists:
                data_filt.read_data_pd()
                data_filt.read_meta_pd()
            else:
                data_filt_pd = filtering_thresholds(names, thresh_sam,
                                                    thresh_feat, data.data[0])
                if harsh_filtering(dat_filt, data_filt_pd):
                    continue
                # write filtered data
                data_filt_pd.to_csv(data_filt.tsv[0], index=True, sep='\t')
                data_filt.data.append(data_filt_pd)
                # import qza
                cmd = run_import(data_filt.tsv[0], data_filt.qza[0],
                                 'FeatureTable[Frequency]')
                self.cmds.setdefault(dat_filt, []).append(cmd)
                # write filtered metadata
                meta_filt_pd = data.metadata[0].loc[
                    data.metadata[0].sample_name.isin(
                        data_filt_pd.columns.tolist())].copy()
                meta_filt_pd.to_csv(data_filt.meta[0], index=False, sep='\t')
                data_filt.metadata.append(meta_filt_pd)
            data_filt.phylo = data.phylo
            data_filt.features = get_gids(data.features, data_filt.data[0])
            project_filt[dat_filt] = data_filt
        self.project.datasets.update(project_filt)
        self.register_command('filter')

    def rarefy(self):
        self.project.set_rarefaction(self.config)
        for dat, data in self.project.datasets.items():
            if not self.config.force:
                if self.config.filt_only and dat not in self.project.filt_raw:
                    continue
            sams_sums = data.data[0].sum()
            for dx, depth_ in enumerate(data.raref_depths[1]):
                depth = get_digit_depth(depth_, sams_sums)
                remaining_samples = sams_sums[sams_sums >= depth].index
                meta_raref_pd = data.metadata[0].loc[
                    data.metadata[0].sample_name.isin(remaining_samples), :]
                meta_raref_pd.to_csv(data.meta[dx+1], index=False, sep='\t')
                data.metadata.append(meta_raref_pd)
                if self.config.force or not isfile(data.tsv[dx+1]):
                    cmd = write_rarefy(data.qza[0], data.qza[dx+1], depth)
                    cmd += run_export(data.qza[dx+1], data.tsv[dx+1],
                                      'FeatureTable[Frequency]')
                    self.cmds.setdefault(dat, []).append(cmd)
                if isfile(data.tsv[dx+1]):
                    tab_filt_pd = pd.read_csv(
                        data.tsv[dx+1], index_col=0, header=0, sep='\t')
                    data.data.append(tab_filt_pd)
                else:
                    data.data.append('raref')
                    data.metadata.append(str(depth))
        self.register_command('rarefy')

    def taxonomy(self):
        method = 'sklearn'
        self.project.set_taxonomy_paths(self.config, method)
        for dat, data in self.project.datasets.items():
            if dat in Datasets.filt_raw:
                continue
            cmd = get_taxonomy_command(dat, self.config, data)
            if cmd:
                self.cmds.setdefault(dat, []).append(cmd)
            if isfile(data.tax[2]):
                cmd = get_edit_taxonomy_command(data)
                if cmd:
                    self.cmds.setdefault(dat, []).append(cmd)
        self.register_command('taxonomy')

    def shear_tree(self):
        self.project.set_tree_paths(self.config)
        if len(Data.wols):
            i_wol_tree = get_wol_tree(self.config.i_wol_tree)
            wol = TreeNode.read(i_wol_tree)
            for dat, data in self.project.datasets.items():
                if dat in Datasets.filt_raw:
                    continue
                if data.phylo and data.phylo[0] == 'wol':
                    if self.config.force or not isfile(data.tree[1]):
                        wol_features = wol.shear(list(data.features.keys()))
                        for tip in wol_features.tips():
                            tip.name = data.features[tip.name]
                        wol_features.write(data.tree[2])
                        cmd = run_import(
                            data.tree[2], data.tree[1], "Phylogeny[Rooted]")
                        self.cmds.setdefault(dat, []).append(cmd)
        self.register_command('wol')

    def sepp(self):
        if len(Data.dnas):
            self.project.set_seqs_paths(self.config)
            ref_tree_qza = get_sepp_tree(self.config.i_sepp_tree)
            for dat, data in self.project.datasets.items():
                if data.phylo and data.phylo[0] == 'amplicon':
                    if self.config.force or not isfile(data.seqs[0]):
                        cmd = write_seqs_fasta(
                            data.seqs[1], data.seqs[0], data.data[0])
                        self.cmds.setdefault(dat, []).append(cmd)
                    if self.config.force or not isfile(data.tree[1]):
                        cmd = write_fragment_insertion(
                            data.seqs[0], ref_tree_qza, data.tree[1],
                            data.qza[0], data.tree[0])
                        self.cmds.setdefault(dat, []).append(cmd)
        self.register_command('sepp')

    def collapse(self):
        collapse_taxo = read_yaml_file(self.config.collapse_taxo)
        collapse_taxo = dict(
            (Datasets.raw_filt[dat], x) for dat, x in collapse_taxo.items()
            if dat in Datasets.raw_filt and dat in self.project.datasets)
        project_coll = {}
        for dat, levels in collapse_taxo.items():
            data = self.project.datasets[dat]
            split_levels, empties = get_split_levels(levels, data.tax_split[0])
            data.collapsed = split_levels
            for tax, level in split_levels.items():
                dat_tax = '%s_tx-%s' % (dat, tax)
                data_tax = Data(dat_tax, self.config.i_datasets_folder)
                for idx, tsv in enumerate(data.tsv):
                    tax_tsv = '%s_tx-%s.tsv' % (splitext(tsv)[0], tax)
                    tax_qza = '%s.qza' % splitext(tax_tsv)[0]
                    tax_meta = '%s_tx-%s.tsv' % (
                        splitext(data.meta[idx])[0], tax)
                    coll_dat = splitext(tax_tsv)[0].split('/tab_')[-1]
                    if isfile(tax_tsv) and isfile(tax_meta):
                        coll_pd = pd.read_csv(tax_tsv, index_col=0,
                                              header=0, sep='\t')
                        if coll_pd.shape[0] < 5:
                            continue
                        cmd = fix_collapsed_data(
                            empties[tax], coll_pd, tax_tsv, tax_qza, tax_meta)
                        if cmd:
                            self.cmds.setdefault(dat, []).append(cmd)
                        Datasets.coll_raw[coll_dat] = dat
                        Datasets.raw_coll.setdefault(dat, []).append(coll_dat)
                        if idx:
                            data_tax.tsv.append(tax_tsv)
                            data_tax.qza.append(tax_qza)
                            data_tax.meta.append(tax_meta)
                        data_tax.data.append(coll_pd)
                        data_tax.metadata.append(coll_pd)
                        data_tax.rarefs.append(data.rarefs[idx])
                    else:
                        cmd = write_collapse_taxo(
                            data.qza[idx], data.tax[1], tax_qza, tax_tsv,
                            data.meta[idx], tax_meta, level, empties[tax])
                        if cmd:
                            self.cmds.setdefault(dat, []).append(cmd)
                data_tax.phylo = ('', 0)
                project_coll[dat_tax] = data_tax
        self.project.datasets.update(project_coll)
        self.register_command('collapse')

    def register_command(self, analysis):
        AnalysisPrep.analyses_commands[analysis] = self.cmds
        self.cmds = {}

