# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import os
import glob
import pandas as pd
from os.path import dirname, isdir, isfile, splitext
from skbio.tree import TreeNode
import pkg_resources
import itertools as its

from microbiome_analyzer.datasets import Datasets, Data
from microbiome_analyzer._cmds import (
    run_import, run_export, write_rarefy, write_fasta, write_collapse,
    write_sepp, write_alpha, write_filter,
    write_barplots, write_tabulate, write_alpha_correlation,
    write_alpha_rarefaction, write_volatility, write_beta, write_deicode,
    write_pcoa, write_biplot, write_emperor, write_emperor_biplot,
    write_empress, write_permanova_permdisp, write_adonis,
    write_tsne, write_umap, write_procrustes, write_mantel, write_phate,
    write_sourcetracking)
from microbiome_analyzer._filter import (
    no_filtering, get_dat_filt, filtering_thresholds,
    harsh_filtering, filter_3d)
from microbiome_analyzer._rarefy import get_digit_depth
from microbiome_analyzer._io import (
    add_q2_type, subset_meta, get_wol_tree, get_sepp_tree)
from microbiome_analyzer._taxonomy import (
    get_taxonomy_command, get_edit_taxonomy_command, get_gids,
    get_split_levels, fix_collapsed_data, find_matching_features)

RESOURCES = pkg_resources.resource_filename(
    "routine_qiime2_analyses", "resources")


class AnalysisPrep(object):
    """
    """
    analyses_commands = {}
    analyses_provenances = {}
    analyses_procrustes = {}
    analyses_mantels = {}
    analyses_nestedness = {}
    analyses_dm_decay = {}

    def __init__(self, config, project) -> None:
        self.config = config
        self.project = project
        self.analysis = {}
        self.cmds = {}
        self.messages = set()

    def get_output(self, dat, cohort=''):
        out = '%s/qiime/%s/%s' % (self.config.folder, self.analysis, dat)
        if cohort:
            out = (out + '/' + cohort).rstrip('/')
        if not isdir(out):
            os.makedirs(out)
        return out

    @staticmethod
    def make_vc(meta_pd, labs, vc_tsv):
        vcs = []
        for lab in labs:
            vc = meta_pd[lab].fillna('NaN').value_counts().reset_index()
            vc.columns = ['factor', 'number']
            vc['variable'] = lab
            vcs.append(vc)
        vc_pd = pd.concat(vcs)
        vc_pd.to_csv(vc_tsv, index=False, sep='\t')

    def import_datasets(self):
        self.analysis = 'import'
        for dat, data in self.project.datasets.items():
            qza = '%s/data/tab_%s.qza' % (self.config.folder, dat)
            # data.qza.append(qza)
            data.qza[''] = qza
            # data.biom.append('%s.biom' % splitext(qza)[0])
            data.biom[''] = '%s.biom' % splitext(qza)[0]
            cmd = run_import(
                # data.tsv[0], qza, 'FeatureTable[Frequency]')
                data.tsv[''], qza, 'FeatureTable[Frequency]')
            self.register_provenance(dat, (qza,), cmd)
            if self.config.force or not isfile(qza):
                self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def explore_filtering(self):
        def collect_filt(name, val):
            if float(val) >= 1:
                defaults[name]['count'].add(int(val))
            else:
                defaults[name]['percent'].add(float(val))
        defaults = {'preval': {'count': {0, 1, 2, 3, 5, 10, 20, 30},
                               'percent': {0, .01, .02, .05, .1}},
                    'abund': {'count': {0, 1, 2, 5, 10, 100, 1000},
                              'percent': {0, .001, .01, .02, .05, .1}}}
        defaults.update(self.config.filt3d)
        targeted = {}
        for dat, params in self.config.__dict__.get('filter', {}).items():
            if 'features' in params:
                collect_filt('abund', params['features'])
                targeted.setdefault(dat, []).append(('0', params['features']))
        self.analysis = 'explore_filters'
        for arg in self.config.__dict__.keys():
            k = arg[:-3]
            if arg.endswith('_fp') and self.config.__dict__.get(k):
                for dat in self.config.__dict__[k].get('filtering', []):
                    preval_abund = self.config.__dict__[k]['filtering'][dat]
                    for name, (preval, abund) in preval_abund.items():
                        collect_filt('preval', preval)
                        collect_filt('abund', abund)
                        targeted.setdefault(dat, []).append((preval, abund))
        self.analysis = 'filter_3D'
        for dat, data in self.project.datasets.items():
            o_dir = self.get_output(dat)
            biom = data.data['']
            for pv, ab in its.product(*[defaults['preval'], defaults['abund']]):
                filter_3d(dat, pv, ab, defaults, biom, targeted, o_dir)

    def filter(self):
        self.analysis = 'filter'
        project_filt = {}
        for dat, data in self.project.datasets.items():
            if dat not in self.config.filter:
                continue
            names = self.config.filter[dat]['names']
            thresh_sam = self.config.filter[dat]['samples']
            thresh_feat = self.config.filter[dat]['features']
            if no_filtering(dat, thresh_sam, thresh_feat):
                continue
            dat_filt = get_dat_filt(dat, names, thresh_sam, thresh_feat)
            Datasets.filt_raw[dat_filt] = dat
            Datasets.raw_filt[dat] = dat_filt
            # register the filtered dataset as an additional dataset
            data_filt = Data(dat_filt)
            data_filt.tsv = {'': data.tsv[''].replace(dat, dat_filt)}
            data_filt.biom = {'': data.biom[''].replace(dat, dat_filt)}
            data_filt.qza = {'': data.qza[''].replace(dat, dat_filt)}
            data_filt.meta = data.meta.replace(dat, dat_filt)
            data_filt.filts = dat
            data_filt.filt = dat_filt.split('%s_' % dat)[-1]
            data_filt.source = data.source
            qza_exists = isfile(data_filt.qza[''])
            meta_exists = isfile(data_filt.meta)
            if not self.config.force and qza_exists and meta_exists:
                data_filt.read_biom()
                data_filt.read_meta_pd()
            else:
                data_filt_biom, flt_cmds = filtering_thresholds(
                    names, thresh_sam, thresh_feat, data.data[''])
                if harsh_filtering(dat_filt, data_filt_biom):
                    continue
                data_filt_pd = data_filt_biom.to_dataframe(dense=True)
                data_filt_pd.index.name = '#OTU ID'
                data_filt_pd.to_csv(data_filt.tsv[''], index=True, sep='\t')
                data_filt.data.append(data_filt_biom)
                cmd = run_import(data_filt.tsv[''], data_filt.qza[''],
                                 'FeatureTable[Frequency]')
                flt_cmds += cmd
                self.register_provenance(
                    dat, (data_filt.tsv[''], data_filt.qza['']), flt_cmds)
                self.cmds.setdefault(dat_filt, []).append(cmd)
                data_filt.metadata = data.metadata
            data_filt.phylo = data.phylo
            data_filt.features = get_gids(data.features, data_filt.data[''])
            project_filt[dat_filt] = data_filt
        self.project.datasets.update(project_filt)
        self.register_command()

    def rarefy(self):
        self.analysis = 'rarefy'
        self.project.set_rarefaction()
        for dat, data in self.project.datasets.items():
            if not data.raref_depths:
                continue
            if not self.config.force:
                if self.config.filt_only and dat not in self.project.filt_raw:
                    continue
            sams_sums = data.data[''].sum()
            o_dir = self.get_output(dat)
            for ddx, depth_ in enumerate(data.raref_depths[1]):
                depth = get_digit_depth(depth_, sams_sums)
                raref = data.rarefs[ddx + 1]
                tsv = '%s/tab_%s%s.tsv' % (o_dir, dat, raref)
                biom = '%s.biom' % splitext(tsv)[0]
                qza = '%s.qza' % splitext(tsv)[0]
                cmd = write_rarefy(data.qza[''], qza, depth)
                cmd += run_export(qza, tsv, 'FeatureTable[Frequency]')
                self.register_provenance(dat, (qza,), cmd)
                if self.config.force or not isfile(tsv):
                    self.cmds.setdefault(dat, []).append(cmd)
                if isfile(tsv) and isfile(biom):
                    data.biom[raref] = biom
                    data.tsv[raref] = tsv
                    data.qza[raref] = qza
                    data.read_biom(raref)
        self.register_command()

    def taxonomy(self):
        method = 'sklearn'
        self.analysis = 'taxonomy'
        self.project.set_taxonomy_paths(method)
        for dat, data in self.project.datasets.items():
            cmd = ''
            if dat in Datasets.filt_raw:
                continue
            cmd += get_taxonomy_command(dat, self.config, data)
            if isfile(data.tax[2]):
                cmd += get_edit_taxonomy_command(data)
            if cmd:
                self.cmds.setdefault(dat, []).append(cmd)
                self.register_provenance(
                    dat, (data.tax[1], data.tax[2],), cmd)
        self.register_command()

    def shear_tree(self):
        self.analysis = 'wol'
        self.project.set_tree_paths()
        if len(Data.wols):
            wol_tree = get_wol_tree(self.config.wol_tree)
            wol = TreeNode.read(wol_tree)
            for dat, data in self.project.datasets.items():
                if dat in Datasets.filt_raw:
                    continue
                if data.phylo and data.phylo[0] == 'wol':
                    cmd = run_import(
                        data.tree[2], data.tree[1], "Phylogeny[Rooted]")
                    cmd_ = '# Shear WOL tree to genomes in data\n'
                    cmd_ += '#  - WOL tree:\t%s\n' % wol_tree
                    cmd_ += '#  - data:\t%s\n' % data.tsv['']
                    if self.config.force or not isfile(data.tree[1]):
                        wol_features = wol.shear(list(data.features.keys()))
                        nms = []
                        for tdx, tip in enumerate(wol_features.tips()):
                            if tip.name != data.features[tip.name]:
                                nms.append('#   %s\t->\t%s\n' % (
                                    tip.name, data.features[tip.name]))
                                tip.name = data.features[tip.name]
                        if nms:
                            cmd_ += '# Rename tips (from  ->  to):\n'
                            cmd_ += ''.join(nms)
                        wol_features.write(data.tree[2])
                        self.cmds.setdefault(dat, []).append(cmd)
                    cmd = cmd_ + cmd
                    self.register_provenance(
                        dat, (data.tree[1], data.tree[2],), cmd)
        self.register_command()

    def sepp(self):
        self.analysis = 'sepp'
        if len(Data.dnas):
            self.project.set_seqs_paths()
            ref_tree_qza = get_sepp_tree(self.config.sepp_tree)
            for dat, data in self.project.datasets.items():
                if data.phylo and data.phylo[0] == 'amplicon':
                    cmd_seqs = write_fasta(data.seqs[1], data.seqs[0],
                                           data.data[''], data.tsv[''])
                    in_qza = data.tree[0]
                    in_tsv = '%s.tsv' % splitext(in_qza)[0]
                    qza_out = in_qza.replace('inTree.qza', 'notInTree.qza')
                    cmd_sepp = write_sepp(data.seqs[0], ref_tree_qza,
                                          data.tree[1], data.qza[''],
                                          in_qza, in_tsv, qza_out)
                    cmd = cmd_seqs + cmd_sepp
                    self.register_provenance(dat, (data.seqs[1], data.seqs[0],
                                                   data.tree[1], in_qza), cmd)
                    if self.config.force or not isfile(data.seqs[0]):
                        self.cmds.setdefault(dat, []).append(cmd_seqs)
                    if self.config.force or not isfile(data.tree[1]):
                        self.cmds.setdefault(dat, []).append(cmd_sepp)
        self.register_command()

    @staticmethod
    def set_taxa(tax, data_tax, data, level):
        taxa_pd = data.taxa[1].iloc[:, :level]
        tax_pd = taxa_pd.apply(lambda x: '; '.join(x), axis=1).reset_index()
        tax_pd.columns = ['Feature ID', 'Taxon']
        tax_tsv = data.tax[2].replace(data.dat, '%s_tx-%s' % (data.dat, tax))
        tax_qza = '%s.qza' % splitext(tax_tsv)[0]
        tax_split = '%s_splitaxa.txt' % splitext(tax_tsv)[0]
        data_tax.tax = ('', tax_qza, tax_tsv)
        data_tax.taxa = (tax_pd, taxa_pd, tax_split)
        if not isfile(tax_tsv):
            if not isdir(dirname(tax_tsv)):
                os.makedirs(dirname(tax_tsv))
            tax_pd.to_csv(tax_tsv, index=False, sep='\t')
        cmd = ''
        if not isfile(tax_qza):
            cmd = run_import(tax_tsv, tax_qza, 'FeatureData[Taxonomy]')
        return cmd

    def collapse_taxa(self):
        self.analysis = 'collapse'
        collapse_taxo = self.config.collapse
        collapse_taxo.update(dict(
            (Datasets.raw_filt[dat], x) for dat, x in collapse_taxo.items()
            if dat in Datasets.raw_filt))
        project_coll = {}
        for dat, data in self.project.datasets.items():
            if data.source not in collapse_taxo:
                continue
            levels = collapse_taxo[data.source]
            data = self.project.datasets[dat]
            split_levels, empties = get_split_levels(levels, data.taxa[1])
            data.collapsed = split_levels
            for tax, level in split_levels.items():
                dat_tax = '%s_tx-%s' % (dat, tax)
                data_tax = Data(dat_tax)
                data_tax.source = data.source
                data_tax.rarefs = data.rarefs
                data_tax.taxon = tax
                data_tax.filt = data.filt
                data_tax.filts = data.filts
                cmd = self.set_taxa(tax, data_tax, data, level)
                for raref, tsv in data.tsv.items():
                    tax_tsv = '%s_tx-%s.tsv' % (splitext(tsv)[0], tax)
                    tax_biom = '%s.biom' % splitext(tax_tsv)[0]
                    tax_qza = '%s.qza' % splitext(tax_tsv)[0]
                    tax_meta = data.meta
                    coll_dat = splitext(tax_tsv)[0].split('/tab_')[-1]
                    if isfile(tax_tsv) and isfile(tax_biom):
                        data_tax.tsv[raref] = tax_tsv
                        data_tax.biom[raref] = tax_biom
                        data_tax.qza[raref] = tax_qza
                        data_tax.read_biom(raref)
                        if raref not in data_tax.data:
                            data_tax.tsv[raref] = ''
                            data_tax.biom[raref] = ''
                            data_tax.qza[raref] = ''
                            continue
                        cmd += fix_collapsed_data(
                            empties[tax], data_tax.data[raref], tax_tsv,
                            tax_qza, tax_meta)
                        Datasets.coll_raw[coll_dat] = dat
                        Datasets.raw_coll.setdefault(dat, []).append(coll_dat)
                        data_tax.meta = tax_meta
                        data_tax.metadata = data.metadata
                    else:
                        cmd += write_collapse(data.qza[raref], data.tax[1],
                                              tax_qza, tax_tsv, data.meta,
                                              tax_meta, level, empties[tax])
                    self.register_provenance(
                        dat, (tax_qza, tax_tsv, tax_meta,), cmd)
                if cmd:
                    self.cmds.setdefault(dat, []).append(cmd)
                data_tax.phylo = ('', 0)
                project_coll[dat_tax] = data_tax
        self.project.datasets.update(project_coll)
        self.register_command()

    def get_metrics(self, user_metrics: tuple) -> list:
        """
        Collect the alpha or beta diversity metrics from a resources file.
        """
        metrics = []
        with open('%s/%s_metrics.txt' % (RESOURCES, self.analysis)) as f:
            for line in f:
                line_strip = line.strip()
                year_version = float(self.config.qiime_env.split('-')[-1])
                if year_version > 2020.6 and line_strip == 'observed_otus':
                    line_strip = 'observed_features'
                if len(line_strip):
                    if user_metrics:
                        if line_strip in user_metrics:
                            metrics.append(line_strip)
                    else:
                        metrics.append(line_strip)
        for user_metric in user_metrics:
            if user_metric not in metrics:
                metrics.append(user_metric)
        return metrics

    def make_pies(self) -> None:
        pies_data = {}
        for dat, data in self.project.datasets.items():
            pies_data[dat] = []
            odir = self.get_output(dat)
            out_pdf = '%s/pies_%s.pdf' % (odir, dat)
            with PdfPages(out_pdf) as pdf:
                split_taxa, split_taxa_fp = split_taxa_pds[dat]
                ranks = split_taxa.columns.tolist()
                ranks_col = []
                for row in split_taxa.values:
                    ranks_col.append(ranks[(len([x for idx, x in enumerate(row)
                                                 if str(x).lstrip(
                            '%s_' % ranks[idx]).strip('_')]) - 1)])
                split_taxa['rank'] = ranks_col
                for idx, (tab_, meta) in enumerate(datasets_read[dat]):
                    nsams = tab_.shape[1]
                    cur_raref = datasets_rarefs[dat][idx]
                    tab_sum = tab_.sum(1)
                    tab_bool = tab_.astype(bool).sum(1)
                    tab = pd.concat([
                        tab_sum,
                        tab_sum / tab_sum.sum(),
                        tab_bool,
                        (tab_bool / nsams) * 100,
                        split_taxa[['rank']]
                    ], axis=1)
                    tab.columns = ['abundance', 'abundance_percent',
                                   'prevalence', 'prevalence_percent', 'rank']
                    pies_data_raref = {}
                    for min_abundance in [1, 2, 5, 10, 100]:
                        tab_min = tab.loc[tab.abundance >= min_abundance].copy()
                        abundance_bins = [int(x) for x in np.logspace(0,
                                                                      np.log10(
                                                                          tab_min[
                                                                              'abundance'].max() + 1),
                                                                      num=16)]
                        tab_min['abundance_bin'] = [
                            '%s-%s' % (abundance_bins[x], abundance_bins[x + 1])
                            for x in
                            np.digitize(tab_min['abundance'],
                                        bins=abundance_bins[1:], right=True)
                        ]

                        prevalence_bins = [1, 2, 5] + list(range(10, 101, 10))
                        tab_min['prevalence_percent_bin'] = [
                            '%s-%s' % (
                            prevalence_bins[x], prevalence_bins[x + 1]) for x in
                            np.digitize(tab_min['prevalence_percent'],
                                        bins=prevalence_bins[1:], right=True)
                        ]

                        abundances = tab_min[
                            ['rank', 'abundance_bin', 'abundance']
                        ].groupby(
                            ['rank', 'abundance_bin']).count().unstack().fillna(
                            0)
                        abundances.columns = abundances.columns.droplevel()
                        pies_data_raref['abundances'] = abundances

                        prevalences = tab_min[
                            ['rank', 'prevalence_percent_bin', 'prevalence']
                        ].groupby(['rank',
                                   'prevalence_percent_bin']).count().unstack().fillna(
                            0)
                        prevalences.columns = prevalences.columns.droplevel()
                        pies_data_raref['prevalences'] = prevalences

                        tab_gb = pd.concat([
                            tab_min.groupby('rank').count().iloc[:, 0],
                            tab_min[['rank', 'abundance']].groupby(
                                'rank').sum(),
                            tab_min[['rank', 'abundance_percent']].groupby(
                                'rank').sum(),
                        ], axis=1)
                        tab_gb.columns = ['count', 'abundance_sum',
                                          'abundance_percent_sum']
                        pies_data_raref['tab_gb'] = tab_gb

                        f = plt.figure(figsize=(6, 6))
                        plt.pie([x[0] for x in tab_gb.values],
                                labels=['%s (%s)\n%s reads' % (
                                r, row.iloc[0], row.iloc[1])
                                        for r, row in tab_gb.iterrows()],
                                autopct='%1.2f',
                                startangle=90)
                        plt.title(
                            "Number of features assigned per taxon level: %s%s\nmin %s reads, %s features" % (
                                dat, cur_raref, min_abundance,
                                tab_min.shape[0]), size=12)
                        plt.close()
                        pdf.savefig(f, bbox_inches='tight')

                        if abundances.shape[0]:
                            cols = sorted(abundances.columns.tolist(),
                                          key=lambda x: int(x.split('-')[0]))
                            ax = abundances.plot(kind='bar', stacked=True)
                            plt.ylabel('Number of features')
                            plt.title(
                                "Features per abundance group: %s%s\nmin %s reads, %s features" % (
                                    dat, cur_raref, min_abundance,
                                    tab_min.shape[0]), size=12)
                            handles, _ = ax.get_legend_handles_labels()
                            plt.legend(handles, cols)
                            pdf.savefig(bbox_inches='tight')
                            plt.close()

                        if prevalences.shape[0]:
                            cols = sorted(prevalences.columns.tolist(),
                                          key=lambda x: int(x.split('-')[0]))
                            ax = prevalences.plot(kind='bar', stacked=True)
                            plt.ylabel('Number of features')
                            plt.title(
                                "Features per prevalence group: %s%s\nmin %s reads, %s features" % (
                                    dat, cur_raref, min_abundance,
                                    tab_min.shape[0]), size=12)
                            handles, _ = ax.get_legend_handles_labels()
                            plt.legend(handles, cols)
                            pdf.savefig(bbox_inches='tight')
                            plt.close()

    def get_features_subsets(self, dat, dat_subset, data,
                             data_subset, feats, raref):
        cmd = ''
        self.analysis = 'subsets'
        o_dir = self.get_output(dat)
        tsv_subset = '%s/tab_%s%s.tsv' % (o_dir, dat_subset, raref)
        biom_subset = '%s.biom' % splitext(tsv_subset)[0]
        qza_subset = '%s.qza' % splitext(tsv_subset)[0]
        data_subset.tsv[raref] = tsv_subset
        data_subset.biom[raref] = biom_subset
        data_subset.qza[raref] = qza_subset
        data_subset.metadata = data.metadata
        data_subset.meta = data.meta
        if self.config.force or not isfile(tsv_subset):
            tmp = '%s.tmp' % splitext(tsv_subset)[0]
            subset_pd = pd.DataFrame({
                'Feature ID': feats,
                'Subset': ['tmpsubsetting'] * len(feats)})
            subset_pd.to_csv(tmp, index=False, sep='\t')
            cmd += write_filter(data.qza[raref], qza_subset, tmp)
            cmd += run_export(qza_subset, tsv_subset, 'FeatureTable')
            cmd += '\nrm %s\n\n' % tmp
            self.register_provenance(
                dat_subset, (qza_subset, tsv_subset,), cmd)
        if cmd:
            self.cmds.setdefault(dat_subset, []).append(cmd)
        if isfile(biom_subset):
            data_subset.read_biom(raref)

    def subset_features(self):
        feature_subsets = {}
        for dat, data in self.project.datasets.items():
            if dat not in self.config.feature_subsets:
                continue
            for subset, regex in self.config.feature_subsets[dat].items():
                feats = find_matching_features(data, regex)
                if not feats:
                    continue
                dat_subset = '%s_%s' % (dat, subset)
                data_subset = Data(dat_subset)
                data_subset.feat_meta = data.feat_meta
                data_subset.phylo = data.phylo
                data_subset.features = set(feats)
                data_subset.rarefs = data.rarefs
                data_subset.raref_depths = data.raref_depths
                data_subset.tax = data.tax
                data_subset.filt = data.filt
                data_subset.filts = data.filts
                data_subset.taxa = data.taxa
                data_subset.tree = data.tree
                data_subset.seqs = data.seqs
                data_subset.source = data.source
                data_subset.subset = subset
                for raref, tab in data.data.items():
                    self.get_features_subsets(
                        dat, dat_subset, data, data_subset, feats, raref)
                feature_subsets[dat_subset] = data_subset
        self.project.datasets.update(feature_subsets)
        self.register_command()

    def edits(self):
        steps = ['filt-', 'tx-', 'sub-']
        edits = pd.DataFrame([[
            d.filt, d.taxon, d.subset] for d in self.project.datasets.values()
        ], columns=steps, index=self.project.datasets.keys())
        edits = edits.loc[:, edits.astype(bool).sum() > 0].T.to_dict()
        for dat, data in self.project.datasets.items():
            data.path = data.source
            for step in steps:
                if step in edits[dat] and edits[dat][step]:
                    data.path += '/%s%s' % (step, edits[dat][step])
        self.project.edits = edits

    def barplot(self):
        self.analysis = 'barplot'
        for dat, data in self.project.datasets.items():
            if not data.taxa:
                continue
            o_dir = self.get_output(data.path)
            for raref, qza in data.qza.items():
                qzv = '%s/bar%s.qzv' % (o_dir, raref)
                if self.config.force or not isfile(qzv):
                    cmd = write_barplots(qza, qzv, data.meta, data.tax[1])
                    self.cmds.setdefault(dat, []).append(cmd)
                    self.register_provenance(dat, (qzv, qza, data.tax[1],), cmd)
        self.register_command()

    def alpha(self):
        self.analysis = 'alpha'
        metrics = self.get_metrics(self.config.alphas)
        for dat, data in self.project.datasets.items():
            for raref, qza_in in data.qza.items():
                alphas = []
                o_dir = self.get_output(data.path)
                for metric in metrics:
                    qza_out = '%s/alpha%s_%s.qza' % (o_dir, raref, metric)
                    tsv_out = '%s.tsv' % splitext(qza_out)[0]
                    cmd = write_alpha(qza_in, qza_out, data.phylo,
                                      data.tree, metric)
                    if not cmd:
                        continue
                    alphas.append([qza_out, qza_in, metric])
                    if self.config.force or not isfile(tsv_out):
                        cmd += run_export(qza_out, tsv_out, '')
                        self.register_provenance(dat, (qza_out, tsv_out,), cmd)
                        self.cmds.setdefault(dat, []).append(cmd)
                data.alpha[raref] = alphas
        self.register_command()

    def alpha_correlations(self):
        self.analysis = 'alpha_correlations'
        for dat, data in self.project.datasets.items():
            for raref, metrics_alphas in data.alpha.items():
                for method in ['spearman', 'pearson']:
                    o_dir = self.get_output('/%s/%s' % (data.path, method))
                    for (qza, _, metric) in metrics_alphas:
                        qzv = '%s/alphacorr%s_%s.qzv' % (o_dir, raref, metric)
                        if self.config.force or not isfile(qzv):
                            cmd = write_alpha_correlation(qza, qzv, method,
                                                          data.meta)
                            self.register_provenance(dat, (qzv,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def merge_alpha(self):
        self.analysis = 'tabulate'
        for dat, data in self.project.datasets.items():
            o_dir = self.get_output(data.path)
            for raref, alphas in data.alpha.items():
                qza_out = '%s/alphas%s.qzv' % (o_dir, raref)
                tsv_out = '%s.tsv' % splitext(qza_out)[0]
                data.alphas.append((tsv_out, raref))
                force = False
                if isfile(tsv_out):
                    with open(tsv_out) as f:
                        for line in f:
                            indices = line.strip().split('\t')[1:]
                            break
                    if len(indices) < len([x[-1] for x in alphas]):
                        force = True
                if self.config.force or force or not isfile(tsv_out):
                    cmd = write_tabulate(qza_out, alphas)
                    cmd += run_export(qza_out, tsv_out, '')
                    self.register_provenance(dat, (qza_out, tsv_out,), cmd)
                    if cmd:
                        self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def merge_metadata(self):
        alphas = {}
        for dat, data in self.project.datasets.items():
            for (alpha, raref) in data.alphas:
                if not isfile(alpha):
                    continue
                pre = [x + str(y) for x, y in self.project.edits[dat].items()]
                alpha_pd = pd.read_table(alpha, dtype=str)
                alpha_pd = alpha_pd.rename(
                    columns={alpha_pd.columns[0]: 'sample_name'})
                alpha_pd = alpha_pd.loc[alpha_pd.sample_name != '#q2:types']
                alpha_pd = alpha_pd.set_index('sample_name')
                alpha_pd.columns = ['alpha-%s%s' % ('__'.join(pre + [x]), raref)
                                    for x in alpha_pd.columns]
                alphas.setdefault(data.source, []).append(alpha_pd)
        for dat, alpha_pds in alphas.items():
            data = self.project.datasets[dat]
            alpha_pd = pd.concat(alpha_pds, axis=1).reset_index().rename(
                columns={'index': 'sample_name'})
            meta_pd = data.metadata[[x for x in data.metadata.columns if x not
                                     in alpha_pd.columns.tolist()[1:]]]
            merged_pd = meta_pd.merge(alpha_pd, on='sample_name', how='left')
            merged_fp = '%s_alphas.tsv' % splitext(data.meta)[0]
            merged_pd.to_csv(merged_fp, sep='\t', index=False)
            print('Written ->', merged_fp)
            for data in self.project.datasets.values():
                if data.source == dat:
                    data.metadata = merged_pd
                    data.meta = merged_fp

    def alpha_rarefaction(self):
        self.analysis = 'alpha_rarefaction'
        for dat, data in self.project.datasets.items():
            o_dir = self.get_output(data.path)
            for raref, alphas in data.alpha.items():
                for (qza, tab, m) in alphas:
                    qzv = '%s/rarcurve%s_%s.qzv' % (o_dir, raref, m)
                    if self.config.force or not isfile(qzv):
                        cmd = write_alpha_rarefaction(tab, qzv, m, data.phylo,
                                                      data.tree, data.meta,
                                                      raref)
                        self.register_provenance(dat, (qzv,), cmd)
                        self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def volatility(self):
        self.analysis = 'volatility'
        host = 'host_subject_id'  # this is subject to change
        for dat, data in self.project.datasets.items():
            if dat != data.source:
                continue
            o_dir = self.get_output(data.path)
            if self.config.longi_column not in set(data.metadata.columns):
                continue
            qzv = '%s/volatility_%s.qzv' % (o_dir, dat)
            if self.config.force or not isfile(qzv):
                cmd = write_volatility(data.meta, qzv, self.config.longi_column,
                                       host)
                self.register_provenance(dat, (qzv,), cmd)
                self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def beta(self):
        self.analysis = 'beta'
        metrics = self.get_metrics(self.config.betas)
        for dat, data in self.project.datasets.items():
            o_dir = self.get_output(data.path)
            for raref, qza in data.qza.items():
                betas = []
                for metric in metrics:
                    dm_qza = '%s/dm%s_%s.qza' % (o_dir, raref, metric)
                    dm_tsv = '%s.tsv' % splitext(dm_qza)[0]
                    if not isfile(qza):
                        continue
                    cmd = write_beta(qza, dm_qza, data.phylo, metric,
                                     data.tree, self.config)
                    if not cmd:
                        continue
                    betas.append([dm_qza, dm_tsv, metric])
                    if self.config.force or not isfile(dm_tsv):
                        cmd += run_export(dm_qza, dm_tsv, '')
                        self.register_provenance(dat, (dm_qza, dm_tsv,), cmd)
                        self.cmds.setdefault(dat, []).append(cmd)
                data.beta[raref] = betas
        self.register_command()

    def check_testing(self, data, cohort, sams) -> list:
        tests = []
        meta = data.metadata.copy()
        meta = meta.loc[meta.sample_name.isin(sams)]
        columns = set(meta.columns)
        for test in self.config.tests:
            if test not in columns:
                self.messages.add('[%s] %s has no variable %s' % (
                    self.analysis, data.source, test))
                continue
            meta_vc = meta[test].value_counts()
            if meta_vc.size > (meta.shape[0] * .8):
                self.messages.add(
                    '[%s] %s "%s" has too many (%s) factors' % (
                        self.analysis, data.source, test, meta_vc.size))
                continue
            if meta_vc.size == 1:
                if 'var-%s' % test not in cohort:
                    self.messages.add('[%s] %s "%s" has 1 factor' % (
                        self.analysis, data.source, test))
                continue
            meta_vc = meta_vc[meta_vc >= 10]
            if not meta_vc.size >= 2:
                self.messages.add(
                    '[%s] %s "%s" has <2 factors with 10 samples' % (
                        self.analysis, data.source, test))
                continue
            tests.append(test)
        return tests

    def permanova(self):
        self.analysis = 'permanova'
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                perms = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    for test in self.check_testing(data, cohort, sams):
                        cv = '%s/cv%s_%s.tsv' % (o_dir, raref, test)
                        meta = '%s/meta%s_%s.tsv' % (o_dir, raref, test)
                        meta_pd = subset_meta(data.metadata, sams, group, test)
                        if add_q2_type(meta_pd, meta, cv, [test]):
                            continue
                        for ddx, (dm, _, me) in enumerate(dms_metrics):
                            if not isfile(dm):
                                continue
                            dm_filt = '%s/dm%s_%s.qza' % (o_dir, raref, me)
                            full_cmd = ''
                            for typ in self.config.beta_type:
                                qzv = '%s/%s%s_%s__%s.qzv' % (
                                    o_dir, typ, raref, me, test)
                                html = '%s.html' % splitext(qzv)[0]
                                perms.setdefault((typ, cohort), []).append(
                                    (cv, test, html, me))
                                cmd = write_permanova_permdisp(
                                    meta, test, typ, dm, dm_filt, qzv, html)
                                self.register_provenance(dat, (qzv, dm,), cmd)
                                if self.config.force or not isfile(html):
                                    if cmd:
                                        full_cmd += cmd
                            if full_cmd:
                                full_cmd += 'rm %s\n\n' % dm_filt
                                self.cmds.setdefault(dat, []).append(full_cmd)
                data.perms[raref] = perms
        self.register_command()
        for message in sorted(self.messages):
            print(message)

    def permanova_r(self):
        self.analysis = 'permanova_r'

    def get_models_stratas(self):
        strata = ''
        yml = self.config.adonis
        if 'strata' in yml and 'global' in yml:
            strata = yml['strata']['global']
        for dat, data in self.project.datasets.items():
            source = data.source
            variables = set(data.metadata.columns)
            if 'models' in yml and source in yml['models']:
                for mod, var in yml['models'][source].items():
                    terms = set(re.split('[*/+-]', var))
                    if sorted(variables & terms) == sorted(terms):
                        stratas = list([strata])
                        if 'strata' in yml and source in yml['strata']:
                            stratas.extend(yml['strata'][source].get(mod, []))
                        stratas = set([x for x in stratas if x in variables])
                        data.adonis[mod] = [var, list(stratas)]

    def adonis(self):
        self.analysis = 'adonis'
        self.get_models_stratas()
        template = open('%s/r_scripts/adonis.R' % RESOURCES).readlines()
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics_ in data.beta.items():
                dms_metrics = [x for x in dms_metrics_ if isfile(x[0])]
                if not dms_metrics:
                    continue
                r2s = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    out = '%s/r2%s.txt' % (o_dir, raref)
                    if self.config.force or not isfile(out):
                        r_scripts = []
                        for model, (formula, stratas) in data.adonis.items():
                            variables = re.split('[*/+-]', formula)
                            terms = list(set(variables + stratas))
                            meta_pd = subset_meta(
                                data.metadata, sams, group, '', terms)
                            cv = '%s/cv_%s%s.tsv' % (o_dir, model, raref)
                            meta = '%s/meta_%s%s.tsv' % (o_dir, model, raref)
                            if add_q2_type(meta_pd, meta, cv, terms, False):
                                continue
                            r2s.setdefault(cohort, []).append((model, out))
                            r = write_adonis(meta, formula, variables, stratas,
                                             dms_metrics, out, template)
                            r_scripts.extend(r)
                        if r_scripts:
                            r_fp = '%s.R' % splitext(out)[0]
                            with open(r_fp, 'w') as o:
                                for line in r_scripts:
                                    o.write(line)
                            cmd = 'R -f %s --vanilla\n\n' % r_fp
                            self.register_provenance(dat, (out,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
                data.r2[raref] = r2s
        self.register_command()

    def deicode(self):
        self.analysis = 'deicode'
        for dat, data in self.project.datasets.items():
            for raref, qza in data.qza.items():
                rpcas = {}
                if not isfile(qza):
                    continue
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    qzv = '%s/biplot%s.qzv' % (o_dir, raref)
                    ordi = '%s/ordination%s.qza' % (o_dir, raref)
                    ordi_tsv = '%s.tsv' % splitext(ordi)[0]
                    rpcas[cohort] = ordi_tsv
                    if self.config.force or not isfile(qzv):
                        meta_pd = subset_meta(data.metadata, sams, group)
                        meta = '%s/meta%s.tsv' % (o_dir, raref)
                        meta_pd.to_csv(meta, index=False, sep='\t')
                        dm_qza = '%s/dm%s.qza' % (o_dir, raref)
                        new_qza = '%s/tab%s.qza' % (o_dir, raref)
                        cmd = write_deicode(qza, meta, new_qza, ordi,
                                            dm_qza, qzv)
                        cmd += run_export(ordi, ordi_tsv, 'pcoa')
                        self.register_provenance(dat, (ordi, qza,), cmd)
                        self.cmds.setdefault(dat, []).append(cmd)
                data.rpca[raref] = rpcas
        self.register_command()

    def tsne(self):
        self.analysis = 'tsne'
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                tsnes = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    for ddx, (dm, _, metric) in enumerate(dms_metrics):
                        if not isfile(dm):
                            continue
                        meta = '%s/meta%s' % (o_dir, raref)
                        meta_met = '%s_%s.tsv' % (meta, metric)
                        dm_filt = '%s/dm%s_%s.qza' % (o_dir, raref, metric)
                        tsne = '%s/tsne%s_%s.qza' % (o_dir, raref, metric)
                        tsne_tsv = '%s.tsv' % splitext(tsne)[0]
                        tsnes.setdefault(cohort, []).append(
                            (tsne, meta, metric))
                        if self.config.force or not isfile(tsne_tsv):
                            meta_pd = subset_meta(data.metadata, sams, group)
                            meta_pd.to_csv(meta_met, index=False, sep='\t')
                            cmd = write_tsne(dm, dm_filt, meta, meta_met,
                                             group, tsne)
                            cmd += run_export(tsne, tsne_tsv, 'pcoa')
                            self.register_provenance(dat, (tsne, dm,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
                data.tsne[raref] = tsnes
        self.register_command()

    def umap(self):
        self.analysis = 'umap'
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                umaps = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    for ddx, (dm, _, metric) in enumerate(dms_metrics):
                        if not isfile(dm):
                            continue
                        meta = '%s/meta%s' % (o_dir, raref)
                        meta_met = '%s_%s.tsv' % (meta, metric)
                        dm_filt = '%s/dm%s_%s.qza' % (o_dir, raref, metric)
                        umap = '%s/umap%s_%s.qza' % (o_dir, raref, metric)
                        umap_tsv = '%s.tsv' % splitext(umap)[0]
                        umaps.setdefault(cohort, []).append(
                            (umap, meta, metric))
                        if self.config.force or not isfile(umap_tsv):
                            meta_pd = subset_meta(data.metadata, sams, group)
                            meta_pd.to_csv(meta_met, index=False, sep='\t')
                            cmd = write_umap(dm, dm_filt, meta, meta_met,
                                             group, umap)
                            cmd += run_export(umap, umap_tsv, 'pcoa')
                            self.register_provenance(dat, (umap, dm,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
                data.umap[raref] = umaps
        self.register_command()

    def pcoa(self):
        self.analysis = 'pcoa'
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                pcoas = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    for ddx, (dm, _, metric) in enumerate(dms_metrics):
                        if not isfile(dm):
                            continue
                        meta = '%s/meta%s' % (o_dir, raref)
                        meta_met = '%s_%s.tsv' % (meta, metric)
                        dm_filt = '%s/dm%s_%s.qza' % (o_dir, raref, metric)
                        pcoa = '%s/pcoa%s_%s.qza' % (o_dir, raref, metric)
                        pcoa_tsv = '%s.tsv' % splitext(pcoa)[0]
                        pcoas.setdefault(cohort, []).append(
                            (pcoa, meta, metric))
                        if self.config.force or not isfile(pcoa_tsv):
                            meta_pd = subset_meta(data.metadata, sams, group)
                            meta_pd.to_csv(meta_met, index=False, sep='\t')
                            cmd = write_pcoa(dm, dm_filt, meta,
                                             meta_met, group, pcoa)
                            cmd += run_export(pcoa, pcoa_tsv, 'pcoa')
                            self.register_provenance(dat, (pcoa, dm,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
                data.pcoa[raref] = pcoas
        self.register_command()

    def get_pairs(self, procrustes_mantel_pairs):
        pairs = {}
        for pair, paired in procrustes_mantel_pairs.items():
            n_dats = len(paired)
            if n_dats != 2:
                print('[%s] Must be two datasets per mmvec pair (found %s)\n'
                      'Exiting\n' % (self.analysis, n_dats))
            else:
                pairs[pair] = paired
        return pairs

    def procrustes_mantel(self, analysis):
        self.analysis = analysis
        pairs = self.get_pairs(self.config.__dict__[self.analysis]['pairs'])
        for pair, (dat1, dat2) in pairs.items():
            if not {dat1, dat2}.issubset(self.project.datasets):
                continue
            data1 = self.project.datasets[dat1]
            data2 = self.project.datasets[dat2]
            for r1, r2 in its.product(*[data1.beta, data2.beta]):
                p1, p2 = data1.path, data2.path
                b1, b2 = data1.beta[r1], data2.beta[r2]
                for cohort, (sams1, group) in data1.subsets[r1].items():
                    if cohort not in data2.subsets[r2]:
                        continue
                    sams2 = data2.subsets[r2][cohort][0]
                    sams = list(set(sams1) & set(sams2))
                    if len(sams) < 10:
                        continue
                    meta = subset_meta(data1.metadata, sams, group)
                    path = p1.replace(dat1, dat1 + r1 + '__' + dat2 + r2)
                    o_dir = self.get_output((pair + '/' + path), cohort)
                    for (d1, _, m1), (d2, _, m2) in its.product(*[b1, b2]):
                        if not isfile(d1) or not isfile(d2):
                            continue
                        meta_fp = '%s/meta' % o_dir
                        meta_me = '%s_%s-%s.tsv' % (meta_fp, m1, m2)
                        d1f = '%s/dm1_%s-%s.qza' % (o_dir, m1, m2)
                        d2f = '%s/dm2_%s-%s.qza' % (o_dir, m1, m2)
                        qzv = '%s/%s_%s-%s.qzv' % (analysis, o_dir, m1, m2)
                        dis = '%s/m2_%s-%s.qza' % (o_dir, m1, m2)
                        if analysis == 'mantel':
                            out = '%s.html' % splitext(qzv)[0]
                        else:
                            out = '%s.tsv' % splitext(dis)[0]
                        AnalysisPrep.analyses_procrustes.setdefault(
                            (path, m1, m2), []).append((qzv, out))
                        if self.config.force or not isfile(out):
                            meta.to_csv(meta_me, index=False, sep='\t')
                            if analysis == 'mantel':
                                cmd = write_mantel(
                                    dat1 + r1, dat2 + r2, meta_fp, meta_me,
                                    d1, d2, d1f, d2f, qzv, out)
                            else:
                                cmd = write_procrustes(
                                    meta_fp, meta_me, d1, d2, d1f, d2f,
                                    qzv, dis, out)
                            self.register_provenance(
                                (dat1, dat2), (qzv, d1, d2), cmd)
                            self.cmds.setdefault(
                                '%s__%s' % (dat1, dat2), []).append(cmd)
        self.register_command()

    def emperor(self):
        self.analysis = 'emperor'
        for dat, data in self.project.datasets.items():
            for (typ, data_ordis) in [('pcoa', data.pcoa),
                                      ('tsne', data.tsne),
                                      ('umap', data.umap)]:
                for raref, ordis in data_ordis.items():
                    for cohort, metrics in ordis.items():
                        o_dir = self.get_output((typ + '/' + data.path), cohort)
                        for (ordi, _, metric) in metrics:
                            if not isfile(ordi):
                                continue
                            qzv = '%s/emperor%s_%s.qzv' % (o_dir, raref, metric)
                            if self.config.force or not isfile(qzv):
                                cmd = write_emperor(ordi, qzv, data.meta)
                                self.register_provenance(dat, (qzv, ordi,), cmd)
                                self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def biplots(self):
        self.analysis = 'biplot'
        for dat, data in self.project.datasets.items():
            tax = data.tax[-1]
            for raref, pcoas in data.pcoa.items():
                biplots = {}
                qza, tsv = data.qza[raref], data.tsv[raref]
                for cohort, metrics in pcoas.items():
                    o_dir = self.get_output(data.path, cohort)
                    for (pcoa, meta, metric) in metrics:
                        if not isfile(pcoa):
                            continue
                        biplot_tax = '%s/tax%s_%s.tsv' % (o_dir, raref, metric)
                        biplot = '%s/biplot%s_%s.qza' % (o_dir, raref, metric)
                        biplots.setdefault(cohort, []).append(
                            (biplot, biplot_tax, metric))
                        if self.config.force or not isfile(biplot):
                            cmd = write_biplot(
                                tsv, qza, meta, pcoa, tax, biplot, biplot_tax)
                            self.register_provenance(dat, (qza, pcoa,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
                data.biplot[raref] = biplots
        self.register_command()

    def emperor_biplot(self):
        self.analysis = 'emperor_biplot'
        for dat, data in self.project.datasets.items():
            tax = data.taxa[1]
            for raref, biplots in data.biplot.items():
                for cohort, metrics in biplots.items():
                    o_dir = self.get_output(data.path, cohort)
                    for (biplot, biplot_tax, metric) in metrics:
                        if not isfile(biplot):
                            continue
                        qzv = '%s/emperor_biplot%s_%s.qzv' % (o_dir, raref,
                                                              metric)
                        if self.config.force or not isfile(qzv):
                            cmd = write_emperor_biplot(
                                biplot, biplot_tax, data.meta, qzv, tax)
                            self.register_provenance(dat, (qzv, biplot,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def empress(self):
        self.analysis = 'empress'
        for dat, data in self.project.datasets.items():
            if not isfile(data.tree[1]):
                continue
            feat_metas = [x for x in [data.tax[-1], data.sb] if x]
            for raref, pcoas in data.pcoa.items():
                qza = data.qza[raref]
                for cohort, metrics in pcoas.items():
                    o_dir = self.get_output(data.path, cohort)
                    for (pcoa, _, metric) in metrics:
                        if not isfile(pcoa):
                            continue
                        qzv = '%s/empress%s_%s.qzv' % (o_dir, raref, metric)
                        if self.config.force or not isfile(qzv):
                            cmd = write_empress(qza, pcoa, qzv, data.meta,
                                                feat_metas, data.tree)
                            self.register_provenance(dat, (qzv, pcoa,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def empress_biplot(self):
        self.analysis = 'empress_biplot'
        for dat, data in self.project.datasets.items():
            if not isfile(data.tree[1]):
                continue
            feat_metas = [x for x in [data.tax[-1], data.sb] if x]
            for raref, biplots in data.biplot.items():
                qza = data.qza[raref]
                for cohort, metrics in biplots.items():
                    o_dir = self.get_output(data.path, cohort)
                    for (biplot, biplot_tax, metric) in metrics:
                        if not isfile(biplot):
                            continue
                        qzv = '%s/empress_biplot%s_%s.qzv' % (
                            o_dir, raref, metric)
                        if self.config.force or not isfile(qzv):
                            cmd = write_empress(
                                qza, biplot, qzv, data.meta,
                                feat_metas, data.tree)
                            self.register_provenance(dat, (qzv, biplot,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
        self.register_command()

    def make_subsets(self, cohort, txt, subsets_update):
        phate_pd = pd.read_csv(txt, sep='\t', dtype={'sample_name': str})
        phate_pd = phate_pd.loc[phate_pd['variable'].str.contains('cluster_k')]
        if len(phate_pd[['knn', 'decay', 't']].drop_duplicates()) > 2:
            print(
                'Warning: PHATE was run for >10 parameters combinations:\n'
                ' --> May be unwise (chose few, desired sets of parameters)')
            return None
        clusters = dict(phate_pd[
            ['sample_name', 'knn', 'decay', 't', 'variable', 'factor']
        ].groupby(
            ['knn', 'decay', 't', 'variable', 'factor']
        ).apply(
            func=lambda x: x.sample_name.tolist())
        )
        for (knn, decay, t, k, cluster), samples in clusters.items():
            if len(samples) >= 40:
                group = '%s/k-%s_d-%s_t-%s/k-%s/clst-%s' % (
                    cohort, knn, decay, t, k.split('cluster_k')[-1], cluster)
                subsets_update[group] = (samples, 'ALL')

    def phate(self):
        self.analysis = 'phate'
        for dat, data in self.project.datasets.items():
            labs = set()
            meta = data.metadata
            if 'labels' in self.config.phate:
                labs = set(self.config.phate['labels']) & set(meta.columns)
            for raref, qza in data.qza.items():
                phates = {}
                subsets_update = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    meta_tsv = '%s/meta%s.tsv' % (o_dir, raref)
                    new_qza = '%s/tab%s.qza' % (o_dir, raref)
                    new_tsv = '%s/tab%s.tsv' % (o_dir, raref)
                    vc_tsv = '%s/cv%s.tsv' % (o_dir, raref)
                    html = '%s/phate%s_xphate.html' % (o_dir, raref)
                    txt = '%s.tsv' % splitext(html)[0]
                    if len(glob.glob('%s/TOO_FEW.*' % o_dir)):
                        continue
                    phates.setdefault(cohort, []).append((html, txt))
                    if self.config.force or not isfile(txt):
                        meta_pd = subset_meta(meta, sams, group, '', list(labs))
                        if not isfile(vc_tsv):
                            self.make_vc(meta_pd, labs, vc_tsv)
                        meta_pd.to_csv(meta_tsv, index=False, sep='\t')
                        cmd = write_phate(self.config, qza, new_qza,
                                          new_tsv, meta_tsv, html)
                        self.cmds.setdefault(dat, []).append(cmd)
                    if isfile(txt) and self.config.phate.get('partition'):
                        self.make_subsets(cohort, txt, subsets_update)
                data.subsets[raref].update(subsets_update)
                data.phate[raref] = phates
        self.register_command()

    def sourcetracking(self):
        self.analysis = 'sourcetracking'
        print(self.config.sourcetracking)
        print(sourcetrackingdsa)
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                sourcetrackings = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    for ddx, (dm, _, metric) in enumerate(dms_metrics):
                        if not isfile(dm):
                            continue
                        meta = '%s/meta%s' % (o_dir, raref)
                        meta_met = '%s_%s.tsv' % (meta, metric)
                        dm_filt = '%s/dm%s_%s.qza' % (o_dir, raref, metric)
                        pcoa = '%s/pcoa%s_%s.qza' % (o_dir, raref, metric)
                        pcoa_tsv = '%s.tsv' % splitext(pcoa)[0]
                        sourcetrackings.setdefault(cohort, []).append(
                            (pcoa, meta, metric))
                        if self.config.force or not isfile(pcoa_tsv):
                            meta_pd = subset_meta(data.metadata, sams, group)
                            meta_pd.to_csv(meta_met, index=False, sep='\t')
                            cmd = write_sourcetracking(
                                dm, dm_filt, meta, meta_met, group, pcoa)
                            cmd += run_export(pcoa, pcoa_tsv, 'pcoa')
                            self.register_provenance(dat, (pcoa, dm,), cmd)
                            self.cmds.setdefault(dat, []).append(cmd)
                data.sourcetracking[raref] = sourcetrackings
        self.register_command()

    def register_provenance(self, dat, outputs, cmd):
        AnalysisPrep.analyses_provenances[(self.analysis, dat, outputs)] = cmd

    def register_command(self):
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.cmds = {}
