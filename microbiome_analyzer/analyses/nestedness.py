# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import numpy as np
import pandas as pd
import pkg_resources

from os.path import basename, isdir, isfile, splitext
from microbiome_analyzer.core.datasets import Datasets
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import (
    write_filter, run_export, run_add_metadata)
from microbiome_analyzer._io_utils import get_cohort, subset_meta
from microbiome_analyzer._scratch import io_update, to_do, rep


RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources/python_scripts")


class Nestedness(object):

    def __init__(self, config, project):
        self.config = config
        self.project = project
        self.dir = project.dir
        self.dirs = project.dirs
        self.out = ''
        self.analysis = 'nestedness'
        self.ios = {}
        self.cmds = {}
        self.subsets = {'ALL': [[]]}
        self.nodfs = []
        self.nodfs_vars = []
        self.colors = {'sample': [], 'feature': []}
        self.nulls = ['equiprobablefixed']
        self.modes = ['betweeneachpairoftypes']
        self.params = {'iterations': 100}
        self.binary = None
        self.res = []
        self.res_pd = pd.DataFrame
        self.graphs = []
        # run things
        self.run()
        self.figure_nodfs()
        self.figure_graphs()
        AnalysisPrep.analyses_nestedness = self.res_pd

    def get_output(self, dat, cohort):
        self.out = '%s/%s/%s' % (self.dir, self.analysis, dat)
        if cohort:
            self.out = (self.out + '/' + cohort).rstrip('/')
        self.out = self.out.replace(' ', '_')
        if not isdir(rep(self.out)):
            os.makedirs(rep(self.out))

    def get_path(self, path_):
        path = path_
        params = self.config.run_params[self.analysis]
        if not self.config.jobs or not params['scratch']:
            path = rep(path_)
        return path

    def get_config(self):
        conf = self.read_config()
        if conf:
            if 'subsets' in conf:
                self.subsets.update(conf['subsets'])
            if 'nodfs' in conf:
                self.nodfs.extend(conf['nodfs'])
            if 'sample_colors' in conf:
                self.colors['sample'].extend(conf['sample_colors'])
            if 'feature_colors' in conf:
                self.colors['feature'].extend(conf['feature_colors'])
            if 'nulls' in conf:
                self.nulls = conf['nulls']
            if 'modes' in conf:
                self.modes = conf['modes']
                if 'overall' in self.modes:
                    self.modes = [
                        x for x in self.modes if x != 'overall'] + ['overall']
            if 'params' in conf:
                self.params.update(conf['params'])

    def read_config(self):
        conf = self.config.nestedness
        if 'soft' not in conf:
            print('Need path to the Nestedness soft')
            return None
        if conf['soft'].endswith('Autocorrelation.jar') and isfile(
                conf['soft']):
            self.binary = conf['soft']
        else:
            self.binary = '%s/bin/Autocorrelation.jar' % conf['soft']
            if not isfile(self.binary):
                print('Need path to the Nestedness soft')
                return None
        return conf

    def get_taxo_level(self, dat, data):
        level = 'feature'
        if dat in Datasets.coll_raw:
            level = data.taxon
        return level

    def write_meta(self, meta, full_meta_pd, meta_pd, nodfs):
        meta_pd = full_meta_pd.loc[full_meta_pd.sample_name.isin(
            meta_pd.sample_name)]
        lat_lon_date = ['latitude', 'longitude', 'datetime']
        cols = set()
        self.nodfs_vars = []
        for col in (nodfs + lat_lon_date):
            if col not in set(meta_pd.columns):
                continue
            col_series = meta_pd[col]
            if col_series.unique().size == 1:
                continue
            if col not in lat_lon_date and min(col_series.value_counts()) == 1:
                continue
            cols.add(col)
            if col in self.nodfs:
                self.nodfs_vars.append(col)
        meta_pd = meta_pd[(['sample_name'] + sorted(cols))]
        meta_pd.columns = ['#SampleID'] + sorted(cols)
        meta_pd = meta_pd.loc[~meta_pd[self.nodfs_vars].isna().any(axis=1)]
        meta_pd.to_csv(rep(meta), index=False, sep='\t')

    def write_fields(self, o_dir, raref):
        fields = '%s/fields%s.txt' % (o_dir, raref)
        with open(fields, 'w') as fields_o:
            for ndx, nodf in enumerate(self.nodfs_vars):
                fields_o.write('%s\n' % nodf)
        return fields

    def write_nestedness_graph(self, biom: str, graph: str) -> str:
        """
        https://github.com/jladau/Nestedness
        """
        cmd = ''
        if to_do(graph):
            cmd += 'java -cp %s' % self.binary
            cmd += 'edu.ucsf.Nestedness.Grapher.GrapherLauncher'
            cmd += ' --sBIOMPath=%s' % biom
            cmd += ' --bCheckRarefied=false'
            cmd += ' --bNormalize=true'
            cmd += ' --bPresenceAbsence=false'
            cmd += ' --sTaxonRank=otu'
            cmd += ' --sOutputPath=%s' % graph
            cmd += ' --rgsSampleMetadataFields=%s\n' % ','.join(self.nodfs_vars)
        return cmd

    def write_comparisons(self, biom: str, comp: str,
                          mode: str, nodf: str) -> str:
        """
        https://github.com/jladau/Nestedness
        """
        cmd = ''
        if self.config.force or to_do(comp):
            cmd += 'java -Xmx5g -cp %s' % self.binary
            cmd += 'edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher'
            cmd += ' --sBIOMPath=%s' % biom
            cmd += ' --sOutputPath=%s' % comp
            cmd += ' --bCheckRarefied=false'
            cmd += ' --bNormalize=true'
            cmd += ' --bPresenceAbsence=false'
            cmd += ' --sTaxonRank=otu'
            if mode in ["betweeneachpairoftypes", "withineachtype"]:
                cmd += ' --sMetadataField=%s' % nodf
            cmd += ' --iRandomSeed=1234'
            cmd += ' --sComparisonMode=%s' % mode
            cmd += ' --iNestednessPairs=1000'
            cmd += ' --sNestednessAxis=sample\n'
        return cmd

    def write_simulations(self, biom: str, simul: str,
                          comp: str, null: str) -> str:
        """
        https://github.com/jladau/Nestedness
        """
        cmd = ''
        if self.config.force or to_do(simul):
            cmd += 'java -cp %s' % self.binary
            cmd += 'edu.ucsf.Nestedness.Calculator.CalculatorLauncher'
            cmd += ' --sBIOMPath=%s' % biom
            cmd += ' --sOutputPath=%s' % simul
            cmd += ' --bCheckRarefied=false'
            cmd += ' --bNormalize=true'
            cmd += ' --bPresenceAbsence=false'
            cmd += ' --sTaxonRank=otu'
            cmd += ' --sComparisonsPath=%s' % comp
            cmd += ' --iNullModelIterations=%s' % str(
                self.params['iterations'])
            cmd += ' --bOrderedNODF=false'
            cmd += ' --sNestednessAxis=sample'
            cmd += ' --sNestednessNullModel=%s' % null
            cmd += ' --bSimulate=true\n'
        return cmd

    def merge_outputs(self, mode_dir: str, level: str,
                      raref: str, cohort: str, nodfs: str):
        com_sta_pds = []
        mode = mode_dir.split('/')[-1].replace('output_', '')
        for sim_fp in glob.glob('%s/simulations_*%s.csv' % (
                rep(mode_dir), raref)):
            sim_pd = pd.read_csv(sim_fp)
            null = splitext(basename(sim_fp))[0].split('_')[-1]
            sim_pd['NULL'] = null
            if mode == 'overall':
                meta_name = np.nan
            else:
                meta_name = '_'.join(basename(sim_fp).split('_')[1:-1])
            sim_pd['METADATA'] = meta_name
            sim_pd['MODE'] = mode
            sim_pd['LEVEL'] = level
            sim_pd['RAREF'] = raref.replace('_raref', '')
            sim_pd['SAMPLE_SUBSET'] = cohort
            bins_lt = ['', '>', '>>', '>>>']
            bins_gt = ['', '<', '<<', '<<<']
            bins_vals = [0, 0.95, 0.97, 0.99]
            ps = []
            for (p, b) in [('PR_LT_OBSERVED', bins_lt),
                           ('PR_GT_OBSERVED', bins_gt)]:
                ps.append(p)
                sim_pd[p] = [
                    b[x - 1] for x in np.digitize(sim_pd[p], bins=bins_vals)]
            sim_pd['PVALUE'] = sim_pd[ps].apply(
                func=lambda x: ''.join(x), axis=1)
            sim_pd = sim_pd.drop(columns=(['PR_ET_OBSERVED', 'NODF_SES'] + ps))
            comp_fp = sim_fp.replace('simulations', 'comparisons')
            comp_pd = pd.read_csv(comp_fp.replace('_%s' % null, ''))
            comp_pd['COMPARISON'] = comp_pd.fillna('overall')[
                ['VERTEX_1_CLASSIFICATION', 'VERTEX_2_CLASSIFICATION']
            ].apply(func=lambda x: '-'.join(list(set(x))), axis=1)
            comp_pd = comp_pd[ ['GRAPH_ID', 'COMPARISON']].drop_duplicates()
            nodf_pd = sim_pd[sim_pd['GRAPH_EDGE_COUNT'] > 5].merge(
                comp_pd, on=['GRAPH_ID'], how='left')
            com_sta_pds.append(nodf_pd)
        if com_sta_pds:
            com_sta_pd = pd.concat(com_sta_pds)
            com_sta_pd.to_csv(rep(nodfs), index=False, sep='\t')

    def init_cmds(self, out_dir, meta: str, data, raref: str):
        qza = '%s/tab%s.qza' % (out_dir, raref)
        biom_ = '%s.biom' % splitext(qza)[0]
        tsv = '%s.tsv' % splitext(qza)[0]
        biom = '%s_w-md.biom' % splitext(qza)[0]
        cmd = ''
        if to_do(biom_):
            cmd += write_filter(data.qza[raref], qza, meta)
            cmd += run_export(qza, tsv, 'FeatureTable')
            cmd += 'rm %s %s\n' % (qza, tsv)
            i_f = []
            if isfile(meta):
                io_update(self, i_f=[data.qza[raref], meta], o_f=biom_,
                          key=data.dat)
        else:
            io_update(self, i_f=biom_, key=data.dat)

        if to_do(biom):
            cmd += run_add_metadata(biom_, biom, meta)
            cmd += 'rm %s\n' % biom_
            if to_do(biom_):
                io_update(self, o_f=biom, key=data.dat)
            else:
                io_update(self, i_f=biom_, o_f=biom, key=data.dat)
        else:
            io_update(self, i_f=biom, key=data.dat)

        return cmd, biom

    def write_graph(self, data, meta_pd, raref) -> (str, str):
        graph = '%s/graphs%s.csv' % (self.out, raref)
        cmd = ''
        if self.config.force or to_do(graph):
            meta = '%s/meta.tsv' % self.out
            self.write_meta(meta, data.metadata, meta_pd, self.nodfs)
            cmd, biom = self.init_cmds(self.out, meta, data, raref)
            cmd += self.write_nestedness_graph(biom, graph)
        return cmd, graph

    def scripts(
            self,
            meta_pd: pd.DataFrame,
            data,
            raref: str,
            level: str,
            cohort: str
    ):
        """

        Parameters
        ----------
        meta_pd
        data
        raref
        level
        cohort

        Returns
        -------

        """
        dat = Datasets.coll_raw.get(data.dat, data.dat)
        cmd, graph = self.write_graph(data, meta_pd, raref)
        self.graphs.append([dat, data.dat, cohort, raref, graph])
        for mdx, mode in enumerate(self.modes):
            print()
            print()
            print()
            print('mode: [',  mdx, ']', mode)
            m_dir = '%s/output_%s%s' % (self.out, mode, raref.replace('_', '/'))
            if not isdir(rep(m_dir)):
                os.makedirs(rep(m_dir))
            nodfs_fp = '%s/nodfs.tsv' % m_dir
            if mode == 'overall':
                nodfs = ['sample_name']
            else:
                nodfs = self.nodfs
            meta = '%s/meta.tsv' % m_dir
            init_cmd, biom = self.init_cmds(m_dir, meta, data, raref)
            nest_cmd = ''
            print('len(self.nodfs_vars):', len(self.nodfs_vars))
            print("    self.nodfs_vars :", self.nodfs_vars)
            for ndx, nodf_var in enumerate(self.nodfs_vars):
                print()
                print()
                print('         nodf_var   : [', ndx, ']', nodf_var)
                comp = '%s/comparisons_%s.csv' % (m_dir, nodf_var)
                nest_cmd += self.write_comparisons(biom, comp, mode, nodf_var)
                if not to_do(biom):
                    io_update(self, i_f=biom, o_f=comp, key=data.dat)
                else:
                    io_update(self, o_f=comp, key=data.dat)
                for null in self.nulls:
                    simul = '%s/simulations_%s_%s.csv' % (m_dir, nodf_var, null)
                    nest_cmd += self.write_simulations(biom, simul, comp, null)
                    i_f = []
                    if not to_do(biom):
                        i_f.append(biom)
                    if not to_do(comp):
                        i_f.append(comp)
                    io_update(self, i_f=i_f, o_f=simul, key=data.dat)
            self.merge_outputs(m_dir, level, raref, cohort, nodfs_fp)
            self.res.append([dat, nodfs_fp, cohort, self.out])
            print()
            print('len(nest_cmd)', len(nest_cmd))
            if nest_cmd:
                self.write_meta(meta, data.metadata, meta_pd, nodfs)
                cmd += init_cmd + nest_cmd
        self.cmds.setdefault(data.dat, []).append(cmd)

    def figure_nodfs(self):
        self.analysis = 'nestedness_nodfs'
        self.res_pd = pd.DataFrame(self.res, columns=['dataset', 'nodfs_fps',
                                                      'cohort', 'o_dir'])
        nodfs_py = '%s/%s.py' % (RESOURCES, self.analysis)
        for dat, data in self.project.datasets.items():
            if data.taxon or not data.collapsed:
                continue
            dat_pd = self.res_pd.loc[self.res_pd.dat == dat]
            o_dir = dat_pd.loc[dat_pd['cohort'] == 'ALL'].o_dir.values[0]
            pdf, py = '%s/nodfs.pdf' % o_dir, '%s/nodfs.py' % o_dir
            if self.config.force or to_do(pdf):
                cmd = 'python3 %s\n' % py
                self.cmds.setdefault(dat, []).append(cmd)
                io_update(self, i_f=([py] + dat_pd['nodfs_fps']),
                          o_d=o_dir, key=dat)
                with open(rep(py), 'w') as o, open(rep(nodfs_py)) as f:
                    for line in f:
                        line_edit = line
                        if '<DAT>' in line:
                            line_edit = line_edit.replace('<DAT>', dat)
                        if '<ODIR>' in line:
                            line_edit = line_edit.replace(
                                '<ODIR>', self.get_path(o_dir))
                        if '<NODFS>' in line:
                            ns = [self.get_path(x) for x in dat_pd['nodfs_fps']]
                            line_edit = line_edit.replace(
                                "'<NODFS>'", str(ns))
                        if '<COLLAPSED>' in line:
                            line_edit = line_edit.replace(
                                "'<COLLAPSED>'", str(data.collapsed))
                        o.write(line_edit)
        self.register_io_command()

    def figure_graphs(self):
        self.analysis = 'nestedness_graphs'
        graphs_pd = pd.DataFrame(self.graphs, columns=[
            'dat', 'dataset', 'cohort', 'raref', 'graph'])
        graphs_py = '%s/%s.py' % (RESOURCES, self.analysis)
        for dat, data in self.project.datasets.items():
            tax = self.project.datasets[data.source].taxa[-1]
            level = self.get_taxo_level(dat, data)
            dat_pd = graphs_pd.loc[
                graphs_pd.dat == dat, ['dat', 'cohort', 'raref', 'graph']]
            for (dat_notax, cohort, raref, graph_fp) in dat_pd.values:
                if to_do(graph_fp):
                    continue
                pdf = '%s.pdf' % splitext(graph_fp)[0]
                txt = '%s.txt' % splitext(graph_fp)[0]
                py = '%s.py' % splitext(graph_fp)[0]
                if self.config.force or to_do(pdf) and to_do(txt):
                    cmd = 'python3 %s\n' % py
                    self.cmds.setdefault(dat, []).append(cmd)
                    io_update(self, i_f=[data.tsv[raref], data.meta, graph_fp],
                              o_f=[pdf, txt], key=dat)
                    with open(rep(py), 'w') as o, open(rep(graphs_py)) as f:
                        for line in f:
                            o.write(line.replace(
                                '<DAT>', dat).replace(
                                # '<RAREF>', data.rarefs[idx]).replace(
                                '<RAREF>', raref
                            ).replace(
                                '<TAB_FP>', self.get_path(data.tsv[raref])
                            ).replace(
                                '<META_FP>', self.get_path(data.meta)
                            ).replace(
                                "'<COL_SAMPLE>'", str(self.colors['sample'])
                            ).replace(
                                "'<COL_FEATURE>'", str(self.colors['feature'])
                            ).replace(
                                '<TAX_DAT>', dat_notax).replace(
                                '<TAXA_FP>', tax).replace(
                                '<LEVEL>', level).replace(
                                '<COHORT>', cohort).replace(
                                "'<COLL>'", str(data.collapsed)
                            ).replace(
                                '<GRAPH_FP>', self.get_path(graph_fp)))
        self.register_io_command()

    def run(self):
        self.get_config()
        for dat, data in self.project.datasets.items():
            level = self.get_taxo_level(dat, data)
            for raref, qza in data.qza.items():
                for group, group_vals in data.sample_subsets.items():
                    for vals in group_vals:
                        cohort = get_cohort('nestedness', group, vals, data)
                        if not cohort:
                            continue
                        self.get_output(data.path, cohort)
                        meta_pd = subset_meta(
                            data.metadata, vals, group, col=group)
                        if meta_pd.shape[0] < 10:
                            continue
                        print()
                        print()
                        print('========================')
                        print('dat:', dat)
                        print('level:', level)
                        print('raref:', raref)
                        print('group:', group)
                        print('vals:', vals)
                        print('meta_pd.shape:', meta_pd.shape)
                        print('- - - - - - - - - - - - - - - - - - - - - - ')
                        print('meta_pd.nunique():\n', meta_pd.nunique())
                        print('- - - - - - - - - - - - - - - - - - - - - - ')
                        print('========================')
                        self.scripts(meta_pd, data, raref, level, cohort)
        self.register_io_command()

    def register_io_command(self):
        AnalysisPrep.analyses_ios[self.analysis] = dict(self.ios)
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.ios = {}
        self.cmds = {}
