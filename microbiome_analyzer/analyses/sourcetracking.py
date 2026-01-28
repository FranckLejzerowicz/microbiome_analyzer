# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import sys
import yaml
import pkg_resources

from os.path import isdir
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import write_sourcetracking
from microbiome_analyzer._io_utils import subset_meta
from microbiome_analyzer._scratch import rep

RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources")


class Sourcetracking(object):

    def __init__(self, config, project):
        self.config = config
        self.project = project
        self.out = ''
        self.analysis = 'sourcetracking'
        self.dir = '%s/%s' % (project.dir, self.analysis)
        self.params_yml = '%s/params.yml' % self.dir
        self.ios = {}
        self.cmds = {}
        self.params = {
            'methods': ['sourcetracker', 'feast'],
            'estimators': [
                'RandomForestClassifier',
                'GradientBoostingClassifier',
                'AdaBoostClassifier',
                'LinearSVC'
            ],
            'times': 1,
            'size': 1,
            'em_iterations': 1000,
            'prevalence': None,
            'abundance': None,
            'order': 'meta-filter',
            'run': True,
            'loo': False,
            'diff_src': False,
            'alpha1': 0.001,
            'alpha2': 0.001,
            'beta': 10,
            'coverage': None,
            'source_rarefaction_depth': 0,
            'sink_rarefaction_depth': 0,
            'restarts': 10,
            'draws_per_restart': 1,
            'burnin': 100,
            'delay': 10,
            'cluster_start_delay': 25
        }
        self.source_sink = {}
        self.valids = 1
        self.run()

    def get_output(self, dat, cohort=''):
        self.out = '%s/%s' % (self.dir, dat)
        if cohort:
            cohort = cohort.replace('(', '').replace(')', '').replace(
                ' ', '_').replace(',', '_')
            self.out = (self.out + '/' + cohort).rstrip('/')
        if not isdir(rep(self.out)):
            os.makedirs(rep(self.out))

    @staticmethod
    def check_numerics(p, v):
        if p in ['prevalence', 'abundance'] and v:
            if str(v).isdigit():
                if not int(v) > 0:
                    sys.exit('[sourcetracking] "%s" invalid int (must >0)' % p)
            else:
                try:
                    if not 1 > float(v) > 0:
                        sys.exit('[sourcetracking] "%s" not in range [0-1)' % p)
                except ValueError:
                    sys.exit('[sourcetracking] "%s" not a float' % p)
        elif p in ['beta', 'source_rarefaction_depth', 'sink_rarefaction_depth',
                   'restarts', 'draws_per_restart', 'burnin', 'delay',
                   'cluster_start_delay', 'iterations', 'coverage'] and v:
            if not str(v).isdigit() or not int(v) > 0:
                sys.exit('[sourcetracking] "%s" invalid int (must be >0)' % p)
        elif p in ['alpha1', 'alpha2']:
            try:
                if not 1 > float(v) > 0:
                    sys.exit('[sourcetracking] "%s" not in range [0-1)' % p)
            except ValueError:
                sys.exit('[sourcetracking] "%s" not a float' % p)

    def get_sourcetracking_config(self):
        if 'params' in self.config.sourcetracking:
            for p, v in self.config.sourcetracking['params'].items():
                if p == 'order':
                    vs = ['meta-filter', 'filter-meta']
                    if v not in vs:
                        sys.exit('[sourcetracking] "%s" not in %s' % (p, vs))
                elif p == 'methods':
                    vs = ['sourcetracker', 'feast', 'q2_classifier']
                    if set(v).difference(vs):
                        sys.exit('[sourcetracking] "%s" not in %s' % (p, vs))
                elif p == 'estimators':
                    vs = ['RandomForestClassifier', 'AdaBoostClassifier',
                          'GradientBoostingClassifier', 'LinearSVC']
                    if set(v).difference(vs):
                        sys.exit('[sourcetracking] "%s" not in %s' % (p, vs))
                elif p in ['loo', 'run', 'diff_src'] and v not in [True, False]:
                    sys.exit('[sourcetracking] "%s" not a boolean')
                else:
                    self.check_numerics(p, v)
                self.params[p] = v
        if 'sourcesink' not in self.config.sourcetracking:
            print('At least one sink for one metadata column must be set '
                  '(no "sourcesink" in %s)' % self.config.sourcetracking_fp)
        else:
            self.source_sink = self.config.sourcetracking['sourcesink']

    def check(self, data, source_sink) -> (str, list, str, dict):
        column = source_sink['column']
        sources = ['']
        if 'source' in source_sink:
            sources = source_sink['source']
        if column not in data.metadata.columns:
            raise IOError('[%s] "%s" not in metadata' % (self.analysis, column))
        sink = source_sink['sink']
        if sink not in set(data.metadata[column]):
            raise IOError('All sinks "%s" not in metadata column "%s"' % (
                sink, column))
        if sources != [''] and not set(sources).issubset(
                set(data.metadata[column].unique())):
            raise IOError('All sources "%s" not in metadata column "%s"' % (
                sources, column))

        reps = {sink: sink.replace(
            '/', '').replace(
            '(', '').replace(
            ')', '').replace(
            ' ', '_')}
        for source in sources:
            reps.update({source: source.replace(
                '/', '').replace(
                '(', '').replace(
                ')', '').replace(
                ' ', '_')})
        sink = reps[sink]
        sources = [reps[source] for source in sources]
        return column, sources, sink, reps

    def get_outs(self, method) -> (str, bool):
        miss = False
        dir_method = self.out + '/' + method
        if method == 'q2_classifier':
            for root, dirs, files in os.walk(dir_method):
                if len(root.split(dir_method)[-1].split('/')) == 4:
                    print(method, root.split(dir_method)[-1].split('/'))
                    if 'predictions.tsv' not in files:
                        print('\n'.join(files))
                        miss = True
            outs = dir_method + '/t0/r*/*/predictions.tsv'
        elif method == 'feast':
            outs = dir_method + '/t0/*/out.r0*'
        elif method == 'sourcetracker':
            if self.params.get('loo'):
                outs = dir_method + '/t0/loo/mixing_proportions.txt'
            else:
                outs = dir_method + '/t0/r0/mixing_proportions.txt'
            for root, dirs, files in os.walk(dir_method):
                if len(root.split(dir_method)[-1].split('/')) == 3:
                    print(method, root.split(dir_method)[-1].split('/'))
                    if 'mixing_proportions.txt' not in files:
                        print('\n'.join(files))
                        miss = True
        else:
            raise IOError('Wring method name ("%s")' % method)
        return outs, miss

    # @staticmethod
    # def get_odir(cohort, name, var, snks, srcs):
    #     var_snk_src = '%s-%s/%s' % (var, snks, '_'.join(srcs))
    #     o_dir = '/'.join([cohort, name, var_snk_src])
    #     o_dir = o_dir.rstrip('/')
    #     return o_dir

    def run(self):
        self.get_sourcetracking_config()
        for dat, data in self.project.datasets.items():
            for raref, qza in data.qza.items():
                sourcetrackings = {}
                for cohort, (sams, variables) in data.subsets[raref].items():
                    for name, srcs_snks in self.source_sink.items():
                        for meth in self.params['methods']:
                            var, srcs, snk, reps = self.check(data, srcs_snks)
                            odir = '/'.join([cohort, name])
                            self.get_output(data.path, odir)
                            new_meta = '%s/meta.tsv' % self.out
                            new_qza = '%s/tab.qza' % self.out
                            new_tsv = '%s/tab.tsv' % self.out
                            meta_pd = subset_meta(
                                data.metadata, sams, variables, var)
                            meta_pd.replace({var: reps}, inplace=True)
                            meta_pd.to_csv(rep(new_meta), index=False, sep='\t')
                            vs = set(meta_pd.columns)
                            if var not in vs or snk not in set(meta_pd[var]):
                                continue
                            outs, miss = self.get_outs(meth)
                            sourcetrackings.setdefault((
                                cohort, name, meth), []).append((outs,))
                            if self.config.force or not glob.glob(outs) or miss:
                                cmd = write_sourcetracking(
                                    self, dat, qza, new_qza, new_tsv, new_meta,
                                    meth, var, snk, srcs)
                                self.cmds.setdefault(dat, []).append(cmd)
                data.sourcetracking[raref] = sourcetrackings
        self.register_io_command()
        self.write_params()

    def write_params(self):
        with open(rep(self.params_yml), 'w') as fpo:
            yaml.dump(self.params, fpo, default_flow_style=False)

    def register_io_command(self):
        AnalysisPrep.analyses_ios[self.analysis] = dict(self.ios)
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.ios = {}
        self.cmds = {}
