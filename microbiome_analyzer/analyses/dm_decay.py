# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import pandas as pd

from os.path import isdir, splitext
from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import run_export
from microbiome_analyzer._io_utils import subset_meta
from microbiome_analyzer._scratch import io_update, to_do, rep

RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources/python_scripts")


class DmDecay(object):

    def __init__(self, config, project):
        self.config = config
        self.dir = project.dir
        self.dirs = project.dirs
        self.out = ''
        self.project = project
        self.analysis = 'dm_decay'
        self.ios = {}
        self.cmds = {}
        self.subsets = {'ALL': [[]]}
        self.modes = {'individual': ['']}
        self.step = 10
        self.iteration = 10
        self.res = []
        self.res_pd = pd.DataFrame
        self.group_source = ''
        self.group_target = ''
        self.val_source = ''
        self.val_target = ''
        # run stuff
        self.run()
        self.plot()
        AnalysisPrep.analyses_dm_decay = self.res_pd

    def get_output(self, dat):
        self.out = '%s/%s/%s' % (self.dir, self.analysis, dat)
        if not isdir(rep(self.out)):
            os.makedirs(rep(self.out))

    def get_config(self):
        dm_decay_config = self.config.dm_decay
        if 'modes' in dm_decay_config:
            self.modes.update(dm_decay_config['modes'])
        if 'params' in dm_decay_config:
            if 'step' in dm_decay_config['params']:
                self.step = dm_decay_config['params']['step']
            if 'iteration' in dm_decay_config['params']:
                self.iteration = dm_decay_config['params']['iteration']

    def write_dm_decay(
            self, dat: str, dm: str, dm_filt: str, qza: str,
            tsv: str, meta_fp: str, mode: str, mode_group: str) -> str:
        cmd = ''
        if to_do(qza):
            cmd += 'qiime diversity filter-distance-matrix '
            cmd += ' --i-distance-matrix %s ' % dm
            cmd += ' --m-metadata-file %s ' % meta_fp
            cmd += ' --o-filtered-distance-matrix %s\n' % dm_filt
            cmd += 'qiime distance-decay %s ' % mode
            cmd += ' --i-distance-matrix %s ' % dm_filt
            if mode != 'individual':
                cmd += ' --m-metadata-file %s ' % meta_fp
                if 'targeted' in mode:
                    cmd += ' --p-source-category %s ' % self.group_source
                    cmd += ' --p-target-category %s ' % self.group_target
                    cmd += ' --p-source-category-value %s ' % self.val_source
                    cmd += ' --p-target-category-value %s ' % self.val_target
                else:
                    cmd += ' --m-category-column %s ' % mode_group
                cmd += ' --p-balance '
            cmd += ' --p-iteration %s ' % self.iteration
            cmd += ' --p-step %s ' % self.step
            params = self.config.run_params['dm_decay']
            cmd += ' --p-jobs %s ' % (int(params['nodes']) *
                                         int(params['cpus']))
            cmd += ' --o-distance-decay %s\n' % qza
            io_update(self, i_f=[dm, meta_fp], key=dat)
            cmd += 'rm %s\n' % dm_filt

        if to_do(tsv):
            cmd += run_export(qza, tsv, 'decay')
            cmd += 'rm %s\n' % qza
            if not to_do(qza):
                io_update(self, i_f=qza, o_f=tsv, key=dat)
            else:
                io_update(self, o_f=tsv, key=dat)

        return cmd

    def dm_decay(
            self,
            dat: str,
            full_meta: pd.DataFrame,
            meta: pd.DataFrame,
            dm: str,
            metric: str,
            raref: str
    ):
        meta = full_meta.set_index('sample_name').loc[meta.sample_name.tolist()]
        for mode, mode_groups in self.modes.items():
            for mode_group in mode_groups:
                if mode == 'individual':
                    m_dir = '%s/%s/%s' % (self.out, mode, metric)
                    mode_pd = meta.iloc[:, :2]
                elif 'targeted' in mode:
                    # soon implemented
                    continue
                else:
                    m_dir = '%s/%s_%s/%s' % (self.out, mode, mode_group, metric)
                    if mode_group not in set(meta.columns):
                        continue
                    if meta[mode_group].unique().size == 1:
                        continue
                    if min(meta[mode_group].value_counts()) < 25:
                        continue
                    if str(meta[mode_group].dtype) != 'object':
                        continue
                    mode_pd = meta[[mode_group]]
                if not isdir(rep(m_dir)):
                    os.makedirs(rep(m_dir))
                meta_fp = '%s_meta%s.tsv' % (m_dir, raref)
                # mode_pd.columns = ['#SampleID'] + mode_pd.columns.tolist()[1:]
                mode_pd.to_csv(rep(meta_fp), sep='\t')
                dm_filt = '%s.qza' % m_dir
                qza = '%s_decay.qza' % m_dir
                tsv = '%s_decay.tsv' % m_dir
                self.res.append([dat, raref, metric, mode, mode_group, tsv])
                if self.config.force or to_do(tsv):
                    cmd = self.write_dm_decay(
                        dat, dm, dm_filt, qza, tsv, meta_fp, mode, mode_group)
                    self.cmds.setdefault(dat, []).append(cmd)

    def run(self):
        self.analysis = 'dm_decay'
        self.get_config()
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                for cohort, (sams, group) in data.subsets[raref].items():
                    self.get_output(data.path)
                    for ddx, (dm, _, me) in enumerate(dms_metrics):
                        if to_do(dm):
                            continue
                        meta = subset_meta(data.metadata, sams, group)
                        if meta.shape[0] < 10:
                            continue
                        self.dm_decay(dat, data.metadata, meta, dm, me, raref)
        self.register_io_command()

    def get_path(self, path_):
        path = path_
        params = self.config.run_params[self.analysis]
        if not self.config.jobs or not params['scratch']:
            path = rep(path_)
        return path

    def plot(self):
        self.analysis = 'dm_decay_plot'
        dm_decay_plot_py = '%s/%s.py' % (RESOURCES, self.analysis)
        self.res_pd = pd.DataFrame(self.res, columns=[
            'dat', 'raref', 'metric', 'mode', 'mode_group', 'tsv'])
        for dat, dat_pd in self.res_pd.groupby('dat'):
            self.get_output(dat)
            dat_fpo = '%s/res.tsv' % self.out
            fig_o = '%s/decays.pdf' % self.out
            dat_pd.to_csv(rep(dat_fpo), index=False, sep='\t')
            py = '%s.py' % splitext(dat_fpo)[0]
            if self.config.force or to_do(fig_o):
                cmd = 'python3 %s\n' % py
                self.cmds.setdefault(dat, []).append(cmd)
                io_update(self, i_f=[py, dat_fpo], o_f=[fig_o], key=dat)
                with open(rep(py), 'w') as o, open(rep(dm_decay_plot_py)) as f:
                    for line in f:
                        o.write(line.replace(
                            '<FIG_O>', self.get_path(fig_o)).replace(
                            '<DATASET>', dat).replace(
                            '<DAT_FP>', self.get_path(dat_fpo)))
        self.register_io_command()

    def register_io_command(self):
        AnalysisPrep.analyses_ios[self.analysis] = dict(self.ios)
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.ios = {}
        self.cmds = {}
