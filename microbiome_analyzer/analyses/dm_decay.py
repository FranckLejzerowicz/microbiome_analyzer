# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import pandas as pd

from os.path import isdir, isfile, splitext
from microbiome_analyzer.analysis import AnalysisPrep
from microbiome_analyzer._cmds import run_export
from microbiome_analyzer._io import subset_meta, get_output

RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources/python_scripts")


class DmDecay(object):

    def __init__(self, config, project):
        self.config = config
        self.project = project
        self.analysis = 'dm_decay'
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
        out = '%s/qiime/%s/%s' % (self.config.folder, self.analysis, dat)
        if not isdir(out):
            os.makedirs(out)
        return out

    def get_config(self):
        dm_decay_config = self.config.dm_decay
        if 'modes' in dm_decay_config:
            self.modes.update(dm_decay_config['modes'])
        if 'params' in dm_decay_config:
            if 'step' in dm_decay_config['params']:
                self.step = dm_decay_config['params']['step']
            if 'iteration' in dm_decay_config['params']:
                self.iteration = dm_decay_config['params']['iteration']

    def write_dm_decay(self, dm: str, dm_filt: str, qza: str, tsv: str,
                       meta_fp: str, mode: str, mode_group: str) -> str:
        cmd = ''
        if not isfile(qza):
            cmd += 'qiime diversity filter-distance-matrix \\\n'
            cmd += '--m-metadata-file %s \\\n' % meta_fp
            cmd += '--i-distance-matrix %s \\\n' % dm
            cmd += '--o-filtered-distance-matrix %s\n' % dm_filt

            cmd += 'qiime distance-decay %s \\\n' % mode
            cmd += '--i-distance-matrix %s \\\n' % dm_filt
            if mode != 'individual':
                cmd += '--m-metadata-file %s \\\n' % meta_fp
                if 'targeted' in mode:
                    cmd += '--p-source-category %s \\\n' % self.group_source
                    cmd += '--p-target-category %s \\\n' % self.group_target
                    cmd += '--p-source-category-value %s \\\n' % self.val_source
                    cmd += '--p-target-category-value %s \\\n' % self.val_target
                else:
                    cmd += '--m-category-column %s \\\n' % mode_group
                cmd += '--p-balance \\\n'
            cmd += '--p-iteration %s \\\n' % self.iteration
            cmd += '--p-step %s \\\n' % self.step
            params = self.config.run_params['dm_decay']
            cmd += '--p-jobs %s \\\n' % (int(params['n_nodes']) *
                                         int(params['n_procs']))
            cmd += '--o-distance-decay %s\n' % qza
            cmd += 'rm %s\n' % dm_filt
        if not isfile(tsv):
            cmd += run_export(qza, tsv, 'decay')
            cmd += 'rm %s\n' % qza
        return cmd

    def dm_decay(self, dat: str, full_meta: pd.DataFrame, meta: pd.DataFrame,
                 dm: str, metric: str, raref: str, o_dir: str):

        meta = full_meta.set_index('sample_name').loc[meta.sample_name.tolist()]
        for mode, mode_groups in self.modes.items():
            for mode_group in mode_groups:
                if mode == 'individual':
                    m_dir = '%s/%s/%s' % (o_dir, mode, metric)
                    mode_pd = meta.iloc[:, :2]
                elif 'targeted' in mode:
                    # soon implemented
                    continue
                else:
                    m_dir = '%s/%s_%s/%s' % (o_dir, mode, mode_group, metric)
                    if mode_group not in set(meta.columns):
                        continue
                    if meta[mode_group].unique().size == 1:
                        continue
                    if min(meta[mode_group].value_counts()) < 25:
                        continue
                    if str(meta[mode_group].dtype) != 'object':
                        continue
                    mode_pd = meta[[mode_group]]
                if not isdir(m_dir):
                    os.makedirs(m_dir)
                meta_fp = '%s_meta%s.tsv' % (m_dir, raref)
                # mode_pd.columns = ['#SampleID'] + mode_pd.columns.tolist()[1:]
                mode_pd.to_csv(meta_fp, sep='\t')
                dm_filt = '%s.qza' % m_dir
                qza = '%s_decay.qza' % m_dir
                tsv = '%s_decay.tsv' % m_dir
                self.res.append([dat, raref, metric, mode, mode_group, tsv])
                if self.config.force or not isfile(tsv):
                    cmd = self.write_dm_decay(
                        dm, dm_filt, qza, tsv, meta_fp, mode, mode_group)
                    self.cmds.setdefault(dat, []).append(cmd)

    def run(self):
        self.analysis = 'dm_decay'
        self.get_config()
        for dat, data in self.project.datasets.items():
            for raref, dms_metrics in data.beta.items():
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path)
                    for ddx, (dm, _, me) in enumerate(dms_metrics):
                        if not isfile(dm):
                            continue
                        meta = subset_meta(data.metadata, sams, group)
                        if meta.shape[0] < 10:
                            continue
                        self.dm_decay(dat, data.metadata, meta,
                                      dm, me, raref, o_dir)
        self.register_command()

    def plot(self):
        self.analysis = 'dm_decay_plot'
        dm_decay_plot_py = '%s/%s.py' % (RESOURCES, self.analysis)
        self.res_pd = pd.DataFrame(self.res, columns=[
            'dat', 'raref', 'metric', 'mode', 'mode_group', 'tsv'])
        for dat, dat_pd in self.res_pd.groupby('dat'):
            o_dir = get_output(self.config.folder, 'dm_decay/%s' % dat)
            dat_fpo = '%s/res.tsv' % o_dir
            fig_o = '%s/decays.pdf' % o_dir
            dat_pd.to_csv(dat_fpo, index=False, sep='\t')
            py = '%s.py' % splitext(dat_fpo)[0]
            if self.config.force or not isfile(fig_o):
                cmd = 'python3 %s\n' % py
                self.cmds.setdefault(dat, []).append(cmd)
                with open(py, 'w') as o, open(dm_decay_plot_py) as f:
                    for line in f:
                        o.write(line.replace(
                            '<FIG_O>', fig_o).replace(
                            '<DATASET>', dat).replace(
                            '<DAT_FP>', dat_fpo))
        self.register_command()

    def register_command(self):
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.cmds = {}
