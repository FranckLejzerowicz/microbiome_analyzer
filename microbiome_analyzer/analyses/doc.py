# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import random
import pkg_resources

from os.path import isdir, isfile
from microbiome_analyzer.analysis import AnalysisPrep
from microbiome_analyzer._cmds import run_export
from microbiome_analyzer._io import subset_meta

RESOURCES = pkg_resources.resource_filename(
    "routine_qiime2_analyses", "resources/r_scripts")


class DOC(object):

    def __init__(self, config, project):
        self.config = config
        self.project = project
        self.analysis = 'doc'
        self.cmds = {}
        self.doc = {}
        self.template = open('%s/doc.R' % RESOURCES).readlines()
        self.do_r = 0
        self.params = {
            'r': 50,
            'subr': 0,
            'mov_avg': 5,
            'ci': ['0.025', '0.5', '0.975'],
            'span': 0.2,
            'degree': 1,
            'family': 'symmetric',
            'iterations': 2,
            'surface': 'direct',
            'nulls': 1,
            'non_zero': 1,
            'null': 1,
            'use_mp': False,
            'replicate': None}
        self.get_doc_config()
        # run
        self.run_python()
        if self.do_r:
            self.run_r()

    def get_output(self, dat, cohort=''):
        out = '%s/qiime/doc/%s' % (self.config.folder, dat)
        if cohort:
            out = (out + '/' + cohort).rstrip('/')
        if not isdir(out):
            os.makedirs(out)
        return out

    def get_doc_config(self):
        if 'params' in self.config.doc:
            for param, value in self.config.doc['params'].items():
                self.params[param] = value
        if 'do_r' in self.config.doc and self.config.doc['do_r']:
            self.do_r = 1

    @staticmethod
    def write_tsv(qza: str, new_meta: str,
                  new_qza: str, new_tsv: str) -> str:
        cmd = ''
        if not isfile(new_tsv):
            cmd += '\nqiime feature-table filter-samples \\\n'
            cmd += '--i-table %s \\\n' % qza
            cmd += '--m-metadata-file %s \\\n' % new_meta
            cmd += '--o-filtered-table %s\n' % new_qza

            cmd += run_export(new_qza, new_tsv, 'FeatureTable')
            cmd += 'rm %s %s\n' % (new_qza, new_qza.replace('.qza', '.biom'))
        return cmd

    def write_doc_python(self, new_meta: str, new_tsv: str, o_dir: str,
                         mp_dir: str, mp_tsv: str) -> str:
        cmd = ''
        if self.params['use_mp']:
            cmd += 'rm -rf %s\n' % mp_dir
            cmd += 'mkdir -p %s\n' % mp_dir
            cmd += 'cp %s %s\n' % (new_tsv, mp_tsv)

        cmd += '\nXDOC \\\n'
        if self.params['use_mp']:
            cmd += '--i-otu %s \\\n' % mp_tsv
            cmd += '--o-outdir %s \\\n' % mp_dir
        else:
            cmd += '--i-otu %s \\\n' % new_tsv
            cmd += '--o-outdir %s \\\n' % o_dir
        n_nodes = self.config.run_params[self.analysis]['n_nodes']
        n_procs = self.config.run_params[self.analysis]['n_procs']
        if self.params['use_mp']:
            cmd += '--p-cpus %s \\\n' % (int(n_nodes) * int(n_procs))
            cmd += '--use_mp \\\n'
        else:
            cmd += '--p-cpus 1 \\\n'
        if self.params['replicate']:
            cmd += '-k %s \\\n' % self.params['replicate']
            cmd += '-m %s \\\n' % new_meta
        cmd += '--p-r %s \\\n' % self.params['r']
        cmd += '--p-subr %s \\\n' % self.params['subr']
        cmd += '--p-mov-avg %s \\\n' % self.params['mov_avg']
        # for ci in self.params['ci']:
        #     cmd += '--p-ci %s \\\n' % ci
        cmd += '--p-span %s \\\n' % self.params['span']
        cmd += '--p-degree %s \\\n' % self.params['degree']
        cmd += '--p-family %s \\\n' % self.params['family']
        cmd += '--p-iterations %s \\\n' % self.params['iterations']
        cmd += '--p-surface %s \\\n' % self.params['surface']
        cmd += '--p-nulls %s \\\n' % self.params['nulls']
        if self.params['non_zero']:
            cmd += '--non-zero \\\n'
        if self.params['null']:
            cmd += '--null \\\n'
        cmd += '--verbose\n\n'

        if self.params['use_mp']:
            cmd = 'rsync -r %s/ %s\n' % (mp_dir, o_dir)
            cmd += 'rm -rf %s\n\n' % mp_dir
        return cmd

    def write_doc_r(self, new_tsv, r_dir, doc, pdf):
        r_script = []
        for line in self.template:
            r_script.append(line.replace(
                '<TAB_TSV>', new_tsv).replace(
                '<DOC_TSV>', doc).replace(
                '<R_DIR>', r_dir).replace(
                '<PDF>', pdf))
        return r_script

    def run_python(self):
        self.analysis = 'doc'
        for dat, data in self.project.datasets.items():
            for raref, qza in data.qza.items():
                docs = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    o_dir = self.get_output(data.path, cohort)
                    meta = '%s/meta.tsv' % o_dir
                    new_qza = '%s/tab.qza' % o_dir
                    new_tsv = '%s/tab.tsv' % o_dir
                    token = [str(random.choice(range(100))) for x in range(3)]
                    mp_dir = '%s/doc_%s' % (self.config.folder, ''.join(token))
                    mp_tsv = '%s/tab.tsv' % mp_dir
                    docs.setdefault(cohort, []).append((o_dir,))
                    if self.config.force or not isfile('%s/DO.tsv' % o_dir):
                        meta_pd = subset_meta(data.metadata, sams, group)
                        meta_pd.to_csv(meta, index=False, sep='\t')
                        cmd = self.write_tsv(qza, meta, new_qza, new_tsv)
                        cmd += self.write_doc_python(new_qza, new_tsv, o_dir,
                                                    mp_dir, mp_tsv)
                        self.cmds.setdefault(dat, []).append(cmd)

    def run_r(self):
        self.analysis = 'doc_R'
        for dat, data in self.project.datasets.items():
            for raref, qza in data.qza.items():
                docs = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    r_dir = self.get_output(data.path, (cohort + '/R'))
                    meta = '%s/meta.tsv' % r_dir
                    new_qza = '%s/tab.qza' % r_dir
                    new_tsv = '%s/tab.tsv' % r_dir
                    log = '%s/log.error' % r_dir
                    pdf = '%s/plot.pdf' % r_dir
                    do = '%s/DO.tsv' % r_dir
                    docs.setdefault(cohort, []).append((r_dir,))
                    if self.config.force or not isfile(pdf):
                        meta_pd = subset_meta(data.metadata, sams, group)
                        meta_pd.to_csv(meta, index=False, sep='\t')
                        r = self.write_doc_r(new_tsv, r_dir, do, pdf)
                        script = '%s/doc.R' % r_dir
                        if r:
                            with open(script, 'w') as o:
                                for line in r:
                                    o.write(line)
                        cmd = self.write_tsv(qza, meta, new_qza, new_tsv)
                        cmd += 'echo "***" >> %s\n' % log
                        cmd += 'R -f %s --vanilla 2>> %s\n' % (script, log)
                        cmd += 'echo "end" >> %s\n' % log
                        self.cmds.setdefault(dat, []).append(cmd)
                data.doc[raref] = docs
        self.register_command()

    def register_command(self):
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.cmds = {}
