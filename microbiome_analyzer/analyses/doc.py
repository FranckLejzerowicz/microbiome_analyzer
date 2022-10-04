# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import random
import pkg_resources
from os.path import isdir

from microbiome_analyzer.core.analysis import AnalysisPrep
from microbiome_analyzer.core.commands import run_export
from microbiome_analyzer._io_utils import subset_meta
from microbiome_analyzer._scratch import io_update, to_do, rep


RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources/r_scripts")


class DOC(object):

    def __init__(self, config, project):
        self.config = config
        self.project = project
        self.dir = project.dir
        self.out = ''
        self.dirs = project.dirs
        self.analysis = 'doc'
        self.ios = {}
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
        self.out = '%s/doc/%s' % (self.dir, dat)
        if cohort:
            self.out = (self.out + '/' + cohort).rstrip('/')
        if not isdir(rep(self.out)):
            os.makedirs(rep(self.out))

    def get_path(self, path_):
        path = path_
        params = self.config.run_params[self.analysis]
        if not self.config.jobs or not params['scratch']:
            path = rep(path_)
        return path

    def get_doc_config(self):
        if 'params' in self.config.doc:
            for param, value in self.config.doc['params'].items():
                self.params[param] = value
        if 'do_r' in self.config.doc and self.config.doc['do_r']:
            self.do_r = 1

    def write_tsv(self, dat: str, qza: str, new_meta: str,
                  new_qza: str, new_tsv: str) -> str:
        cmd = ''
        if to_do(new_tsv):
            cmd += '\nqiime feature-table filter-samples'
            cmd += ' --i-table %s' % qza
            cmd += ' --m-metadata-file %s' % new_meta
            cmd += ' --o-filtered-table %s\n' % new_qza
            cmd += run_export(new_qza, new_tsv, 'FeatureTable')
            cmd += 'rm %s %s\n' % (new_qza, new_qza.replace('.qza', '.biom'))
            io_update(self, i_f=[qza, new_meta], o_f=new_tsv, key=dat)
        else:
            io_update(self, i_f=new_tsv, key=dat)
        return cmd

    def write_doc_python(self, dat: str, new_meta: str, new_tsv: str,
                         mp_dir: str, mp_tsv: str) -> str:
        cmd = ''
        if self.params['use_mp']:
            cmd += 'rm -rf %s\n' % mp_dir
            cmd += 'mkdir -p %s\n' % mp_dir
            cmd += 'cp %s %s\n' % (new_tsv, mp_tsv)
            io_update(self, i_f=new_tsv, key=dat)

        cmd += '\nXDOC'
        if self.params['use_mp']:
            cmd += ' --i-otu %s' % mp_tsv
            cmd += ' --o-outdir %s' % mp_dir
            io_update(self, i_f=mp_tsv, o_d=mp_dir, key=dat)
        else:
            cmd += ' --i-otu %s' % new_tsv
            cmd += ' --o-outdir %s' % self.out
            io_update(self, i_f=new_tsv, o_d=self.out, key=dat)

        nodes = self.config.run_params[self.analysis]['nodes']
        cpus = self.config.run_params[self.analysis]['cpus']
        if self.params['use_mp']:
            cmd += ' --p-cpus %s' % (int(nodes) * int(cpus))
            cmd += ' --use_mp'
        else:
            cmd += ' --p-cpus 1'
        if self.params['replicate']:
            cmd += '-k %s' % self.params['replicate']
            cmd += '-m %s' % new_meta
        cmd += ' --p-r %s' % self.params['r']
        cmd += ' --p-subr %s' % self.params['subr']
        cmd += ' --p-mov-avg %s' % self.params['mov_avg']
        # for ci in self.params['ci']:
        #     cmd += ' --p-ci %s' % ci
        cmd += ' --p-span %s' % self.params['span']
        cmd += ' --p-degree %s' % self.params['degree']
        cmd += ' --p-family %s' % self.params['family']
        cmd += ' --p-iterations %s' % self.params['iterations']
        cmd += ' --p-surface %s' % self.params['surface']
        cmd += ' --p-nulls %s' % self.params['nulls']
        if self.params['non_zero']:
            cmd += ' --non-zero'
        if self.params['null']:
            cmd += ' --null'
        cmd += ' --verbose\n\n'

        if self.params['use_mp']:
            cmd = 'rsync -r %s/ %s\n' % (mp_dir, self.out)
            cmd += 'rm -rf %s\n\n' % mp_dir
        return cmd

    def write_doc_r(self, new_tsv, doc, pdf):
        r_script = []
        for line in self.template:
            r_script.append(line.replace(
                '<TAB_TSV>', self.get_path(new_tsv)).replace(
                '<DOC_TSV>', self.get_path(doc)).replace(
                '<R_DIR>', self.get_path(self.out)).replace(
                '<PDF>', self.get_path(pdf)))
        return r_script

    def run_python(self):
        self.analysis = 'doc'
        for dat, data in self.project.datasets.items():
            for raref, qza in data.qza.items():
                docs = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    self.get_output(data.path, cohort)
                    meta = '%s/meta.tsv' % self.out
                    new_qza = '%s/tab.qza' % self.out
                    new_tsv = '%s/tab.tsv' % self.out
                    token = [str(random.choice(range(100))) for x in range(3)]
                    mp_dir = '%s/doc_%s' % (self.dir, ''.join(token))
                    mp_tsv = '%s/tab.tsv' % mp_dir
                    docs.setdefault(cohort, []).append((self.out,))
                    if self.config.force or to_do('%s/DO.tsv' % self.out):
                        meta_pd = subset_meta(data.metadata, sams, group)
                        meta_pd.to_csv(rep(meta), index=False, sep='\t')
                        cmd = self.write_tsv(dat, qza, meta, new_qza, new_tsv)
                        cmd += self.write_doc_python(
                            dat, new_qza, new_tsv, mp_dir, mp_tsv)
                        self.cmds.setdefault(dat, []).append(cmd)

    def run_r(self):
        self.analysis = 'doc_r'
        for dat, data in self.project.datasets.items():
            for raref, qza in data.qza.items():
                docs = {}
                for cohort, (sams, group) in data.subsets[raref].items():
                    self.get_output(data.path, (cohort + '/R'))
                    meta = '%s/meta.tsv' % self.out
                    new_qza = '%s/tab.qza' % self.out
                    new_tsv = '%s/tab.tsv' % self.out
                    log = '%s/log.error' % self.out
                    pdf = '%s/plot.pdf' % self.out
                    do = '%s/DO.tsv' % self.out
                    docs.setdefault(cohort, []).append((self.out,))
                    if self.config.force or to_do(pdf):
                        meta_pd = subset_meta(data.metadata, sams, group)
                        meta_pd.to_csv(rep(meta), index=False, sep='\t')
                        r = self.write_doc_r(new_tsv, do, pdf)
                        script = '%s/doc.R' % self.out
                        if r:
                            with open(rep(script), 'w') as o:
                                for line in r:
                                    o.write(line)
                        cmd = self.write_tsv(dat, qza, meta, new_qza, new_tsv)
                        cmd += 'echo "***" >> %s\n' % log
                        cmd += 'R -f %s --vanilla 2>> %s\n' % (script, log)
                        cmd += 'echo "end" >> %s\n' % log
                        self.cmds.setdefault(dat, []).append(cmd)
                        io_update(self, i_f=[meta, script],
                                  o_d=self.out, key=dat)
                data.doc[raref] = docs
        self.register_io_command()

    def register_io_command(self):
        AnalysisPrep.analyses_ios[self.analysis] = dict(self.ios)
        AnalysisPrep.analyses_commands[self.analysis] = dict(self.cmds)
        self.ios = {}
        self.cmds = {}
