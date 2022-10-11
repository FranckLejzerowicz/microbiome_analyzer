# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import numpy as np
from os.path import dirname, isdir, splitext
from microbiome_analyzer._scratch import get_roundtrip, rep


class CreateScripts(object):

    def __init__(self, config, project):
        self.config = config
        self.project = project
        self.dir = project.dir
        self.jobs_folders = {}
        self.job_fps = []
        self.job_name = None
        self.jobs_dir = None
        self.tmpdir = None
        self.analysis = None
        self.nlss = None
        self.main = None
        self.to_chunk = {}
        self.params = {}
        self.chunks = {}
        self.run = {}
        self.cmds = {}
        self.ios = {}
        self.shs = {}
        self.sh = None
        self.cmd = []
        self.scripts = []
        self.scheduler = self.get_scheduler()

    def get_scheduler(self):
        if self.config.jobs:
            if self.config.torque:
                return 'qsub'
            else:
                return 'sbatch'
        else:
            return 'sh'

    def make_dirs(self):
        for directory in sorted(self.project.dirs):
            if not isdir(directory):
                os.makedirs(directory)

    def print_message(self, message, sh_pbs, to_run):
        if message:
            print('#', message)
        if os.getcwd().startswith('/panfs') or os.getcwd().startswith('${SCRA'):
            to_run = to_run.replace(os.getcwd(), '')
        if self.config.jobs:
            print(sh_pbs, to_run)
        else:
            print('sh', to_run.replace('.slm', '.sh'))

    def get_nlss_nm(self):
        """
        Get a smaller name for printing in qstat / squeue.
        """
        alpha = 'aeiouy'
        self.nlss = ''.join(x for x in self.analysis if x.lower() not in alpha)
        if self.nlss == '':
            self.nlss = self.analysis

    def scratch(self, dat, cmds):
        key = (dat,)
        if isinstance(dat, tuple):
            key = dat
        if self.params['scratch'] and self.config.jobs and key in self.ios:
            roundtrip = get_roundtrip(self.ios[key])
            extended_cmds = ['\n# Move to SCRATCH_FOLDER']
            extended_cmds += roundtrip['to']
            extended_cmds += ['\n# data: %s' % ' '.join(key)]
            extended_cmds += cmds
            if self.config.move_back:
                extended_cmds += ['\n# Move from SCRATCH_FOLDER']
                extended_cmds += roundtrip['from']
        else:
            extended_cmds = ['\n# data: %s' % ' '.join(key)]
            extended_cmds += [x.replace('${SCRATCH_FOLDER}', '') for x in cmds]
        self.cmds[key] = extended_cmds

    def get_cmds(self, cmds):
        self.cmds = {}
        for dat, dat_cmds in cmds.items():
            self.scratch(dat, dat_cmds)

    def get_chunks(self):
        self.chunks = {}
        n_chunks = self.params['chunks']
        if n_chunks and len(self.cmds) >= n_chunks:
            cmds_d = dict(enumerate(self.cmds))
            chunks = [list(x) for x in np.array_split(list(cmds_d), n_chunks)]
            for cdx, chunk in enumerate(chunks):
                if len(chunk):
                    self.chunks[str(cdx)] = [cmds_d[x] for x in chunk]
        else:
            for key in self.cmds:
                chunkey = '-'.join(key).replace('/', '_')
                self.chunks[chunkey] = [key]

    def get_jobs_dir(self):
        self.jobs_dir = '%s/%s/jobs' % (rep(self.dir), self.analysis)
        self.tmpdir = '%s/%s/jobs/tmp' % (rep(self.dir), self.analysis)
        if not isdir(self.jobs_dir):
            os.makedirs(self.jobs_dir)
        if not isdir(self.tmpdir):
            os.makedirs(self.tmpdir)

    def write_chunks(self, chunk_keys):
        with open(self.sh, 'w') as sh:
            sh.write('TMPDIR=%s\n' % self.tmpdir)
            cleanup = 'cleanup rm -rf ${TMPDIR}'
            if self.params['scratch'] and self.config.jobs:
                cleanup += ' ${SCRATCH_FOLDER}/*'
            sh.write('%s\n' % cleanup)
            for key in chunk_keys:
                sh.write('TMPDIR=%s/%s\n' % (self.tmpdir, key.replace(' ', '')))
                for cmd in self.cmds[key]:
                    sh.write('%s\n' % cmd)

    def get_job_name(self, chunk: str):
        self.job_name = self.nlss + '.' + self.config.prjct_nm
        self.job_name += '.' + chunk.replace('/', '_')

    def write_script(self):
        if self.config.jobs:
            self.prep_script()
            self.call_cmd()
        else:
            self.job_fps.append('%s.sh' % splitext(self.sh)[0])

    def prep_script(self) -> None:
        # mandatory options
        self.cmd = [
            'Xhpc',
            '-i', self.sh,
            '-j', self.job_name,
            '-t', str(self.params['time']),
            '-c', str(self.params['cpus']),
            '-M', str(self.params['mem']), self.params['mem_dim'],
            '--no-stat', '--no-abspath']
        # whether the cpus requests is per node
        if self.params['nodes']:
            self.cmd.extend(['-n', str(self.params['nodes'])])
        # always provide an account
        if self.config.account:
            self.cmd.extend(['-a', self.config.account])
        # get the job script file path and use it
        job_script = '%s.slm' % splitext(self.sh)[0]
        if self.config.torque:
            self.cmd.append('--torque')
            job_script = '%s.pbs' % splitext(self.sh)[0]
        self.job_fps.append(job_script)
        self.cmd.extend(['-o', job_script])
        # machine-specific setup: env activating and slurm partition
        if self.params['machine']:
            self.cmd.append('--%s' % self.params['machine'])
        if self.params['partition']:
            self.cmd.append('--p-partition %s' % self.params['partition'])
        # whether an environment must be used for the current software
        if self.params['env']:
            self.cmd.extend(['-e', self.params['env']])
        # setup the scratch location to be used for the current software
        if isinstance(self.params['scratch'], int):
            self.cmd.extend(['--localscratch', str(self.params['scratch'])])
        elif self.params['scratch'] == 'scratch':
            self.cmd.append('--scratch')
        elif self.params['scratch'] == 'userscratch':
            self.cmd.append('--userscratch')
        self.cmd.append('--quiet')

    def call_cmd(self):
        cmd = ' '.join(self.cmd)
        subprocess.call(cmd.split())
        os.remove(self.sh)

    def write_jobs(self):
        self.get_jobs_dir()
        self.get_nlss_nm()
        for chunk, chunk_keys in self.chunks.items():
            self.sh = '%s/run_%s.sh' % (self.jobs_dir, chunk.replace('/', '_'))
            self.write_chunks(chunk_keys)
            self.get_job_name(chunk)
            self.write_script()

    def get_main_sh(self, local=''):
        self.main = '%s/%s/run%s.sh' % (rep(self.dir), self.analysis, local)
        self.run['%s%s' % (self.analysis, local)] = self.main

    def write_main(self):
        self.get_main_sh()
        with open(self.main, 'w') as o:
            if len(self.job_fps):
                job_dir = dirname(self.job_fps[0])
                o.write('mkdir -p %s/output\n' % job_dir)
                o.write('cd %s/output\n' % job_dir)
                for jdx, job_fp in enumerate(self.job_fps):
                    o.write('%s %s\n' % (self.scheduler, job_fp))
                self.job_fps = []

    def writing(self, prep):
        m = max(len(x) for x in prep.analyses_commands.keys()) + 1
        for sdx, (analysis, cmds) in enumerate(prep.analyses_commands.items()):
            if not len(cmds):
                continue
            self.analysis = analysis
            self.params = self.config.run_params[analysis]
            self.ios = prep.analyses_ios[analysis]
            self.get_cmds(cmds)
            self.get_chunks()
            self.write_jobs()
            self.write_main()

    def display(self):
        print()
        soft_print = '========== microbiome analysis scripts ========== '
        print(soft_print)
        self.scripts.append(soft_print)
        for name, main in self.run.items():
            main_print = '>%s\nsh %s' % (name, main)
            print(main_print)
            self.scripts.append(main_print)
