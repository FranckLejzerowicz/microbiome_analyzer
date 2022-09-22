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
from os.path import isdir, splitext


class CreateScripts(object):

    def __init__(self, config):
        self.config = config
        self.jobs_folders = {}
        self.to_chunk = {}
        self.shs = {}
        self.sh = None

    def xpbs_call(self, sh, pbs_slm, job, params):
        # mandatory options
        cmd = [
            'Xhpc',
            '-i', sh,
            '-o', pbs_slm,
            '-j', job,
            '-e', self.config.qiime_env,
            '-t', params['time'],
            '-n', params['n_nodes'],
            '-c', params['n_procs'],
            '-M', params['mem_num'], params['mem_dim'],
            '--no-stat', '--quiet'
        ]
        # always provide an account
        if self.config.account:
            cmd.extend(['-a', self.config.account])
        # get the job script file path and use it
        if self.config.torque:
            cmd.append('--torque')
        # machine-specific setup: env activating and slurm partition
        if params['machine']:
            cmd.append('--%s' % params['machine'])
        if params['partition']:
            cmd.append('--p-partition %s' % params['partition'])
        # setup the scratch location to be used for the current software
        if isinstance(params['scratch'], int):
            cmd.extend(['--localscratch', str(params['scratch'])])
        elif params['scratch'] == 'scratch':
            cmd.append('--scratch')
        elif params['scratch'] == 'userscratch':
            cmd.append('--userscratch')

        subprocess.call(cmd)

        if os.getcwd().startswith('/panfs'):
            pbs_lines = open(pbs_slm).readlines()
            with open(pbs_slm, 'w') as pbs:
                for pbs_line in pbs_lines:
                    pbs.write(pbs_line.replace(os.getcwd(), ''))

    def run_xpbs(self, sh, pbs_slm, nlss, params, dat,
                 idx='None', single=False, main_o=None):

        if os.getcwd().startswith('/panfs'):
            sh_lines = open(sh).readlines()
            with open(sh, 'w') as o:
                if not self.config.jobs:
                    o.write('conda activate %s\n' % self.config.qiime_env)
                for sh_line in sh_lines:
                    o.write(sh_line.replace(os.getcwd(), ''))

        job = '%s.%s.%s%s' % (self.config.prjct_nm, nlss,
                              dat, self.config.filt_raref)
        if idx != 'None':
            job = '%s_%s' % (job, idx)
        if self.config.jobs:
            self.xpbs_call(sh, pbs_slm, job, params)

        if single:
            if self.config.torque:
                launcher = 'qsub'
            else:
                launcher = 'sbatch'
            if os.getcwd().startswith('/panfs'):
                pbs_slm = pbs_slm.replace(os.getcwd(), '')
            if not main_o:
                if self.config.jobs:
                    print('%s' % launcher, pbs_slm)
                else:
                    print('sh', sh)
            else:
                if self.config.jobs:
                    main_o.write('%s %s\n' % (launcher, pbs_slm))
                else:
                    main_o.write('sh %s\n' % sh)

    def do_chunks(self, job_folder, to_chunk, nlss, params, main_o):
        chunks = {}
        if len(to_chunk) > self.config.chunkit:
            array_split = np.array_split(to_chunk, self.config.chunkit)
            for idx, keys in enumerate(array_split):
                head_sh = '%s_%s_chunk%s.sh' % (
                    job_folder, self.config.prjct_nm, idx)
                chunks[(head_sh, idx)] = sorted(keys)
        else:
            chunks = dict((('%s_%s_chunk%s.sh' % (
                job_folder, self.config.prjct_nm, idx), idx), [x])
                          for idx, x in enumerate(to_chunk))

        for (sh, idx), commands in chunks.items():
            if commands:
                with open(sh, 'w') as o:
                    for command in commands:
                        o.write(command)
                if self.config.torque:
                    pbs_slm = '%s.pbs' % splitext(sh)[0]
                    launcher = 'qsub'
                else:
                    pbs_slm = '%s.slm' % splitext(sh)[0]
                    launcher = 'sbatch'
                if self.config.jobs:
                    self.run_xpbs(sh, pbs_slm, nlss, params, 'chnk', str(idx))
                    if os.getcwd().startswith('/panfs'):
                        pbs_slm = pbs_slm.replace(os.getcwd(), '')
                    main_o.write('%s %s\n' % (launcher, pbs_slm))
                else:
                    main_o.write('sh %s\n' % sh)

    def print_message(self, message, sh_pbs, to_run):
        if message:
            print('#', message)
        if os.getcwd().startswith('/panfs') or os.getcwd().startswith('${SCRA'):
            to_run = to_run.replace(os.getcwd(), '')
        if self.config.jobs:
            print(sh_pbs, to_run)
        else:
            print('sh', to_run.replace('.slm', '.sh'))

    def get_prjct_anlss_nm(self, project_name: str) -> str:
        """
        Get a smaller name for printing in qstat / squeue.

        Parameters
        ----------
        project_name : str
            Command-line passed project name.

        Returns
        -------
        prjct_nm : str
            Same name without the vows ("aeiouy").
        """
        alpha = 'aeiouy'
        prjct_nm = ''.join(x for x in project_name if x.lower() not in alpha)
        if prjct_nm == '':
            prjct_nm = project_name
        return prjct_nm

    def write_scripts(self, analyses_commands):
        self.get_jobs_folders(analyses_commands)
        self.get_shs(analyses_commands)
        for analysis, shs in self.shs.items():
            nlss = self.get_prjct_anlss_nm(analysis)
            params = self.config.run_params.get(
                'import', self.config.run_params['default'])
            main_sh = self.jobs_folders[analysis][0]
            with open(main_sh, 'w') as main_o:
                if not self.config.chunkit:
                    for sh, dats_commands in shs.items():
                        with open(sh, 'w') as o:
                            for dat, cmd in dats_commands:
                                # o.write('echo "%s"\n' % cmd.replace('"', ''))
                                o.write('%s\n' % cmd)
                        if self.config.torque:
                            pbs_slm = '%s.pbs' % splitext(sh)[0]
                        else:
                            pbs_slm = '%s.slm' % splitext(sh)[0]
                        self.run_xpbs(sh, pbs_slm, nlss, params,
                                      dat, 'None', True, main_o)
                else:
                    job_folder = self.jobs_folders[analysis][1]
                    to_chunk = self.to_chunk[analysis]
                    self.do_chunks(job_folder, to_chunk, nlss, params, main_o)
            self.print_message(analysis, 'sh', main_sh)

    def get_shs(self, analyses_commands):
        for analysis, dats_commands in analyses_commands.items():
            job_folder = self.jobs_folders[analysis][1]
            for dat, commands in dats_commands.items():
                sh = '%s_%s.sh' % (job_folder, dat)
                for idx, command in enumerate(commands):
                    if analysis not in self.shs:
                        self.shs[analysis] = {}
                        self.to_chunk[analysis] = []
                    self.shs[analysis].setdefault(sh, []).append((dat, command))
                    self.to_chunk[analysis].append(command)

    def get_jobs_folders(self, analyses_commands):
        for analysis in analyses_commands:
            analysis_dir = '%s/%s' % (self.config.output_folder, analysis)
            main_sh = '%s/run_%s%s.sh' % (analysis_dir, self.config.prjct_nm,
                                          self.config.filt_raref)
            jobs_dir = '%s/jobs' % analysis_dir
            job_chunk = '%s/chunks_%s%s' % (jobs_dir, self.config.prjct_nm,
                                           self.config.filt_raref)
            if not isdir(jobs_dir):
                os.makedirs(jobs_dir)
            self.jobs_folders[analysis] = [main_sh, job_chunk]
