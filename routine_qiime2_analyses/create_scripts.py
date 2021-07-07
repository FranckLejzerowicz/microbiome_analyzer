# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import numpy as np
from os.path import splitext
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_prjct_anlss_nm, get_job_folder
)


class CreateScripts(object):

    def __init__(self, config, prjct_nm, run_params, filt_raref, slurm):
        self.config = config
        self.prjct_nm = prjct_nm
        self.run_params = run_params
        self.filt_raref = filt_raref
        self.slurm = slurm
        self.jobs_folders = {}
        self.shs = {}
        self.to_chunk = {}

    def xpbs_call(self, sh, pbs_slm, job, params):
        cmd = [
            'Xpbs',
            '-i', sh,
            '-o', pbs_slm,
            '-j', job,
            '-e', self.config.qiime_env,
            '-t', params['time'],
            '-n', params['n_nodes'],
            '-p', params['n_procs'],
            '-M', params['mem_num'], params['mem_dim'],
            '-c', self.config.chmod,
            '--noq'
        ]
        if self.slurm:
            cmd.append('--slurm')
        # if tmp:
        #     cmd.extend(['-T', tmp])
        if not self.config.loc:
            cmd.append('--no-loc')
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

        job = '%s.%s.%s%s' % (self.prjct_nm, nlss, dat, self.filt_raref)
        if idx != 'None':
            job = '%s_%s' % (job, idx)
        if self.config.jobs:
            self.xpbs_call(sh, pbs_slm, job, params)

        if single:
            if self.slurm:
                launcher = 'sbatch'
            else:
                launcher = 'qsub'
            if os.getcwd().startswith('/panfs'):
                pbs_slm = pbs_slm.replace(os.getcwd(), '')
            if not main_o:
                if self.config.jobs:
                    print('[TO RUN] %s' % launcher, pbs_slm)
                else:
                    print('[TO RUN] sh', sh)
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
                head_sh = '%s_chunk%s.sh' % (job_folder, idx)
                chunks[(head_sh, idx)] = sorted(keys)
        else:
            chunks = dict((('%s_chunk%s.sh' % (job_folder, idx), idx), [x])
                          for idx, x in enumerate(to_chunk))

        for (sh, idx), commands in chunks.items():
            if commands:
                with open(sh, 'w') as o:
                    for command in commands:
                        o.write(command)
                if self.slurm:
                    pbs_slm = '%s.slm' % splitext(sh)[0]
                    launcher = 'sbatch'
                else:
                    pbs_slm = '%s.pbs' % splitext(sh)[0]
                    launcher = 'qsub'
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
        if os.getcwd().startswith('/panfs'):
            to_run = to_run.replace(os.getcwd(), '')
        if self.config.jobs:
            print('[TO RUN]', sh_pbs, to_run)
        else:
            print('[TO RUN]', 'sh', to_run.replace('.pbs', '.sh'))

    def write_scripts(self, analyses_commands):
        self.get_jobs_folders(analyses_commands)
        self.get_shs(analyses_commands)
        for analysis, shs in self.shs.items():
            nlss = get_prjct_anlss_nm(analysis)
            if 'import' in analysis:
                params = self.run_params['import']
            else:
                params = self.run_params[analysis]
            main_sh = self.jobs_folders[analysis][0]
            with open(main_sh, 'w') as main_o:
                if not self.config.chunkit:
                    for sh, dats_commands in shs.items():
                        with open(sh, 'w') as o:
                            for dat, cmd in dats_commands:
                                o.write('echo "%s"\n' % cmd.replace('"', ''))
                                o.write('%s\n' % cmd)
                        if self.slurm:
                            pbs_slm = '%s.slm' % splitext(sh)[0]
                        else:
                            pbs_slm = '%s.pbs' % splitext(sh)[0]
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
            self.jobs_folders[analysis] = [
                '%s/run_%s_%s%s.sh' % (
                    get_job_folder(self.config.i_datasets_folder, analysis),
                    analysis, self.prjct_nm, self.filt_raref),
                '%s/run_%s_%s%s' % (
                    get_job_folder(
                        self.config.i_datasets_folder, '%s/chunks' % analysis),
                    analysis, self.prjct_nm, self.filt_raref),
            ]
