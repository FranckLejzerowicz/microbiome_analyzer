# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
from os.path import isfile
from typing import TextIO


def run_xpbs(out_sh: str, out_pbs: str, job_name: str,
             qiime_env: str, time: str, n_nodes: str,
             n_procs: str, mem_num: str, mem_dim: str,
             written: int, single: str, o: TextIO = None) -> None:
    if written:
        xpbs_call(out_sh, out_pbs, job_name, qiime_env,
                  time, n_nodes, n_procs, mem_num, mem_dim)
        if single:
            if not o:
                print(single)
                print('[TO RUN] qsub', out_pbs)
            else:
                o.write('qsub %s\n' % out_pbs)
    else:
        os.remove(out_sh)
        if isfile(out_pbs):
            os.remove(out_pbs)


def xpbs_call(run_sh: str, run_pbs: str, prjct_nm: str,
              qiime_env: str, time: str, n_nodes: str,
              n_procs: str, mem_num: str, mem_dim: str) -> None:
    """
    Call the subprocess to run Xpbs on the current bash script.

    :param run_sh: input current bash script.
    :param run_pbs: output script with directives.
    :param prjct_nm: project/job name.
    :param qiime_env: conda environment.
    :param time: walltime in hours.
    :param n_nodes: number of nodes to use.
    :param n_procs: number of processors to use.
    :param mem_num: memory in number.
    :param mem_dim: memory dimension to the number.
    :return:
    """
    cmd = [
        'Xpbs', '-i', run_sh,
        '-j', prjct_nm,
        '-o', run_pbs,
        '-e', qiime_env,
        '-t', time,
        '-n', n_nodes,
        '-p', n_procs,
        '-M', mem_num, mem_dim,
        '--noq'
    ]
    subprocess.call(cmd)
