# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess


def xpbs_call(run_sh, run_pbs, prjct_nm, qiime_env,
              time, n_nodes, n_procs, mem_num, mem_dim):
    cmd = [
        'Xpbs', '-i', run_sh,
        '-j', '%s.tree.mprt' % prjct_nm,
        '-o', run_pbs, '-e', qiime_env,
        '-t', str(time), '-n', str(n_nodes), '-p', str(n_procs),
        '-M', str(mem_num), str(mem_dim), '--noq'
    ]
    # print(' '.join(cmd))
    subprocess.call(cmd)
