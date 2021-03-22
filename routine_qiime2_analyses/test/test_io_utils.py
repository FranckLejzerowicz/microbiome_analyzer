# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest
import pkg_resources
from routine_qiime2_analyses._routine_q2_io_utils import (
    get_prjct_nm,
    read_yaml_file,
    get_run_params
)

TEST = pkg_resources.resource_filename(
    'routine_qiime2_analyses', 'test')
RESOURCES = pkg_resources.resource_filename(
    'routine_qiime2_analyses', 'resources')


class InitChecks(unittest.TestCase):

    def setUp(self) -> None:
        self.dummy_yaml_fp = '%s/dummy_yaml.yml' % TEST
        with open(self.dummy_yaml_fp, 'w') as o:
            o.write('a:\n  b:\n    - "x"\n    - "y"')
        self.no_dummy_yaml_fp = '%s/no_dummy_yaml.yml' % TEST

    def test_get_prjct_nm(self):

        project_name = 'abracadabra'
        expected = 'brcdbr'
        observed = get_prjct_nm(project_name)
        self.assertEqual(expected, observed)

        project_name = 'aeiouy'
        expected = 'aeiouy'
        observed = get_prjct_nm(project_name)
        self.assertEqual(expected, observed)

    def test_read_yaml_file(self):
        expected = {'a': {'b': ['x', 'y']}}
        observed = read_yaml_file(self.dummy_yaml_fp)
        self.assertEqual(expected, observed)

        expected = {}
        observed = read_yaml_file(self.no_dummy_yaml_fp)
        self.assertEqual(expected, observed)

    def tearDown(self) -> None:
        os.remove(self.dummy_yaml_fp)


class RunParams(unittest.TestCase):

    def setUp(self) -> None:
        # this self.run_params certainly contains defaults:
        self.run_params_fp = '%s/run_params.yml' % RESOURCES
        self.run_params = read_yaml_file(self.run_params_fp)
        # import:
        #   time: "4"
        #   n_nodes: "1"
        #   n_procs: "1"
        #   mem_num: "10"
        #   mem_dim: "gb"
        #   env: "qiime2-2020.2"
        self.conda_envs = {'qiime2-2020.2', 'a_conda_env'}

        self.update_time_param_fp = '%s/update_time_param_fp.yml' % TEST
        with open(self.update_time_param_fp, 'w') as o:
            o.write('import:\n  time: "10"\n')
        self.update_env_param_fp = '%s/update_env_params_fp.yml' % TEST
        with open(self.update_env_param_fp, 'w') as o:
            o.write('import:\n  time: "10"\n  env: "a_conda_env"\n')
        self.update_not_a_env_fp = '%s/update_not_a_env_fp.yml' % TEST
        with open(self.update_not_a_env_fp, 'w') as o:
            o.write('import:\n  env: "not_a_conda_env"\n')
        self.update_not_mem_dim_fp = '%s/update_not_mem_dim_fp.yml' % TEST
        with open(self.update_not_mem_dim_fp, 'w') as o:
            o.write('import:\n  mem_dim: "tb"\n')

    def test_get_run_params(self):
        expected = {'time': "10", 'n_nodes': "1", 'n_procs': "1",
                    'mem_num': "10", 'mem_dim': "gb", 'env': "qiime2-2020.2"}
        observed = get_run_params(self.update_time_param_fp, self.conda_envs)
        self.assertEqual(observed['import'], expected)

        expected = {'time': "10", 'n_nodes': "1", 'n_procs': "1",
                    'mem_num': "10", 'mem_dim': "gb", 'env': "a_conda_env"}
        observed = get_run_params(self.update_env_param_fp, self.conda_envs)
        self.assertEqual(observed['import'], expected)

        with self.assertRaises(SystemExit) as cm:
            get_run_params(self.update_not_a_env_fp, self.conda_envs)
        self.assertEqual(cm.exception.code, 1)

        expected = {'time': "4", 'n_nodes': "1", 'n_procs': "1",
                    'mem_num': "10", 'mem_dim': "gb", 'env': "qiime2-2020.2"}
        observed = get_run_params(self.update_not_mem_dim_fp, self.conda_envs)
        self.assertEqual(observed['import'], expected)

    def tearDown(self) -> None:
        os.remove(self.update_env_param_fp)
        os.remove(self.update_time_param_fp)
        os.remove(self.update_not_a_env_fp)
        os.remove(self.update_not_mem_dim_fp)


if __name__ == '__main__':
    unittest.main()
