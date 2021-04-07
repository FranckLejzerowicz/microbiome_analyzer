# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import isfile
import re
import unittest
import pandas as pd
import pkg_resources
from pathlib import Path
from pandas.testing import assert_frame_equal
from routine_qiime2_analyses._routine_q2_io_utils import (
    check_input,
    get_prjct_nm,
    get_filt_raref_suffix,
    read_yaml_file,
    get_run_params,
    get_data_paths,
    get_corresponding_meta,
    read_data_table,
    read_meta_pd,
    get_feature_sample_col,
    check_if_not_dna,
    collect_genome_id,
    check_features,
    genome_id_or_dna,
)

TEST = pkg_resources.resource_filename(
    'routine_qiime2_analyses', 'test')
RESOURCES = pkg_resources.resource_filename(
    'routine_qiime2_analyses', 'resources')


class InitChecks(unittest.TestCase):

    def setUp(self) -> None:

        self.no_dummy_folder = '%s/no_dummy_folder' % TEST
        self.dummy_folder = '%s/dummy_folder' % TEST
        os.makedirs(self.dummy_folder)
        self.dummy_folder_fp = '%s/dummy_folder_fp' % TEST
        with open(self.dummy_folder_fp, 'w') as _:
            pass

        self.dummy_yaml_fp = '%s/dummy_yaml.yml' % TEST
        with open(self.dummy_yaml_fp, 'w') as o:
            o.write('a:\n  b:\n    - "x"\n    - "y"')
        self.no_dummy_yaml_fp = '%s/no_dummy_yaml.yml' % TEST

    def test_check_input(self):
        with self.assertRaises(SystemExit) as cm:
            check_input(self.no_dummy_folder)
        self.assertEqual(cm.exception.code, 1)

        with self.assertRaises(SystemExit) as cm:
            check_input(self.dummy_folder_fp)
        self.assertEqual(cm.exception.code, 1)

        observed = check_input(self.dummy_folder)
        expected = '%s/dummy_folder' % TEST
        self.assertEqual(observed, expected)

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
        observed = read_yaml_file(None)
        self.assertEqual(expected, observed)

    def test_get_filt_raref_suffix(self):
        observed = get_filt_raref_suffix(None, False)
        self.assertEqual(observed, '')
        observed = get_filt_raref_suffix('', False)
        self.assertEqual(observed, '')
        observed = get_filt_raref_suffix('any string', False)
        self.assertEqual(observed, '_flt')
        observed = get_filt_raref_suffix(None, True)
        self.assertEqual(observed, '_rrf')
        observed = get_filt_raref_suffix('', True)
        self.assertEqual(observed, '_rrf')
        observed = get_filt_raref_suffix('any string', True)
        self.assertEqual(observed, '_flt_rrf')

    def tearDown(self) -> None:
        os.remove(self.dummy_yaml_fp)
        os.remove(self.dummy_folder_fp)
        os.rmdir(self.dummy_folder)


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


class DiscoverDatasets(unittest.TestCase):

    def setUp(self) -> None:
        self.i_datasets_folder = '%s/dataset_discovery' % TEST
        self.data = '%s/dataset_discovery/data' % TEST
        self.metadata = '%s/dataset_discovery/metadata' % TEST
        os.makedirs(self.i_datasets_folder)
        os.makedirs(self.data)
        os.makedirs(self.metadata)
        self.tsv = '%s/dataset_discovery/data/tab_a.tsv' % TEST
        self.biom = '%s/dataset_discovery/data/tab_a.biom' % TEST
        self.meta_tsv = '%s/dataset_discovery/metadata/meta_a.tsv' % TEST
        self.meta_txt = '%s/dataset_discovery/metadata/meta_a.txt' % TEST

    def test_get_data_paths_tsv(self):
        expected = {'a': self.tsv}
        Path(self.tsv).touch()
        observed = get_data_paths(('a',), self.i_datasets_folder)
        os.remove(self.tsv)
        self.assertEqual(observed, expected)

    def test_get_data_paths_biom(self):
        expected = {'a': self.biom}
        Path(self.biom).touch()
        observed = get_data_paths(('a',), self.i_datasets_folder)
        os.remove(self.biom)
        self.assertEqual(observed, expected)

    def test_get_data_paths_empty(self):
        with self.assertRaises(SystemExit) as cm:
            get_data_paths(('a',), self.i_datasets_folder)
        self.assertEqual(cm.exception.code, 0)

    def test_get_corresponding_meta(self):
        Path(self.tsv).touch()
        Path(self.meta_tsv).touch()
        observed = get_corresponding_meta(self.tsv)
        os.remove(self.meta_tsv)
        self.assertEqual(self.meta_tsv, observed)

        Path(self.meta_txt).touch()
        observed = get_corresponding_meta(self.tsv)
        os.remove(self.meta_txt)
        self.assertEqual(self.meta_txt, observed)

        with self.assertRaises(SystemExit) as cm:
            get_corresponding_meta(self.tsv)
        self.assertEqual(cm.exception.code, 0)
        os.remove(self.tsv)

    def tearDown(self) -> None:
        os.rmdir(self.data)
        os.rmdir(self.metadata)
        os.rmdir(self.i_datasets_folder)


class ReadDatasets(unittest.TestCase):

    def setUp(self) -> None:
        self.i_datasets_folder = '%s/read_datasets' % TEST
        self.metadata = '%s/read_datasets/metadata' % TEST
        self.data = '%s/read_datasets/data' % TEST
        os.makedirs(self.i_datasets_folder)
        os.makedirs(self.metadata)
        os.makedirs(self.data)
        self.tsv = '%s/read_datasets/data/tab_a.tsv' % TEST
        self.biom = '%s/read_datasets/data/tab_a.biom' % TEST
        self.meta_tsv = '%s/read_datasets/metadata/meta_a.tsv' % TEST
        self.meta_txt = '%s/read_datasets/metadata/meta_a.txt' % TEST

    def test_get_data_paths_tsv(self):
        expected = {'a': self.tsv}
        Path(self.tsv).touch()
        observed = get_data_paths(('a',), self.i_datasets_folder)
        os.remove(self.tsv)
        self.assertEqual(observed, expected)

    def test_get_data_paths_biom(self):
        expected = {'a': self.biom}
        Path(self.biom).touch()
        observed = get_data_paths(('a',), self.i_datasets_folder)
        os.remove(self.biom)
        self.assertEqual(observed, expected)

    def test_get_data_paths_empty(self):
        with self.assertRaises(SystemExit) as cm:
            get_data_paths(('b',), self.i_datasets_folder)
        self.assertEqual(cm.exception.code, 0)

    def test_get_corresponding_meta(self):
        Path(self.tsv).touch()
        Path(self.meta_tsv).touch()
        observed = get_corresponding_meta(self.tsv)
        self.assertEqual(self.meta_tsv, observed)
        os.remove(self.meta_tsv)

        Path(self.meta_txt).touch()
        observed = get_corresponding_meta(self.tsv)
        self.assertEqual(self.meta_txt, observed)
        os.remove(self.meta_txt)

        with self.assertRaises(SystemExit) as cm:
            get_corresponding_meta(self.tsv)
        self.assertEqual(cm.exception.code, 0)
        os.remove(self.tsv)

    def tearDown(self) -> None:
        os.rmdir(self.data)
        os.rmdir(self.metadata)
        os.rmdir(self.i_datasets_folder)


class ReadDatasets(unittest.TestCase):

    def setUp(self) -> None:
        self.i_datasets_folder = '%s/read_datasets' % TEST
        self.metadata = '%s/read_datasets/metadata' % TEST
        self.data = '%s/read_datasets/data' % TEST
        os.makedirs(self.i_datasets_folder)
        os.makedirs(self.metadata)
        os.makedirs(self.data)
        self.tsv = '%s/read_datasets/data/tab_a.tsv' % TEST
        self.biom = '%s/read_datasets/data/tab_a.biom' % TEST
        self.meta_tsv = '%s/read_datasets/metadata/meta_a.tsv' % TEST
        self.meta_txt = '%s/read_datasets/metadata/meta_a.txt' % TEST

    def test_get_data_paths_tsv(self):
        expected = {'a': self.tsv}
        Path(self.tsv).touch()
        observed = get_data_paths(('a',), self.i_datasets_folder)
        os.remove(self.tsv)
        self.assertEqual(observed, expected)

    def test_get_data_paths_biom(self):
        expected = {'a': self.biom}
        Path(self.biom).touch()
        observed = get_data_paths(('a',), self.i_datasets_folder)
        os.remove(self.biom)
        self.assertEqual(observed, expected)

    def test_get_data_paths_empty(self):
        with self.assertRaises(SystemExit) as cm:
            get_data_paths(('b',), self.i_datasets_folder)
        self.assertEqual(cm.exception.code, 0)

    def test_get_corresponding_meta(self):
        Path(self.tsv).touch()
        Path(self.meta_tsv).touch()
        observed = get_corresponding_meta(self.tsv)
        self.assertEqual(self.meta_tsv, observed)
        os.remove(self.meta_tsv)

        Path(self.meta_txt).touch()
        observed = get_corresponding_meta(self.tsv)
        self.assertEqual(self.meta_txt, observed)
        os.remove(self.meta_txt)

        with self.assertRaises(SystemExit) as cm:
            get_corresponding_meta(self.tsv)
        self.assertEqual(cm.exception.code, 0)
        os.remove(self.tsv)

    def tearDown(self) -> None:
        os.rmdir(self.data)
        os.rmdir(self.metadata)
        os.rmdir(self.i_datasets_folder)


class ReadDataTable(unittest.TestCase):

    def setUp(self) -> None:
        self.data = '%s/read_data/data' % TEST
        os.makedirs(self.data)

        self.tsv_a = '%s/read_data/data/tab_a.tsv' % TEST
        with open(self.tsv_a, 'w') as o:
            o.write('index\tA\tB\n')
            o.write('s1\t1\t2\n')
        self.tsv_b = '%s/read_data/data/tab_b.tsv' % TEST
        with open(self.tsv_b, 'w') as o:
            o.write('first_columns\tA\tB\n')
            o.write('s1\t1\t2\n')
        self.tsv_c = '%s/read_data/data/tab_c.tsv' % TEST
        Path(self.tsv_c).touch()

    def test_read_data_table(self):
        expected = pd.DataFrame({'A': [1], 'B': [2]}, index=['s1'])
        expected.index.name = '#OTU ID'
        observed = read_data_table(self.tsv_a)
        assert_frame_equal(observed, expected)

    def test_read_meta_pd(self):
        expected = pd.DataFrame({'sample_id': ['s1'], 'A': [1], 'B': [2]})
        observed = read_meta_pd(self.tsv_a, 'sample_id')
        assert_frame_equal(observed, expected)

        expected = pd.DataFrame({'any_id': ['s1'], 'A': [1], 'B': [2]})
        observed = read_meta_pd(self.tsv_a, 'any_id')
        assert_frame_equal(observed, expected)

    def test_get_feature_sample_col(self):
        observed = get_feature_sample_col(self.tsv_a)
        self.assertEqual(observed, 'index')

        observed = get_feature_sample_col(self.tsv_b)
        self.assertEqual(observed, 'first_columns')

        with self.assertRaises(SystemExit) as cm:
            get_feature_sample_col(self.tsv_c)
        self.assertEqual(cm.exception.code, 0)

    def tearDown(self) -> None:
        os.remove(self.tsv_a)
        os.remove(self.tsv_b)
        os.remove(self.tsv_c)
        os.rmdir(self.data)


class FeatureDetection(unittest.TestCase):

    def setUp(self) -> None:
        self.not_dna = re.compile('[^ACGTN].*?')
        self.tsv = '%s/a.tsv' % TEST
        self.genome_ids = {'G000000001': 'taxon1', 'G000000002': 'taxon2'}
        self.datasets_features = {'dat': ['G000000001', 'G000000002']}
        self.tab = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]},
            index=['a;1;2__G000000001', 'b; 1; 2__G000000002'])
        self.tab1 = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]},
            index=['tax[a;1;2]', 'b;1;2'])
        self.tab2 = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]},
            index=['a; 1; 2', 'b; 1; 2'])

    def test_check_if_not_dna(self):
        observed = check_if_not_dna('ACGT', True, self.not_dna)
        self.assertTrue(observed)
        observed = check_if_not_dna('ACXGT', True, self.not_dna)
        self.assertFalse(observed)
        observed = check_if_not_dna('ACXGT', False, self.not_dna)
        self.assertFalse(observed)
        # at this point (False), the sequence is not checked
        observed = check_if_not_dna('ACGT', False, self.not_dna)
        self.assertFalse(observed)

    def test_collect_genome_id_false(self):
        datasets_features, datasets_read, datasets_phylo = {}, {}, {}
        collect_genome_id(
            'dat', self.tsv, self.tab1, self.genome_ids,
            datasets_features, datasets_read, datasets_phylo, False)
        self.assertEqual(datasets_features, {'dat': self.genome_ids})
        self.assertEqual(datasets_read, {})
        self.assertEqual(datasets_phylo, {'dat': ('wol', 0)})

    def test_collect_genome_id_1(self):
        tab = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]},
            index=['tax[a|1|2]', 'b|1|2'])
        tab_written = pd.DataFrame(
            {'index': ['tax[a|1|2]', 'b|1|2'], 'A': [1, 2], 'B': [3, 4]})
        datasets_features, datasets_phylo = {}, {}
        datasets_read = {'dat': [pd.DataFrame(), None]}
        collect_genome_id(
            'dat', self.tsv, self.tab1, self.genome_ids,
            datasets_features, datasets_read, datasets_phylo, True)
        self.assertEqual(datasets_features, {'dat': self.genome_ids})
        assert_frame_equal(datasets_read['dat'][0], tab)
        self.assertEqual(datasets_phylo, {'dat': ('wol', 1)})
        assert_frame_equal(self.tab1, tab)
        assert_frame_equal(pd.read_table(self.tsv), tab_written)
        os.remove(self.tsv)

    def test_collect_genome_id_2(self):
        tab = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]}, index=['a|1|2', 'b|1|2'])
        tab_written = pd.DataFrame(
            {'index': ['a|1|2', 'b|1|2'], 'A': [1, 2], 'B': [3, 4]})
        datasets_features, datasets_phylo = {}, {}
        datasets_read = {'dat': [pd.DataFrame(), None]}
        collect_genome_id(
            'dat', self.tsv, self.tab2, self.genome_ids,
            datasets_features, datasets_read, datasets_phylo, True)
        self.assertEqual(datasets_features, {'dat': self.genome_ids})
        assert_frame_equal(datasets_read['dat'][0], tab)
        self.assertEqual(datasets_phylo, {'dat': ('wol', 1)})
        assert_frame_equal(self.tab2, tab)
        assert_frame_equal(pd.read_table(self.tsv), tab_written)
        os.remove(self.tsv)

    def test_check_features(self):
        data_pd = pd.DataFrame(
            {'A': [1, 2]},
            index=['AAA', 'CCC'])
        obs_gids, obs_dna, obs_correct = check_features(data_pd)
        self.assertEqual(obs_gids, {})
        self.assertTrue(obs_dna)
        self.assertFalse(obs_correct)

        data_pd.index = ['AAA', 'BBB']
        obs_gids, obs_dna, obs_correct = check_features(data_pd)
        self.assertEqual(obs_gids, {})
        self.assertFalse(obs_dna)
        self.assertFalse(obs_correct)

        data_pd.index = ['AAA__G000000001', 'CCC__G000000002']
        expected_gids = {'G000000001': 'AAA__G000000001',
                         'G000000002': 'CCC__G000000002'}
        obs_gids, obs_dna, obs_correct = check_features(data_pd)
        self.assertEqual(obs_gids, expected_gids)
        self.assertFalse(obs_dna)
        self.assertFalse(obs_correct)

        data_pd.index = ['A;A;A__G000000001', 'C; C; C__G000000002']
        expected_gids = {'G000000001': 'A|A|A__G000000001',
                         'G000000002': 'C|C|C__G000000002'}
        obs_gids, obs_dna, obs_correct = check_features(data_pd)
        self.assertEqual(obs_gids, expected_gids)
        self.assertFalse(obs_dna)
        self.assertTrue(obs_correct)

    def test_genome_id_or_dna_is_numeric(self):
        num_dtype_index = pd.DataFrame({'A': [1, 2], 'B': [3, 4]})
        datasets_features, datasets_phylo, datasets_read = {}, {}, {}
        genome_id_or_dna('dat', self.tsv, num_dtype_index, datasets_read,
                         datasets_features, datasets_phylo)
        self.assertEqual(datasets_features, {})
        self.assertEqual(datasets_phylo, {})
        self.assertEqual(datasets_read, {})

    def test_genome_id_or_dna_is_dna(self):
        dna = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]},
            index=['ACGT', 'TGCA'])
        datasets_features, datasets_phylo, datasets_read = {}, {}, {}
        genome_id_or_dna('dat', self.tsv, dna, datasets_read,
                         datasets_features, datasets_phylo)
        self.assertEqual(datasets_features, {})
        self.assertEqual(datasets_phylo, {'dat': ('amplicon', 0)})
        self.assertEqual(datasets_read, {})
        self.assertFalse(isfile(self.tsv))

    def test_genome_id_or_dna_is_gid(self):
        expected = pd.DataFrame(
            {'A': [1, 2], 'B': [3, 4]},
            index=['a|1|2__G000000001', 'b|1|2__G000000002'])
        expected_written = pd.DataFrame({
            'index': ['a|1|2__G000000001', 'b|1|2__G000000002'],
            'A': [1, 2], 'B': [3, 4]})
        datasets_features, datasets_phylo = {}, {}
        datasets_read = {'dat': [pd.DataFrame(), None]}
        genome_id_or_dna('dat', self.tsv, self.tab, datasets_read,
                         datasets_features, datasets_phylo)
        assert_frame_equal(datasets_read['dat'][0], expected)
        self.assertEqual(datasets_phylo, {'dat': ('wol', 1)})
        assert_frame_equal(self.tab , expected)
        assert_frame_equal(pd.read_table(self.tsv), expected_written)
        os.remove(self.tsv)


if __name__ == '__main__':
    unittest.main()
