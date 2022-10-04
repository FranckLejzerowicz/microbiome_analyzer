# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pkg_resources
from microbiome_analyzer.core.config import AnalysesConfig

OUT = pkg_resources.resource_filename("microbiome_analyzer", "tests")
RESOURCE = pkg_resources.resource_filename("microbiome_analyzer", "resources")



class ConfigTests(unittest.TestCase):

    def setUp(self):

        config_dict = {
            'datasets': ('dat1', 'dat2',),
            'datasets_folder': '',
            'project_name': 'test',
            'qiime_env': 'qiime2-2021.11',
            'classifier': '%s/gg-13-8-99-nb-classifier.qza' % RESOURCE,
            'sepp_tree': '%s/sepp-refs-gg-13-8.qza' % RESOURCE,
            'wol_tree': '%s/wol/wol_tree.nwk' % RESOURCE,
            'filter_fp': 'example/filtering.yml',
            'filt3d_fp': 'example/filtering_3d.yml',
            'rarefs_fp': 'example/rarefactions.yml',
            'feature_subsets_fp': 'example/feature_subsets.yml',
            'sample_subsets_fp': 'example/sample_subsets.yml',
            'collapse_fp': 'example/collapse.yml',
            'procrustes_fp': 'example/procrustes.yml',
            'mantel_fp': 'example/procrustes.yml',
            'dm_decay_fp': 'example/dm_decay.yml',
            'geo_decay_fp': False,
            'run_params_fp': None,
            'qemistree': None,
            'adonis_fp': 'example/adonis.yml',
            'beta_type': ('permanova', 'permdisp',),
            'tests': ('sex', 'age_cat',),
            'nestedness_fp': 'example/nestedness.yml',
            'phate_fp': 'example/phate.yml',
            'doc_fp': 'example/doc.yml',
            'sourcetracking_fp': 'example/sourcetracking.yml',
            'train_test_fp': 'example/train_test.yml',
            'diff_models_fp': 'example/songbird.yml',
            'mmvec_pairs_fp': 'example/mmvec.yml',
            'mmvec_highlights_fp': 'example/mmvec_highlight.yml',
            'xmmvec_fp': 'example/mmbird.yml',
            'longi_column': 'timepoint',
            'alphas': None,
            'betas': None,
            'biplot': True,
            'chmod': '664',
            'chunks': None,
            'filt_only': False,
            'force': True,
            'gpu': False,
            'jobs': False,
            'loc': True,
            'raref': True,
            'skip': None,
            'slurm': False}
        config = AnalysesConfig(**config_dict)
        print(config.__dict__)

    def test_check_input(self):
        pass

    def test_get_prjct_anlss_nm(self):
        pass

    def test_parse_yamls(self):
        pass

    def test_get_run_params(self):
        pass

    def test_get_filt_raref_suffix(self):
        pass

    def test_get_train_test_dict(self):
        pass

    def test_check_xpbs_install(self):
        pass

    def test_get_workflow(self):
        pass


if __name__ == '__main__':
    unittest.main()
