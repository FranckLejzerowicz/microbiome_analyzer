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
    get_train_test_dict,
)

TEST = pkg_resources.resource_filename(
    'routine_qiime2_analyses', 'test')
RESOURCES = pkg_resources.resource_filename(
    'routine_qiime2_analyses', 'resources')


class GetTrainTest(unittest.TestCase):

    def test_get_train_test_dict(self):
        config_fp = '%s/config.yml' % TEST
        observed = get_train_test_dict(config_fp)
        self.assertEqual(observed, {'train': 0.7})

        with open(config_fp, 'w') as o:
            o.write('datasets:\n')
            o.write('  name:\n')
            o.write('    split_on:\n')
            o.write('      - "sex"\n')
        observed = get_train_test_dict(config_fp)
        self.assertEqual(observed, {
            'train': 0.7, 'datasets': {'name': {'split_on': ['sex']}}
        })

        with open(config_fp, 'w') as o:
            o.write('train: 0.2\n')
        observed = get_train_test_dict(config_fp)
        self.assertEqual(observed, {'train': 0.2})

        with open(config_fp, 'w') as o:
            o.write('train: -0.2\n')
        observed = get_train_test_dict(config_fp)
        self.assertEqual(observed, {'train': 0.7})

        with open(config_fp, 'w') as o:
            o.write('train: 1.2\n')
        observed = get_train_test_dict(config_fp)
        self.assertEqual(observed, {'train': 0.7})

        os.remove(config_fp)


if __name__ == '__main__':
    unittest.main()
