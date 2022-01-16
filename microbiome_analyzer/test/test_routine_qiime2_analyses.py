# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest


class MyTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_email_address(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
