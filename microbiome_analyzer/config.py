# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from microbiome_analyzer.core.config import PrepConfig


def configure(**kwargs):
    """
    Main qiime2 functions writer.
    """
    configs = PrepConfig(**kwargs)
    configs.create_config()
    configs.show_config()
