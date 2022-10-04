# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from microbiome_analyzer import __version__

from microbiome_analyzer.scripts.modules import run
from microbiome_analyzer.scripts.modules import config
# from microbiome_analyzer.scripts.modules import export


@click.group(help="microbiome_analyzer commands")
@click.version_option(__version__, prog_name="microbiome_analyzer")
def standalone_analyzer():
    pass


standalone_analyzer.add_command(run.run)
standalone_analyzer.add_command(config.config)
# standalone_analyzer.add_command(export.export)


if __name__ == '__main__':
    standalone_analyzer()
