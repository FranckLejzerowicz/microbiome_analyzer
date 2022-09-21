# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re, ast
from setuptools import find_packages, setup

classes = """
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""

classifiers = [s.strip() for s in classes.split('\n') if s]

description = (
    "microbiome_analyzer is a command line tool that writes commands to run one"
    "by one to perform analyses (mainly using qiime2) on a HPC running Slurm,"
    "Torque or none of these."
)

with open("README.md") as f:
    long_description = f.read()

_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open("microbiome_analyzer/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
    version = str(ast.literal_eval(hit))

standalone = ['microbiome_analyzer=microbiome_analyzer.scripts._standalone_analyzer:standalone_analyzer']

setup(
    name="microbiome_analyzer",
    version=version,
    license="BSD",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Franck Lejzerowicz",
    author_email="franck.lejzerowicz@gmail.com",
    maintainer="Franck Lejzerowicz",
    maintainer_email="franck.lejzerowicz@gmail.com",
    url="https://github.com/FranckLejzerowicz/microbiome_analyzer",
    packages=find_packages(),
    install_requires=[
        "click==7.1.2",
        "pandas",
        "scikit-bio==0.5.6",
        "numpy",
        "scipy",
        "pyyaml",
        "plotly==4.8.2",
        "phate",
        "biom-format",
        "seaborn"
    ],
    classifiers=classifiers,
    entry_points={'console_scripts': standalone},
    package_data={
        'routine_qiime2_analyses': [
            'test/*/*/*',
            'resources/run_params.yml',
            'resources/beta_metrics.txt',
            'resources/alpha_metrics.txt',
            'resources/wol/wol_tree.nwk',
            'resources/wol/g2lineage.txt',
            'resources/r_scripts/doc.R',
            'resources/r_scripts/adonis.R',
            'resources/sh_scripts/nestedness.sh',
            'resources/sh_scripts/spatial_autocorrelation_modeling.sh',
            'resources/python_scripts/nestedness_nodfs.py',
            'resources/python_scripts/nestedness_graphs.py',
            'resources/python_scripts/summarize_permanovas.py',
            'resources/python_scripts/mmvec_pre_paired-heatmaps.py',
            'examples/*'
        ],
    },
    include_package_data=True,
    python_requires='>=3.6',
)
