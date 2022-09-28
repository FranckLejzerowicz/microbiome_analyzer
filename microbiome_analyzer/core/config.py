# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import subprocess
import pkg_resources
from os.path import abspath, exists, isfile, isdir
from microbiome_analyzer._inputs import read_yaml_file

RESOURCES = pkg_resources.resource_filename("microbiome_analyzer", "resources")


class AnalysesConfig(object):
    """
    Collect the data associated with each dataset passed but the user
    """
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.prjct_nm = ''
        self.filt_raref = ''
        self.subsets = {}
        self.analyses = {}
        self.run_params = {}
        self.train_test_dict = {}
        self.conda_envs = set()

    def init(self):
        self.check_input()
        self.parse_yamls()
        self.get_project_name()  # project name for jobs
        self.get_conda_envs()
        self.get_analyses()
        self.get_run_params()  # default (e.g. memory)
        self.get_filt_raref_suffix()  # job suffix (e.g. _flt)
        self.get_train_test_dict()
        if self.jobs:
            self.check_xpbs_install()  # Xpbs must be installed to prepare jobs

    def parse_yamls(self):
        for arg in list(self.__dict__):
            if arg.endswith('_fp'):
                yaml = read_yaml_file(self.__dict__[arg])
                if 'subsets' in yaml:
                    yaml['subsets']['ALL'] = [[]]
                setattr(self, arg[:-3], yaml)

                for var, vals in yaml.get('subsets', {}).items():
                    vals_set = set(map(tuple, vals))
                    if var in self.subsets:
                        self.subsets[var].update(vals_set)
                    else:
                        self.subsets[var] = vals_set
        self.subsets['ALL'] = [[]]

    def check_input(self):
        """Check that the main input folder exists and that it is indeed a
        folder, not a file. Returns its absolute path for unambiguous use in
        the rest of the code.
        """
        if not exists(self.dir):
            sys.exit('%s does exist' % self.dir)
        if not isdir('%s/data' % self.dir):
            sys.exit('%s/data must be a folder' % self.dir)
        if not isdir('%s/metadata' % self.dir):
            sys.exit('%s/metadata must be a folder' % self.dir)
        self.dir = '${SCRATCH_FOLDER}%s' % abspath(self.dir)

    def get_project_name(self):
        """
        Get a smaller name for printing in qstat / squeue.
        """
        alpha = 'aeiouy'
        self.prjct_nm = ''.join(
            x for x in self.project_name if x.lower() not in alpha)
        if self.prjct_nm == '':
            self.prjct_nm = self.project_name

    def check_env(self, env: str = None) -> None:
        """Checks that the qiime2 environment
        exists (i.e., is in the set if discovered
        conda environments.

        Parameters
        ----------
        conda_envs : set
            Conda environments.
        env : str
            An environment.
        """
        conda_env = self.qiime_env
        if env:
            conda_env = env
        if conda_env not in self.conda_envs:
            print('%s is not an existing conda environment' % env)
            sys.exit(1)

    def get_conda_envs(self):
        """Get the names of the conda environments.

        Returns
        -------
        conda_envs : set
            Conda environments.
        """
        for conda_env in subprocess.getoutput('conda env list').split('\n'):
            if len(conda_env) and conda_env[0] != '#':
                self.conda_envs.add(conda_env.split()[0])
        # check that qiime2 environment exists
        self.check_env()

    def check_run_params(
            self,
            run_params_user: dict
    ):
        """Verify that the keys/values passed for
        an analysis are valid.

        Parameters
        ----------
        run_params_user : dict
            Passed run parameters.
        """
        mem_dim = ['kb', 'mb', 'gb']
        integer_values = ['time', 'nodes', 'cpus', 'mem']
        for analysis in list(run_params_user.keys()):
            for param in list(run_params_user[analysis].keys()):
                param_val = run_params_user[analysis][param]
                if param in integer_values:
                    if str(param_val).isdigit():
                        self.run_params[analysis][param] = param_val
                elif param == 'mem_dim':
                    if param_val in mem_dim:
                        self.run_params[analysis][param] = param_val
                elif param == 'env':
                    # check that qiime2 environment exists
                    self.check_env(param_val)
                    self.run_params[analysis][param] = param_val
                else:
                    self.run_params[analysis][param] = param_val

    def get_run_params(self):
        """Get the run parameters based on the default
        values, that are possibly updated by the user.
        """
        self.set_default_params()
        self.set_global_scratch()
        self.set_global_chunks()
        self.set_user_params()

    def set_default_params(self):
        # read default run parameters (one set for each analysis)
        run_params_fp = '%s/run_params.yml' % RESOURCES
        print('* Collecting defaults parameters from', run_params_fp)
        self.run_params = read_yaml_file(run_params_fp)

    def set_global_scratch(self):
        for (_, analysis) in self.analyses:
            if self.localscratch:
                self.run_params[analysis]['scratch'] = self.localscratch
            elif self.scratch:
                self.run_params[analysis]['scratch'] = 'scratch'
            elif self.userscratch:
                self.run_params[analysis]['scratch'] = 'userscratch'

    def set_global_chunks(self):
        print(self.__dict__.keys())
        for (_, analysis) in self.analyses:
            if self.chunks:
                self.run_params[analysis]['chunks'] = self.chunks

    def set_user_params(self):
        # update these run parameters is file is passed
        if self.run_params_fp and isfile(self.run_params_fp):
            print('* Updating with user-defined parameters', self.run_params_fp)
            run_params_update = read_yaml_file(self.run_params_fp)
            # update these run parameters is file is passed
            self.check_run_params(run_params_update)

    def get_filt_raref_suffix(self):
        """Create a suffix based on passed config to denote
         whether the run is for filtered and/or rarefied data.
        """
        self.filt_raref = ''
        if self.filter_fp:
            self.filt_raref += '_flt'
        if self.raref:
            self.filt_raref += '_rrf'

    def get_train_test_dict(self):
        self.train_test_dict = self.train_test
        if 'train' not in self.train_test:
            self.train_test_dict['train'] = 0.7
        elif float(self.train_test['train']) < 0:
            self.train_test_dict['train'] = 0.7
        elif float(self.train_test['train']) > 1:
            self.train_test_dict['train'] = 0.7

    @staticmethod
    def check_xpbs_install():
        """Try to get the install path of third party tool
        Xhpc (https://github.com/FranckLejzerowicz/Xhpc).
        If it exists, nothing happens and the code proceeds.
        Otherwise, the code ends and tells what to do.
        """
        ret_code, ret_path = subprocess.getstatusoutput('which Xhpc')
        if ret_code:
            print('Xhpc is not installed (and make sure '
                  'to edit its config.txt)\nExiting...')
            sys.exit(1)
        else:
            with open(ret_path) as f:
                for line in f:
                    break
            if line.startswith('$HOME'):
                print('Xhpc is installed but its config.txt '
                      'need editing!\nExiting...')
                sys.exit(1)

    def get_analyses(self):
        self.analyses = [('Import data table to Qiime2', 'import')]

        if self.filt3d:
            self.analyses.append(('Explore filtering (in 3D)', 'filter3d'))
        if self.filter:
            self.analyses.append(('Rare sequence filtering', 'filter'))
        if self.raref:
            self.analyses.append(('Rarefaction filtering', 'rarefy'))

        if self.qemistree and 'qemistree' not in self.skip:
            self.analyses.append(('Chemically-informed trees for metabolites',
                                  'qemistree'))

        if 'taxonomy' not in self.skip:
            self.analyses.append(('Taxonomic assignment', 'taxonomy'))
            if 'barplot' not in self.skip:
                self.analyses.append(('Draw barplots', 'barplot'))

        if 'wol' not in self.skip:
            self.analyses.append(('Shear Web-Of-Life tree', 'wol'))
        if 'sepp' not in self.skip:
            self.analyses.append(('Phylogenetic read placement (SEPP)', 'sepp'))

        if 'pies' not in self.skip:
            self.analyses.append(('Draw pie charts', 'pies'))

        if self.collapse and 'collapse' not in self.skip:
            self.analyses.append(('Collapse taxonomy', 'collapse'))

        if self.feature_subsets and 'feature_subsets' not in self.skip:
            self.analyses.append(('Subset feature groups', 'feature_subsets'))

        if 'alpha' not in self.skip:
            self.analyses.append(('Alpha diversity measures', 'alpha'))
            if 'alpha_merge' not in self.skip:
                self.analyses.append(('Merge alpha diversity', 'alpha_merge'))
            if 'alpha_rarefaction' not in self.skip:
                self.analyses.append(('Alpha diversity rarefaction',
                                      'alpha_rarefaction'))
            if 'alpha_correlations' not in self.skip:
                self.analyses.append(('Alpha diversity correlations',
                                      'alpha_correlations'))
            if 'alpha_group_significance' not in self.skip:
                if 'alpha_kw' not in self.skip:
                    self.analyses.append(('Kruskal-Wallis for alpha diversity',
                                          'alpha_group_significance'))
            if self.longi_column and 'volatility' not in self.skip:
                self.analyses.append(('Alpha diversity volatility',
                                      'volatility'))

        if self.phate and 'phate' not in self.skip:
            self.analyses.append(('Potential of Heat-diffusion for Affinity-'
                                  'based Trajectory Embedding (PHATE)',
                                  'phate'))

        if 'beta' not in self.skip:
            self.analyses.append(('Beta diversity measures', 'beta'))
            if 'deicode' not in self.skip:
                self.analyses.append(('Robust Aitchison PCA', 'deicode'))
            if 'pcoa' not in self.skip:
                self.analyses.append(('Principal coordinate analysis', 'pcoa'))
            if 'tsne' not in self.skip:
                self.analyses.append(
                    ('t-Distributed Stochastic Neighbor Embedding', 'tsne'))
            if 'umap' not in self.skip:
                self.analyses.append(
                    ('Uniform Manifold Approximation and Projection', 'umap'))
            if 'emperor' not in self.skip:
                self.analyses.append(('EMPeror plot', 'emperor'))
            if 'empress' not in self.skip:
                self.analyses.append(('EMPress plot', 'empress'))
            if self.biplot and 'biplot' not in self.skip:
                self.analyses.append(('Principal coordinate analysis (biplot)',
                                      'biplot'))
                if 'emperor_biplot' not in self.skip:
                    self.analyses.append(('EMPeror biplot', 'emperor_biplot'))
                if 'empress_biplot' not in self.skip:
                    self.analyses.append(('EMPress biplot', 'empress_biplot'))

            if self.tests and 'permanova' not in self.skip:
                self.analyses.append(
                    ('Permutational Multivariate Analysis of Variance (Qiime2)',
                     'permanova'))
            if self.adonis and 'adonis' not in self.skip:
                self.analyses.append(
                    ('Permutational Multivariate Analysis of Variance (R)',
                     'adonis'))

            if self.procrustes and 'procrustes' not in self.skip:
                self.analyses.append(('Procrustes analysis', 'procrustes'))
                self.analyses.append(('Procrustes analysis (R)',
                                      'procrustes_r'))
            if self.mantel and 'mantel' not in self.skip:
                self.analyses.append(('Mantel tests', 'mantel'))

            if self.nestedness and 'nestedness' not in self.skip:
                self.analyses.append(('Nestedness analysis', 'nestedness'))
                self.analyses.append(('Nestedness Over Decreasing Fill (NODF)',
                                      'nestedness_nodfs'))
                self.analyses.append(('Nestedness graphs', 'nestedness_graphs'))

            if self.dm_decay and 'dm_decay' not in self.skip:
                self.analyses.append(('Beta-distance decay', 'dm_decay'))
                self.analyses.append(('Beta-distance decay (plot)',
                                      'dm_decay_plot'))
            if self.geo_decay and 'geo_decay' not in self.skip:
                self.analyses.append(('Geographic-distance decay', 'geo_decay'))

        if self.sourcetracking and 'sourcetracking' not in self.skip:
            self.analyses.append(('Source-tracking', 'sourcetracking'))
        if self.doc and 'doc' not in self.skip:
            self.analyses.append(('Dissimilarity-Overlap Curves', 'doc'))

        if self.mmvec_pairs and 'mmvec' not in self.skip:
            self.analyses.append(('Microbe-metabolite co-occurrence estimation',
                                  'mmvec'))
        if self.diff_models and 'songbird' not in self.skip:
            self.analyses.append(('Log-fold change with multinomial regression',
                                  'songbird'))
            if 'qurro' not in self.skip:
                self.analyses.append(('Log-fold change visualization',
                                      'qurro'))
        elif self.mmvec_pairs and 'mmbird' not in self.skip:
            self.analyses.append(('Integrating interactions and differentials',
                                  'mmbird'))

