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
from os.path import abspath, dirname, exists, isfile, isdir
from microbiome_analyzer._io import read_yaml_file

RESOURCES = pkg_resources.resource_filename("microbiome_analyzer", "resources")


class AnalysesConfig(object):
    """
    Collect the data associated with each dataset passed but the user
    """
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.folder = ''
        self.prjct_nm = ''
        self.filt_raref = ''
        self.subsets = {}
        self.workflow = {}
        self.run_params = {}
        self.train_test_dict = {}

    def init(self):
        self.check_input_output()
        self.parse_yamls()
        self.get_prjct_anlss_nm()  # project name for jobs
        self.get_run_params()  # default (e.g. memory)
        self.get_filt_raref_suffix()  # job suffix (e.g. _flt)
        self.get_train_test_dict()
        if self.jobs:
            self.check_xpbs_install()  # Xpbs must be installed to prepare jobs
        self.get_workflow()

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

    def check_input_output(self):
        """Check that the main input folder exists
            and that it is indeed a folder, not a file.
        Returns its absolute path for unambiguous use
            in the rest of the code.
        """
        if not exists(self.data_folder):
            sys.exit('%s does exist\nExiting...' % self.data_folder)
        elif not isdir(self.data_folder):
            sys.exit('%s must be a folder\nExiting...' % self.data_folder)
        if not exists(self.metadata_folder):
            sys.exit('%s does exist\nExiting...' % self.metadata_folder)
        elif not isdir(self.metadata_folder):
            sys.exit('%s must be a folder\nExiting...' % self.metadata_folder)
        self.data_folder = abspath(self.data_folder)
        self.metadata_folder = abspath(self.metadata_folder)

        if '/' in self.output_folder:
            if not isdir(dirname(self.output_folder)):
                sys.exit('output folder "%s" must exist\nExiting...' %
                         self.output_folder)
        self.output_folder = abspath(self.output_folder)

    def get_prjct_anlss_nm(self):
        """
        Get a smaller name for printing in qstat / squeue.
        """
        alpha = 'aeiouy'
        self.prjct_nm = ''.join(
            x for x in self.project_name if x.lower() not in alpha)
        if self.prjct_nm == '':
            self.prjct_nm = self.project_name

    def check_env(self, conda_envs: set, env: str = None) -> None:
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
        if conda_env not in conda_envs:
            print('%s is not an existing conda environment' % env)
            sys.exit(1)

    def get_conda_envs(self) -> list:
        """Get the names of the conda environments.

        Returns
        -------
        conda_envs : set
            Conda environments.
        """
        conda_envs = []
        for conda_env in subprocess.getoutput('conda env list').split('\n'):
            if len(conda_env) and conda_env[0] != '#':
                conda_envs.append(conda_env.split()[0])
        # check that qiime2 environment exists
        self.check_env(set(conda_envs))
        return conda_envs

    def check_run_params(
            self,
            run_params_user: dict,
            conda_envs: set):
        """Verify that the keys/values passed for
        an analysis are valid.

        Parameters
        ----------
        run_params_user : dict
            Passed run parameters.
        conda_envs: set
            Conda environments.
        """
        mem_dim = ['kb', 'mb', 'gb']
        integer_values = ['time', 'n_nodes', 'n_procs', 'mem_num']
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
                    self.check_env(conda_envs, param_val)
                    self.run_params[analysis][param] = param_val

    def get_run_params(self):
        """Get the run parameters based on the default
        values, that are possibly updated by the user.
        """
        conda_envs = self.get_conda_envs()
        # read default run parameters (one set for each analysis)
        run_params_fp = '%s/run_params.yml' % RESOURCES
        self.run_params = read_yaml_file(run_params_fp)
        # update these run parameters is file is passed
        if self.run_params_fp and isfile(self.run_params_fp):
            run_params_update = read_yaml_file(self.run_params_fp)
            # update these run parameters is file is passed
            self.check_run_params(run_params_update, set(conda_envs))
        else:
            print('Using run parameters from', run_params_fp)

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

    def get_workflow(self):
        self.workflow = [('Import',),]
        if self.filter:
            self.workflow.append(('Filter',))
        if self.raref:
            self.workflow.append(('Rarefy',))
        if self.filt3d:
            self.workflow.append(('Filter 3D',))
        # if self.i_qemistree and 'qemistree' not in self.skip:
        #     workflow.append(('Qemistree',))
        if 'taxonomy' not in self.skip:
            self.workflow.append(('Taxonomy',))
            if 'barplot' not in self.skip:
                self.workflow.append(('edit taxonomy',))
            self.workflow.append(('Edit taxonomy',))
        if 'wol' not in self.skip:
            self.workflow.append(('Shear WOL tree',))
        if self.sepp_tree and 'sepp' not in self.skip:
            self.workflow.append(('SEPP reads placement',))
        if self.filt_only:
            self.workflow.append(('Delete non-filtered',))
        if 'do_pies' in self.skip:
            self.workflow.append(('Make taxonomy pie charts',))
        if self.collapse and 'collapse' not in self.skip:
            self.workflow.append(('Collapse taxonomy',))
        if 'alpha' not in self.skip:
            self.workflow.append(('Alpha diversity',))
            if 'merge_alpha' not in self.skip:
                self.workflow.append(('Merge alpha diversity',))
                if 'export_alpha' not in self.skip:
                    self.workflow.append(('Export alpha diversity to metadata',))
            if 'alpha_rarefaction' not in self.skip:
                self.workflow.append(('Alpha diversity rarefaction',))
            if 'alpha_correlations' not in self.skip:
                self.workflow.append(('Alpha diversity correlations',))
            if self.longi_column and 'volatility' not in self.skip:
                self.workflow.append(('Alpha diversity volatility',))
            if 'alpha_group_significance' not in self.skip:
                if 'alpha_kw' not in self.skip:
                    self.workflow.append(('Alpha-diversity Kruskal-Wallis tests',))
        if 'beta' not in self.skip:
            self.workflow.append(('Beta diversity',))
            if 'export_beta' not in self.skip:
                self.workflow.append(('Export beta diversity',))
            if 'pcoa' not in self.skip:
                self.workflow.append(('Principal coordinate analysis',))
                if 'emperor' not in self.skip:
                    self.workflow.append(('EMPeror plot',))
                if 'empress' not in self.skip:
                    self.workflow.append(('EMPress plot',))
            if self.biplot and 'biplot' not in self.skip:
                self.workflow.append(('Principal coordinate analysis (biplot)',))
                if 'emperor_biplot' not in self.skip:
                    self.workflow.append(('EMPeror biplot',))
                if 'empress_biplot' not in self.skip:
                    self.workflow.append(('EMPress biplot',))
            if 'deicode' not in self.skip:
                self.workflow.append(('DEICODE',))
            if self.tests and 'permanova' not in self.skip:
                self.workflow.append(('PERMANOVA',))
            if self.adonis and 'adonis' not in self.skip:
                self.workflow.append(('Adonis',))
            if self.procrustes and 'procrustes' not in self.skip:
                self.workflow.append(('Procrustes analysis',))
            if self.mantel and 'mantel' not in self.skip:
                self.workflow.append(('Mantel tests',))
            if self.nestedness and 'nestedness' not in self.skip:
                self.workflow.append(('Nestedness analysis',))
            if self.dm_decay and 'dm_decay' not in self.skip:
                self.workflow.append(('Beta-distance decay',))
            if self.geo_decay and 'geo_decay' not in self.skip:
                self.workflow.append(('Geo-distance decay',))
        if self.phate and 'phate' not in self.skip:
            self.workflow.append(('phate',))
        if 'doc' not in self.skip and self.doc:
            self.workflow.append(('Dissimilarity-Overlap Curves',))
        if self.mmvec_pairs and 'mmvec' not in self.skip:
            self.workflow.append(('MMVEC',))
            if 'mmbird' not in self.skip:
                self.workflow.append(('MMBird',))
        if self.diff_models and 'songbird' not in self.skip:
            self.workflow.append(('Songbird',))
