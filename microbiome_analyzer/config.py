# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import subprocess
import pkg_resources
from os.path import abspath, exists, isfile
from microbiome_analyzer._io import read_yaml_file

RESOURCES = pkg_resources.resource_filename(
    "microbiome_analyzer", "resources")


class AnalysesConfig(object):
    """
    Collect the data associated with each dataset passed but the user
    """
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.folder = self.check_input()
        self.subsets = {}
        self.parse_yamls()
        self.prjct_nm = self.get_prjct_anlss_nm()  # project name for jobs
        self.run_params = self.get_run_params()  # default (e.g. memory)
        self.filt_raref = self.get_filt_raref_suffix()  # job suffix (e.g. _flt)
        self.train_test_dict = self.get_train_test_dict()
        if self.jobs:
            self.check_xpbs_install()  # Xpbs must be installed to prepare jobs
        self.workflow = self.get_workflow()

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

    def check_input(self) -> str:
        """Check that the main input folder exists
            and that it is indeed a folder, not a file.
        Returns its absolute path for unambiguous use
            in the rest of the code.
        """
        if not exists(self.datasets_folder):
            print('%s does exist\nExiting...' % self.datasets_folder)
            sys.exit(1)
        folder = abspath(self.datasets_folder)
        if isfile(folder):
            print('%s is a file. Needs a folder as input\nExiting...' %
                  folder)
            sys.exit(1)
        return folder

    def get_prjct_anlss_nm(self) -> str:
        """
        Get a smaller name for printing in qstat / squeue.

        Returns
        -------
        prjct_nm : str
            Same name without the vows ("aeiouy").
        """
        alpha = 'aeiouy'
        prjct_nm = ''.join(
            x for x in self.project_name if x.lower() not in alpha)
        if prjct_nm == '':
            prjct_nm = self.project_name
        return prjct_nm

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

    def check_run_params(self, run_params: dict, run_params_user: dict,
                         conda_envs: set) -> dict:
        """Verify that the keys/values passed for
        an analysis are valid.

        Parameters
        ----------
        run_params : dict
            Default run parameters.
        run_params_user : dict
            Passed run parameters.
        conda_envs: set
            Conda environments.

        Returns
        -------
        run_params : dict
            Validated run parameters.
        """
        integer_values = ['time', 'n_nodes', 'n_procs', 'mem_num']
        mem_dim = ['kb', 'mb', 'gb']
        for analysis in list(run_params_user.keys()):
            for param in list(run_params_user[analysis].keys()):
                param_val = run_params_user[analysis][param]
                if param in integer_values:
                    if str(param_val).isdigit():
                        run_params[analysis][param] = param_val
                elif param == 'mem_dim':
                    if param_val in mem_dim:
                        run_params[analysis][param] = param_val
                elif param == 'env':
                    # check that qiime2 environment exists
                    self.check_env(conda_envs, param_val)
                    run_params[analysis][param] = param_val
        return run_params

    def get_run_params(self) -> dict:
        """Get the run parameters based on the default
        values, that are possibly updated by the user.

        Returns
        -------
        run_params: dict
            Valid run parameters, which for each analysis could be
                'time': an integer
                'n_nodes': an integer
                'n_procs': an integer
                'mem_num': an integer
                'mem_dim': either on of ['kb', 'mb', 'gb']
                'env': an existing conda environment
        """
        conda_envs = self.get_conda_envs()
        # read default run parameters (one set for each analysis)
        run_params_fp = '%s/run_params.yml' % RESOURCES
        run_params = read_yaml_file(run_params_fp)
        # update these run parameters is file is passed
        if self.run_params_fp and isfile(self.run_params_fp):
            run_params_update = read_yaml_file(self.run_params_fp)
            # update these run parameters is file is passed
            run_params = self.check_run_params(
                run_params, run_params_update, set(conda_envs))
        else:
            print('Using run parameters from', run_params_fp)
        return run_params

    def get_filt_raref_suffix(self) -> str:
        """Create a suffix based on passed config to denote
         whether the run is for filtered and/or rarefied data.

        Returns
        -------
        filt_raref : str
            Suffix for output scripts denoting whether the
            run is for filtered and/or rarefied data
        """
        filt_raref = ''
        if self.filter_fp:
            filt_raref += '_flt'
        if self.raref:
            filt_raref += '_rrf'
        return filt_raref

    def get_train_test_dict(self) -> dict:
        train_test_dict = read_yaml_file(self.train_test)
        if 'train' not in train_test_dict:
            train_test_dict['train'] = 0.7
        elif float(train_test_dict['train']) < 0:
            train_test_dict['train'] = 0.7
        elif float(train_test_dict['train']) > 1:
            train_test_dict['train'] = 0.7
        return train_test_dict

    @staticmethod
    def check_xpbs_install():
        """Try to get the install path of third party tool
        Xpbs (https://github.com/FranckLejzerowicz/Xpbs).
        If it exists, nothing happens and the code proceeds.
        Otherwise, the code ends and tells what to do.
        """
        ret_code, ret_path = subprocess.getstatusoutput('which Xpbs')
        if ret_code:
            print('Xpbs is not installed (and make sure '
                  'to edit its config.txt)\nExiting...')
            sys.exit(1)
        else:
            with open(ret_path) as f:
                for line in f:
                    break
            if line.startswith('$HOME'):
                print('Xpbs is installed but its config.txt '
                      'need editing!\nExiting...')
                sys.exit(1)

    def get_workflow(self):
        workflow = [('Import',),]
        if self.filter:
            workflow.append(('Filter',))
        if self.raref:
            workflow.append(('Rarefy',))
        if self.filt3d:
            workflow.append(('Filter 3D',))
            return workflow
        # if self.i_qemistree and 'qemistree' not in self.skip:
        #     workflow.append(('Qemistree',))
        if 'taxonomy' not in self.skip:
            workflow.append(('Taxonomy',))
            if 'barplot' not in self.skip:
                workflow.append(('edit taxonomy',))
            workflow.append(('Edit taxonomy',))
        if 'wol' not in self.skip:
            workflow.append(('Shear WOL tree',))
        if self.sepp_tree and 'sepp' not in self.skip:
            workflow.append(('SEPP reads placement',))
        if self.filt_only:
            workflow.append(('Delete non-filtered',))
        if 'do_pies' in self.skip:
            workflow.append(('Make taxonomy pie charts',))
        if self.collapse and 'collapse' not in self.skip:
            workflow.append(('Collapse taxonomy',))
        if 'alpha' not in self.skip:
            workflow.append(('Alpha diversity',))
            if 'merge_alpha' not in self.skip:
                workflow.append(('Merge alpha diversity',))
                if 'export_alpha' not in self.skip:
                    workflow.append(('Export alpha diversity to metadata',))
            if 'alpha_rarefaction' not in self.skip:
                workflow.append(('Alpha diversity rarefaction',))
            if 'alpha_correlations' not in self.skip:
                workflow.append(('Alpha diversity correlations',))
            if self.longi_column and 'volatility' not in self.skip:
                workflow.append(('Alpha diversity volatility',))
            if 'alpha_group_significance' not in self.skip:
                if 'alpha_kw' not in self.skip:
                    workflow.append(('Alpha-diversity Kruskal-Wallis tests',))
        if 'beta' not in self.skip:
            workflow.append(('Beta diversity',))
            if 'export_beta' not in self.skip:
                workflow.append(('Export beta diversity',))
            if 'pcoa' not in self.skip:
                workflow.append(('Principal coordinate analysis',))
                if 'emperor' not in self.skip:
                    workflow.append(('EMPeror plot',))
                if 'empress' not in self.skip:
                    workflow.append(('EMPress plot',))
            if self.biplot and 'biplot' not in self.skip:
                workflow.append(('Principal coordinate analysis (biplot)',))
                if 'emperor_biplot' not in self.skip:
                    workflow.append(('EMPeror biplot',))
                if 'empress_biplot' not in self.skip:
                    workflow.append(('EMPress biplot',))
            if 'deicode' not in self.skip:
                workflow.append(('DEICODE',))
            if self.tests and 'permanova' not in self.skip:
                workflow.append(('PERMANOVA',))
            if self.adonis and 'adonis' not in self.skip:
                workflow.append(('Adonis',))
            if self.procrustes and 'procrustes' not in self.skip:
                workflow.append(('Procrustes analysis',))
            if self.mantel and 'mantel' not in self.skip:
                workflow.append(('Mantel tests',))
            if self.nestedness and 'nestedness' not in self.skip:
                workflow.append(('Nestedness analysis',))
            if self.dm_decay and 'dm_decay' not in self.skip:
                workflow.append(('Beta-distance decay',))
            if self.geo_decay and 'geo_decay' not in self.skip:
                workflow.append(('Geo-distance decay',))
        if self.phate and 'phate' not in self.skip:
            workflow.append(('phate',))
        if 'doc' not in self.skip and self.doc:
            workflow.append(('Dissimilarity-Overlap Curves',))
        if self.mmvec_pairs and 'mmvec' not in self.skip:
            workflow.append(('MMVEC',))
            if 'mmbird' not in self.skip:
                workflow.append(('MMBird',))
        if self.diff_models and 'songbird' not in self.skip:
            workflow.append(('Songbird',))
        return workflow
