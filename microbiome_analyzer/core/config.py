# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import time
import yaml
import subprocess

import pandas as pd
import pkg_resources
import numpy as np
from datetime import datetime as dt
from os.path import abspath, dirname, exists, isfile, isdir
from microbiome_analyzer._inputs import read_yaml_file, read_meta_pd
from microbiome_analyzer._formats import check_format

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
        self.subsets['ALL'] = {}
        for arg in list(self.__dict__):
            if arg.endswith('_fp'):
                yaml = read_yaml_file(self.__dict__[arg])
                setattr(self, arg[:-3], yaml)
                if arg == 'sample_subsets_fp':
                    check_format(yaml, 'sample_subsets')
                    for name, var_vals in yaml.items():
                        self.subsets[name] = {}
                        for var, vals in var_vals.items():
                            self.subsets[name][var] = tuple(vals)

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
        self.check_env()

    def check_run_params(self, run_params_user: dict):
        """Verify that the keys/values passed for an analysis are valid.

        Parameters
        ----------
        run_params_user : dict
            Set of steps' parameters defined by the user
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
        """For every analysis, set the default scratch space to use to is
        None, and here it is update to the scatch space tols to be used at
        command line, in the following order of priority:
            - localscracth (because it is the most efficient)
            - scratch (created and erased for each job)
            - userscratch (not erased at the end of the job)
        In practice, each analysis' "scratch" key is given as value to name
        of the scratch location, except for localscratch, which will be a
        number corresponding to the amount of storage to use on the node.
        """
        # for each analysis
        message = '* Set default scratch location to use to:'
        for (_, analysis) in self.analyses:
            # prioritize is `--localscratch #` is given
            if self.localscratch:
                self.run_params[analysis]['scratch'] = self.localscratch
                print('%s localscratch (with %sG storage)' % message)
            # or use `--scratch` boolean
            elif self.scratch:
                self.run_params[analysis]['scratch'] = 'scratch'
                print('%s scratch' % message)
            # or use `--userscratch` boolean
            elif self.userscratch:
                self.run_params[analysis]['scratch'] = 'userscratch'
                print('%s userscratch' % message)

    def set_global_chunks(self):
        """For every analysis, collect the number of chunks (will be the
        number of jobs) to make for the command lines of this analyis,
        if more than zero chunk is given in command or in the user-params.
        """
        # for each analysis
        for (_, analysis) in self.analyses:
            # if at least one chunk is told to be used
            if self.chunks:
                # then collect this number of chunks/jobs for the key "chunks"
                print('* Set default number of jobs: %s per step' % self.chunks)
                self.run_params[analysis]['chunks'] = self.chunks

    def set_user_params(self):
        """If the user defined alternative paramters for one of more analysis
        steps, the corresponding yaml file contaiing these parameters is read
        and the parameters are used to update the default parameters
        previously loaded.
        """
        # update these run parameters is file is passed
        if self.run_params_fp and isfile(self.run_params_fp):
            print('* Updating parameters form user config:', self.run_params_fp)
            run_params_update = read_yaml_file(self.run_params_fp)
            # update these run parameters is file is passed
            self.check_run_params(run_params_update)

    def get_filt_raref_suffix(self):
        """Create a suffix based on passed config to denote whether the run
        is for filtered and/or rarefied data.
        """
        self.filt_raref = ''
        if self.filter_fp:
            self.filt_raref += '_flt'
        if self.rarefy:
            self.filt_raref += '_rrf'

    def get_train_test_dict(self):
        """Set the percent of samples to use for training in ML-based steps."""
        # default the train-test rates to the train_test (None by default)
        self.train_test_dict = self.train_test
        # set to 70% of sample for training if
        if 'train' not in self.train_test:
            # the train-test parameters do not specify a % of training samples
            self.train_test_dict['train'] = 0.7
        elif float(self.train_test['train']) < 0:
            # the value is below 0
            self.train_test_dict['train'] = 0.7
        elif float(self.train_test['train']) > 1:
            # the value is higher that 1
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
        if self.rarefy:
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

        if self.time_subject and 'volatility' not in self.skip:
            self.analyses.append(('Volatility analysis', 'volatility'))
        if self.time_subject and 'ctf' not in self.skip:
            self.analyses.append(('Compositional Tensor Factorisation analysis',
                                  'ctf'))
        if 'rpca' not in self.skip:
            self.analyses.append(('Robust-Principal Componenz Analysis',
                                  'rpca'))
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

            if self.permanova and 'permanova' not in self.skip:
                self.analyses.append(
                    ('Permutational Multivariate Analysis of Variance (Qiime2)',
                     'permanova'))
            if self.adonis and 'adonis' not in self.skip:
                self.analyses.append(
                    ('Permutational Multivariate Analysis of Variance (R)',
                     'adonis'))
            if self.mantel and 'mantel' not in self.skip:
                self.analyses.append(('Mantel tests', 'mantel'))
            if self.procrustes and 'procrustes' not in self.skip:
                self.analyses.append(('Procrustes analysis', 'procrustes'))
                self.analyses.append(('Procrustes analysis: R','procrustes_r'))
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
        if self.diff_abund and 'songbird' not in self.skip:
            self.analyses.append(('Log-fold change with multinomial regression',
                                  'songbird'))
            if 'qurro' not in self.skip:
                self.analyses.append(('Log-fold change visualization',
                                      'qurro'))
        elif self.mmvec_pairs and 'mmbird' not in self.skip:
            self.analyses.append(('Integrating interactions and differentials',
                                  'mmbird'))


class PrepConfig(object):
    """
    Create config files for various data analyses
    """
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.empty = True
        self.meta = {}
        self.vars = {}
        self.sams = set()
        self.configs = {}
        self.config = []
        self.config_fp = ''
        self.analysis = ''
        self.time = dt.now().strftime("%d-%m-%Y_%H") + 'h'
        self.out = '%s/configs/%s' % (self.analysis_folder, self.time)
        self.check_io()
        self.get_vars()

    def setup(self):
        self.prep = '[Config: %s]' % self.analysis
        self.config_fp = '%s/%s.yml' % (self.out, self.analysis)
        sep = '\n%s\n' % ('=' * (len(self.prep) + 10))
        print('%s     %s%sPreparing template: %s' % (
            sep, self.prep, sep, self.config_fp))
        print('<Write QUIT, STOP or EXIT to avoid the config>')
        print('<Press ENTER to skip params>\n')

    def check_io(self):
        if not isdir(self.analysis_folder):
            sys.exit('Input folder not found: "%s"' % self.analysis_folder)
        if not isdir('%s/metadata' % self.analysis_folder):
            sys.exit('"%s/metadata" not found' % self.analysis_folder)
        for dat in self.datasets:
            meta = '%s/metadata/%s/metadata.tsv' % (self.analysis_folder, dat)
            if not isfile(meta):
                sys.exit('No metadata file "%s"' % meta)
            self.meta[dat] = read_meta_pd(meta)
        if not isdir(self.out):
            os.makedirs(self.out)

    def get_vars(self):
        vars = self.collect_vars()
        for var, vals in vars.items():
            for rm in [np.nan, 'not applicable',' not collected',
                       'not provided', 'restricted access']:
                if rm in vals:
                    vals.remove(rm)
            integer = True
            for val in vals:
                if not str(val).isdigit():
                    integer = False
                    try:
                        float(val)
                    except ValueError:
                        self.vars[var] = (str, set(sorted(vals)))
                        break
            else:
                if integer:
                    self.vars[var] = (int, set(sorted(vals)))
                else:
                    self.vars[var] = (float, set(sorted(vals)))

    def collect_vars(self):
        vars = {}
        for dat, meta in self.meta.items():
            self.sams.update(set(meta['sample_name']))
            for var, vals in meta.apply(lambda x: x.unique).to_dict().items():
                if str(vals().dtype).startswith('int'):
                    self.vars[var] = (int, set(sorted(vals())))
                elif str(vals().dtype).startswith('float'):
                    self.vars[var] = (float, set(sorted(vals())))
                elif str(vals().dtype) == 'bool':
                    self.vars[var] = (bool, set(sorted(vals())))
                else:
                    vars.setdefault(var, set()).update(set(vals()))
        return vars

    def check_user_input(self, inp, digit=False, num=False) -> bool:
        if inp in ['QUIT', 'STOP', 'EXIT']:
            print('\n\t <<< Skipping config >>>')
            return True
        if not inp:
            sys.exit('%s\tError: empty input' % self.prep)
        elif ' ' in inp:
            sys.exit('%s\tError: space in input' % self.prep)
        elif not digit and inp[0].isdigit():
            sys.exit('%s\tError: input starts with number' % self.prep)
        elif digit:
            if num:
                try:
                    float(inp)
                except ValueError:
                    if not str(inp).isdigit():
                        sys.exit('%s\tError: input not int/float' % self.prep)
            elif not str(inp).isdigit():
                sys.exit('%s\tError: input is not a digit' % self.prep)
        self.empty = False
        return False

    def prep_filtering(self):
        for dat in self.datasets:
            if self.all_datasets:
                dat = 'All datasets'
                sams = self.sams
            else:
                sams = set(self.meta[dat].sample_name)

            writer = '%s:\n' % dat
            names = set()
            mess = '%s(%s) Samples to remove: ' % (self.prep, dat)
            inp = input(mess)
            while inp:
                if self.check_user_input(inp):
                    return
                if inp not in sams:
                    print('%s\tWarning: No sample "%s"' % (self.prep, inp))
                else:
                    names.add(inp)
                inp = input(mess)
            if names:
                writer += '  names:\n'
                for sam in names:
                    writer += '  - "%s"\n' % sam
            inp = input('%s(%s) Feature threshold: ' % (self.prep, dat))
            if inp:
                if self.check_user_input(inp, True, True):
                    return
                writer += '  features: %s\n' % inp
            inp = input('%s(%s) Sample threshold: ' % (self.prep, dat))
            if inp:
                if self.check_user_input(inp, True, True):
                    return
                writer += '  samples: %s\n' % inp
            if self.all_datasets:
                self.write_all_datasets(writer)
                break
            else:
                self.config.append(writer)

    def prep_rarefactions(self):
        for dat in self.datasets:
            if self.all_datasets:
                dat = 'All datasets'
            writer = '%s:\n' % dat
            mess = '%s(%s) Enter rarefaction depth: ' % (self.prep, dat)
            inp = input(mess)
            while inp:
                if inp == 'min':
                    writer += '  - "min"\n'
                else:
                    if self.check_user_input(inp, True):
                        return
                    writer += '  - "%s"\n' % inp
                inp = input(mess)
            if self.all_datasets:
                self.write_all_datasets(write)
                break
            else:
                self.config.append(writer)

    def prep_feature_subsets(self):
        for dat in self.datasets:
            writer = '%s:\n' % dat
            mess = '%s(%s) Enter subset name: ' % (self.prep, dat)
            inp = input(mess)
            while inp:
                if self.check_user_input(inp):
                    return
                writer += '  %s:\n' % inp
                mess2 = '%s(%s) Feature name/regex: ' % (self.prep, dat)
                inp2 = input(mess2)
                while input(inp2):
                    if self.check_user_input(inp2):
                        return
                    writer += '    - %s\n' % inp2
                    inp2 = input(mess2)
                inp = input(mess)

            if self.all_datasets:
                self.write_all_datasets(write)
                break
            else:
                self.config.append(writer)

    def get_str_bool_vals(self, inp, vals):
        groups = []
        valx = {str(x): y for x, y in enumerate(vals)}
        p = '\n'.join(['%s\t: %s' % (x, y) for x, y in valx.items()])
        mess = '%s\tFactors to group for "%s":' % (self.prep, inp)
        i = input('%s\n%s\nspace-separated indices: ' % (mess, p))
        j = set(i.strip().split())
        off = j.difference(set(valx))
        if off:
            sys.exit('%s\tWrong index selection: %s' % (self.prep, off))
        groups.extend([valx[k] if x else valx[k] for x, k in enumerate(j)])
        self.empty = False
        return groups

    def get_int_float_vals(self, inp, vals, typ):
        groups = []
        desc = str(pd.Series(list(vals), dtype=typ).describe()).split('dtyp')[0]
        mess = '%s\tBounding values to make group for "%s":' % (self.prep, inp)
        i = input('%s\n%s\n("<" or ">" + value): ' % (mess, desc))
        if i[0] not in '<>':
            sys.exit('%s\tInput must start with ">" or "<"' % self.prep)
        if ' ' in i:
            sys.exit('%s\tInput must not contain any space' % self.prep)
        val = i[1:].strip()
        try:
            float(val)
        except ValueError:
            if not str(val).isdigit():
                sys.exit('%s\tInput must be numeric' % self.prep)
        groups.append(i)
        self.empty = False
        return groups

    def get_subsets(self, dat):
        mess = '%s(%s) Subset name (free text / no space): ' % (self.prep, dat)
        name = input(mess)
        subsets = {}
        while name:
            if self.check_user_input(name):
                return
            var = '%s(%s) Enter metadata variable name: ' % (self.prep, dat)
            inp = input(var)
            subset = {}
            while inp:
                if self.check_user_input(inp):
                    return
                if inp not in self.vars:
                    print('Variable "%s" not in metadata files' % inp)
                else:
                    typ, vals = self.vars[inp]
                    if typ in [str, bool]:
                        groups = self.get_str_bool_vals(inp, vals)
                    else:
                        groups = self.get_int_float_vals(inp, vals, typ)
                    if groups:
                        subset[inp] = groups
                inp = input(var)
            subsets[name] = subset
            name = input(mess)
        return subsets

    def prep_sample_subsets(self):
        groups = {}
        dat = 'All datasets'
        subsets = self.get_subsets(dat)
        if subsets:
            self.config.append(yaml.dump(subsets))

    def write_all_datasets(self, writer, head=''):
        if head:
            writer = head + writer
            # self.config.append(head)
        for dat in self.datasets:
            self.config.append(writer.replace('All datasets', dat))

    def prep_time_subject(self):
        for ddx, dat in enumerate(self.datasets):
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)

            time = ''
            i = input('%s> Variable for TIME (Press ENTER to stop):' % h)
            while i:
                if self.check_user_input(i):
                    return
                if i in self.vars:
                    time = i
                    break
                print('No variable "%s" in metadata (ignored)' % i)
                i = input('%s> Variable for TIME (Press ENTER to stop):' % h)

            subject = ''
            i = input('%s> Variable for SUBJECT (Press ENTER to stop):' % h)
            while i:
                if self.check_user_input(i):
                    return
                if i in self.vars:
                    subject = i
                    break
                print('No variable "%s" in metadata (ignored)' % i)
                i = input('%s> Variable for SUBJECT (Press ENTER to stop):' % h)

            if time and subject:
                head = '%s:\n' % dat
                writer = '  time: %s\n' % time
                writer += '  subject: %s\n' % subject
                if self.all_datasets:
                    self.write_all_datasets(writer, head)
                    break
                else:
                    self.config.extend([head, writer])

    def prep_permanova(self):
        for ddx, dat in enumerate(self.datasets):
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)
            i = input('%s> Variable name (Press ENTER to stop): ' % h)
            written = set()
            while i:
                if self.check_user_input(i):
                    return
                if i not in self.vars:
                    print('No variable "%s" in metadata (ignored)' % i)
                else:
                    if i in writer:
                        print('("%s" already in the selection)' % i)
                    written.add(i)
                i = input('%s> Variable name (Press ENTER to stop): ' % h)
            if written:
                head = '%s:\n' % dat
                writer = '\n'.join(['- "%s"' % x for x in written])
                if self.all_datasets:
                    self.write_all_datasets(writer, head)
                    break
                else:
                    self.config.extend([head, writer])

    def prep_adonis(self):
        stratas = {}
        for ddx, dat in enumerate(self.datasets):
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)
            writer = '  %s:\n' % dat
            print('%s* Defining models:' % h)
            i = input('%sModel name: ' % h)
            while i:
                if self.check_user_input(i):
                    return
                form = input('%s> Formula: ' % h)
                if form:
                    writer += '    %s: "%s"\n' % (i, form)
                else:
                    write += '    %s: "FORMULA_NEEDED_HERE"\n' % i
                if dat not in stratas:
                    stratas[dat] = {}
                stratas[dat][i] = set()
                i = input('%sModel name: ' % h)

            if self.all_datasets:
                self.write_all_datasets(writer, 'models:\n')
                break
            else:
                if not ddx:
                    self.config.append('models:\n' + writer)
                self.config.append(writer)

        for dat in self.datasets:
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)
            mess = '%s\n%s* Define models stratas: '% (h, h)
            print(mess)
            for model in stratas[dat]:
                i = input('%sModel "%s":\n%s> Metadata variable: ' % (
                    h, model, h))
                while i:
                    if self.check_user_input(i):
                        return
                    if i not in self.vars:
                        print('No variable "%s" in metadata files' % i)
                    else:
                        stratas[dat].setdefault(model, set()).add(i)
                    i = input('%s> Metadata variable: ' % h)
            if self.all_datasets:
                break

        global_strata = ''
        i = input('%s Global strata? (all datasets/models)\n%s > Metadata '
                  'variable: ' % (self.prep, self.prep))
        while i:
            if self.check_user_input(i):
                return
            if i not in self.vars:
                print('No variable "%s" in metadata files (ignored)' % i)
            else:
                global_strata = i
                break
            i = input('%s> Metadata variable: ' % h)

        if stratas:
            for ddx, dat in enumerate(self.datasets):
                if self.all_datasets:
                    dat = 'All datasets'
                if stratas.get(dat):
                    writer = '  %s:\n' % dat
                    for model, vs in stratas[dat].items():
                        for v in vs:
                            writer += '    %s: "%s"\n' % (model, v)
                    if self.all_datasets:
                        self.write_all_datasets(writer, 'strata:\n')
                        if global_strata:
                            self.config.append('  global: "%s"\n' % global_strata)
                        break
                    else:
                        if not ddx:
                            self.config.append('strata:\n')
                        self.config.append(writer)

        if not self.all_datasets and global_strata:
            self.config.append('  global: "%s"\n' % global_strata)

    def prep_nestedness(self):
        h = '%s(All samples) ' % self.prep
        binary = input('%s* Path to "Nestedness" java binary: ' % h)
        if binary in ['QUIT', 'STOP', 'EXIT']:
            print('\t <<< Skipping config >>>')
            return True
        if binary:
            self.config.append('soft: "%s"\n' % binary)
        else:
            self.config.append('soft: "BINARY_NEEDED_HERE"\n')

        subsets = self.get_subsets('All samples')
        if subsets:
            self.config.append('subsets:\n')
            for inp, groups in subsets.items():
                self.config.append('  %s:\n  %s\n' % (inp, '\n  '.join(groups)))

        print('%s* Samples coloring: ' % h)
        mess = '%s  Metadata variable (categorical): ' % h
        inp = input(mess)
        cols = set()
        while inp:
            if self.check_user_input(inp):
                return
            if inp not in self.vars:
                print('No variable "%s" in metadata files' % inp)
            else:
                typ, vals = self.vars[inp]
                if typ != str:
                    print('No variable "%s" not categorical' % inp)
                else:
                    cols.add(inp)
            inp = input(mess)
        if cols:
            self.config.append(
                'sample_colors:\n- "%s"\n' % ('"\n- "'.join(sorted(cols))))

        print('%s* Features coloring: ' % h)
        mess = '%s  Taxonomic level (number from 1 to <max_level>): ' % h
        inp = input(mess)
        cols = set()
        while inp:
            if inp not in map(str, range(0, 50)):
                print('Not a numeric: "%s" (using default)' % inp)
            else:
                cols.add(str(inp))
            inp = input(mess)
        if cols:
            self.empty = False
            self.config.append(
                'feature_colors:\n- %s\n' % ('\n- '.join(sorted(cols))))

        print('%s* NODF variables: ' % h)
        mess = '%s  Metadata variable (categorical): ' % h
        inp = input(mess)
        cols = set()
        while inp:
            if self.check_user_input(inp):
                return
            if inp not in self.vars:
                print('No variable "%s" in metadata files' % inp)
            else:
                typ, vals = self.vars[inp]
                if typ != str:
                    print('No variable "%s" not categorical' % inp)
                else:
                    cols.add(inp)
            inp = input(mess)
        if cols:
            self.config.append('nodfs:\n- %s\n' % ('\n- '.join(sorted(cols))))

        print('%s* Comparison modes: ' % h)
        modes = set()
        for m in ["overall", "withineachtype", "betweeneachpairoftypes"]:
            inp = input('%s  Use mode "%s": <[y]/n>' % (h, m))
            if not inp or inp == 'y':
                modes.add(m)
            elif inp != 'n':
                print('Must be "y" or "n"')
        if modes:
            self.config.append(
                'modes:\n- "%s"\n' % ('"\n- "'.join(sorted(modes))))

        print('%s* Null models: ' % h)
        nulls = set()
        for n in ["equiprobablefixed", "fixedequiprobable"]:
            inp = input('%s  Use null "%s": <[y]/n>' % (h, n))
            if not inp or inp == 'y':
                nulls.add(n)
            elif inp != 'n':
                print('Must be "y" or "n"')
        if nulls:
            self.config.append(
                'nulls:\n- "%s"\n' % ('"\n- "'.join(sorted(nulls))))

        iterations = input('%s* Number of iterations (default=20): ' % h)
        if not iterations.isdigit():
            print('Not a numeric: "%s" (using default)' % inp)
        elif iterations:
            self.config.append('params:\n  iterations: %s\n' % iterations)
        else:
            self.config.append('params:\n  iterations: 20\n')

    def get_pairs(self):
        sams_idx = {str(x): y for x, y in enumerate(self.datasets)}
        sams = '\n'.join(['%s\t: %s' % s for s in sams_idx.items()])

        pairs = {}
        h = '%s(sample pairs) ' % self.prep
        mess = '%sSamples to pair:\n%s\n' % (h, sams)
        i = input('%s%sTwo space-separated indices: ' % (mess, h))
        if i in ['QUIT', 'STOP', 'EXIT']:
            print('\t <<< Skipping config >>>')
            return pairs

        while i:
            j = set(i.strip().split())
            if len(j) != 2:
                sys.exit('Please enter just two indices')
            off = j.difference(set(sams_idx))
            if off:
                sys.exit('Wrong index selection: %s' % off)
            name = input('%sNickname for this pair: ' % h)
            while 1:
                if name:
                    if self.check_user_input(name):
                        return
                    pairs[name] = [sams_idx[x] for x in j]
                    self.empty = False
                    break
                name = input('%sNickname for this pair: ' % h)
            i = input('%s%sTwo space-separated indices: ' % (mess, h))
        return pairs

    def prep_pairs(self):
        pairs = self.get_pairs()
        if pairs:
            self.config.append('pairs:\n')
            for pair, sams in pairs.items():
                self.config.append('  %s:\n' % pair)
                for sam in sams:
                    self.config.append('  - "%s"\n' % sam)

    def prep_dm_decay(self):
        h = '%s(All samples) ' % self.prep
        modes = {}
        print('%s* Comparison modes: ' % h)
        for m in ['between', 'within']:
            inp = input('%s  Use mode "%s": <[y]/n>' % (h, m))
            if not inp or inp == 'y':
                mess = '%s  Metadata variable (categorical): ' % h
                v = input(mess)
                cols = set()
                while v:
                    if self.check_user_input(v):
                        return
                    if v not in self.vars:
                        print('No variable "%s" in metadata files' % v)
                    else:
                        typ, vals = self.vars[v]
                        if typ != str:
                            print('No variable "%s" not categorical' % v)
                        else:
                            modes.setdefault(m, set()).add(v)
                    v = input(mess)
            elif inp != 'n':
                print('Must be "y" or "n"')

        if modes:
            self.config.append('modes:\n')
            print(modes)
            for m, vs in modes.items():
                self.config.append(
                    '  %s:\n  - "%s"\n' % (m, '"\n  - "'.join(sorted(vs))))

        params = {}
        for p, v in [('iteration', '10'), ('step', '10')]:
            inp = input('%s* Number of %s (default=%s): ' % (h, p, v))
            params[p] = v
            if not inp.isdigit():
                print('Not a numeric: "%s" (using default)' % inp)
                params[p] = v
            elif inp:
                params[p] = inp
        if params:
            self.config.append('params:\n')
            for p, v in params.items():
                self.config.append('  %s: %s\n' % (p, v))

    def prep_geo_decay(self):
        print('<<< NOT AVAILABLE YET >>>')

    def prep_collapse(self):
        for dat in self.datasets:
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)
            collapses = {}
            inp = input('%sEnter a collapse name: ' % h)
            while inp:
                if inp:
                    if self.check_user_input(inp):
                        return
                    name = input('%s> Taxon identifier (or index): ' % h)
                    while 1:
                        if name:
                            collapses[inp] = name
                            break
                        name = input('%sNickname for this pair: ' % h)
                inp = input('%sEnter a collapse name: ' % h)
            if collapses:
                writer = '%s:\n' % dat
                for i, j in collapses.items():
                    writer += '  %s: "%s"\n' % (i, j)

                if self.all_datasets:
                    self.write_all_datasets(writer)
                    break
                else:
                    self.config.append(writer)
            elif self.all_datasets:
                break

    def prep_train_test(self):
        trains = {}
        for dat in self.datasets:
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)
            print('%s* Set training/testing splits:' % h)
            i = input('%s  Enter name of a split: ' % h)
            while i:
                if self.check_user_input(i):
                    return
                inp = input('%s  Metadata variable: ' % h)
                while inp:
                    if self.check_user_input(inp):
                        return
                    if inp not in self.vars:
                        print('No variable "%s" in metadata files' % inp)
                    else:
                        if dat not in trains:
                            trains[dat] = {}
                        trains[dat].setdefault(i, set()).add(inp)
                    inp = input('%s  Metadata variable: ' % h)
                i = input('%s  Enter name of a split: ' % h)
            if self.all_datasets:
                break

        if trains:
            for ddx, dat in enumerate(trains):
                writer = '  %s:\n' % dat
                for k, vs in trains[dat].items():
                    writer += '    %s:\n' % k
                    for v in vs:
                        writer += '    - "%s"\n' % v
                if self.all_datasets:
                    self.write_all_datasets(writer, 'datasets:\n')
                    break
                else:
                    if not ddx:
                        self.config.append('datasets:\n')
                    self.config.append(writer)

        inp = input('%s* Fraction for training (default=0.7): ' % h)
        try:
            float(inp)
            self.config.append('train: %s\n' % inp)
        except ValueError:
            print('Training input must be a float (using 0.7)')
            self.config.append('train: 0.7\n')

    def prep_phate(self):
        h = '%s(All samples) ' % self.prep
        print('%s* Labels for coloring: ' % h)
        inp = input('%s  - Metadata variable: ' % h)
        cols = set()
        while inp:
            if self.check_user_input(inp):
                return
            if inp not in self.vars:
                print('No variable "%s" in metadata files' % inp)
            else:
                cols.add(inp)
            inp = input('%s  - Metadata variable: ' % h)
        if cols:
            self.config.append(
                'labels:\n- "%s"\n' % ('"\n- "'.join(sorted(cols))))

        print('%s* Setup PHATE key parameters: ' % h)
        o.write('params:\n')
        for p, (s, e, t) in {'t': [2, 10, 5],
                             'decay': [1, 10, 2],
                             'knn': [1, 20, 5]}.items():
            self.config.append('  %s:\n' % p)
            d = '%s (start=%s; end=%s; step=%s)' % (
                str(list(range(s, (e + 1), t))), s, e, t)
            print('%s  - Parameters for "%s" (defaults: %s): ' % (h, p, d))
            for step, val in [('start', s), ('end', e), ('step', t)]:
                inp = input('%s    * %s [default=%s]: ' % (h, step, val))
                if not inp.isdigit():
                    print('Not a numeric: "%s" (using default)' % inp)
                    self.config.append('  - %s\n' % val)
                elif inp:
                    self.config.append('  - %s\n' % inp)

    def prep_sourcetracking(self):
        print('<<< NOT AVAILABLE YET >>>')

    def prep_doc(self):
        h = '%s(All samples) ' % self.prep
        print('%s* Setup DOC parameters: ' % h)
        params = []
        for p, v in {'r': 2, 'subr': 50, 'mov_avg': 5, 'degree': 1,
                           'iterations': 2, 'nulls': 1, 'non_zero': 1,
                           'null': 1}.items():
            mess = '%s  - Param for "%s" (INT; default=%s): ' % (h, p, v)
            inp = input(mess)
            if inp in ['QUIT', 'STOP', 'EXIT']:
                print('\t <<< Skipping config >>>')
                return True
            if not inp.isdigit():
                if inp:
                    print('Not an acceptable numeric: %s' % inp)
                params.append('  %s: %s\n' % (p, v))
            elif inp:
                params.append('  %s: %s\n' % (p, inp))
            self.empty = False
        if params:
            self.config.append('params:\n')
            self.config.extent(params)

        for p, v in {'span': 0.2}.items():
            mess = '%s  - Param for "%s" (FLOAT; default=%s): ' % (h, p, v)
            inp = input(mess)
            if inp in ['QUIT', 'STOP', 'EXIT']:
                print('\t <<< Skipping config >>>')
                return True
            if not inp.isdigit():
                if inp:
                    print('Not an acceptable numeric: %s' % inp)
                self.config.append('  %s: %s\n' % (p, v))
            elif inp:
                self.config.append('  %s: %s\n' % (p, inp))
            self.empty = False

        for p, v in {'family': 'symmetric', 'use_mp': 'false',
                     'surface': 'interpolate'}.items():
            mess = '%s  - Param for "%s" (TEXT default=%s): ' % (h, p, v)
            inp = input(mess)
            if inp in ['QUIT', 'STOP', 'EXIT']:
                print('\t <<< Skipping config >>>')
                return True
            if inp:
                self.config.append('  %s: %s\n' % (p, inp))
            else:
                self.config.append('  %s: %s\n' % (p, v))
            self.empty = False

        print('%s  - Param for "ci" (BOUNDS): ')
        self.config.append('  ci:\n')
        for n, v in [('low', '0.025'), ('mid', '0.5'), ('top', '0.975')]:
            inp = input('%s    * %s [default=%s]: ' % (h, n, v))
            if inp in ['QUIT', 'STOP', 'EXIT']:
                print('\t <<< Skipping config >>>')
                return True
            try:
                float(inp)
                self.config.append('  - "%s"\n' % inp)
            except ValueError:
                print('Not a float: "%s"' % inp)
                self.config.append('  - "%s"\n' % v)
            self.empty = False

        inp = input('%s* Use also the R version? <[y]/n>: ' % h)
        if inp in ['QUIT', 'STOP', 'EXIT']:
            print('\t <<< Skipping config >>>')
            return True
        if not inp or inp == 'y':
            self.config.append('do_r: 1\n')
        else:
            self.config.append('do_r: 0\n')
        self.empty = False

    def prep_xmmvec(self):
        print('<<< NOT AVAILABLE YET >>>')

    def prep_mmvec_highlights(self):
        print('<<< NOT AVAILABLE YET >>>')

    def prep_mmvec_pairs(self):
        print('<<< NOT AVAILABLE YET >>>')

    def prep_diff_abund(self):
        baselines = {}
        for ddx, dat in enumerate(self.datasets):
            if self.all_datasets:
                dat = 'All datasets'
            h = '%s(%s) ' % (self.prep, dat)
            writer = '  %s:\n' % dat
            print('%s* Defining models:' % h)
            i = input('%sModel name: ' % h)
            while i:
                if self.check_user_input(i):
                    return
                form = input('%s> Formula: ' % h)
                if form:
                    writer += '    %s: "%s"\n' % (i, form)
                else:
                    writer += '    %s: "FORMULA_HERE"\n' % i
                if dat not in baselines:
                    baselines[dat] = {}
                baselines[dat][i] = {}
                i = input('%sModel name: ' % h)

            if self.all_datasets:
                self.write_all_datasets(writer, 'models:\n')
                break
            else:
                if not ddx:
                    self.config.append('models:\n' + writer)
                self.config.append(writer)

        for ddx, dat in enumerate(baselines):
            mess = '%s\n%s* Define models baselines: '% (h, h)
            print(mess)
            for model in sorted(baselines[dat]):
                print('%s  - Model "%s":' % (h, model))
                name = input('%s    > Baseline name: ' % h)
                while name:
                    if name:
                        if self.check_user_input(name):
                            return
                        form = input('%s      Formula: ' % h)
                        if form:
                            baselines[dat][model][name] = form
                        else:
                            baselines[dat][model][name] = "FORMULA_HERE"
                    name = input('%s    > Baseline name: ' % h)

        if baselines:
            for ddx, dat in enumerate(baselines):
                vs = False
                writer = '  %s:\n' % dat
                for model in baselines[dat]:
                    writer += '    %s:\n' % model
                    for k, v in baselines[dat][model].items():
                        if v:
                            vs = True
                            writer += '      %s: "%s"\n' % (k, v)
                if vs:
                    if self.all_datasets:
                        self.write_all_datasets(writer, 'baselines:\n')
                        break
                    else:
                        if not ddx:
                            self.config.append('baselines:\n')
                        self.config.append(writer)

        o.write('params:\n')
        for p, v in {'batches': '20', 'epochs': "2000",
                     'thresh_feats': '0', 'thresh_samples': '0'}.items():
            mess = '%s  - Param for "%s" (INT; default=%s): ' % (h, p, v)
            inp = input(mess)
            self.config.append('  %s:\n' % p)
            if not inp:
                self.config.append('  - %s\n' % v)
            else:
                while inp:
                    if not inp.isdigit():
                        if inp:
                            print('Not an acceptable numeric: %s' % inp)
                    elif inp:
                        self.config.append('  - %s\n' % inp)
                    inp = input(mess)

        for p, v in {'learns': "1e-2", 'diff_priors': "0.5"}.items():
            mess = '%s  - Param for "%s" (FLOAT; default=%s): ' % (h, p, v)
            inp = input(mess)
            self.config.append('  %s:\n' % p)
            if not inp:
                self.config.append('  - %s\n' % v)
            else:
                while inp:
                    try:
                        float(inp)
                        self.config.append('  - "%s"\n' % inp)
                    except ValueError:
                        print('Not a float: "%s"' % inp)
                    inp = input(mess)

        for p, v in {'train': '0.7'}.items():
            mess = '%s  - Param for "%s" (FLOAT/METADATA; default=%s): ' % (
                h, p, v)
            inp = input(mess)
            self.config.append('  %s:\n' % p)
            if not inp:
                self.config.append('  - %s\n' % v)
            else:
                while inp:
                    try:
                        float(inp)
                        self.config.append('  - "%s"\n' % inp)
                    except ValueError:
                        if inp not in self.vars:
                            print('"%s" not float/metadata column' % inp)
                        else:
                            if self.vars[inp][0] in [str, bool]:
                                self.config.append('  - "%s"\n' % inp)
                    inp = input(mess)

    def prep_filt3d(self):
        filts = {'prevalCount': ['0', '1', '2', '3', '5', '10', '20', '30'],
                 'prevalPercent': ['0', '0.01', '0.02', '0.03', '0.05',
                                   '0.1', '0.2', '0.3'],
                 'abundCount': ['0', '1', '2', '5', '10', '50', '100',
                                '500', '1000'],
                 'abundPercent': ['0', '0.0001', '0.0005', '0.001', '0.005',
                                  '0.01', '0.02', '0.03', '0.05', '0.1']}
        h = '%s(All samples) ' % self.prep
        print('%sDefining filtering thresholds, based on:')
        config = []
        for fil, threshs in filts.items():
            print('%s* %s (defaults=%s): ' % (h, fil, ', '.join(threshs)))
            inp = input('%s  Use default? <[y]/n>: ' % h)
            if inp in ['QUIT', 'STOP', 'EXIT']:
                print('\t <<< Skipping config >>>')
                return True
            self.empty = False
            config.append('%s:\n' % fil)
            if not inp or inp == 'y':
                config.append('- %s\n' % '\n- '.join(map(str, threshs)))
            else:
                mess = '%s  Enter percentage threshold: ' % h
                if 'Count' in fil:
                    mess = '%s  Enter count threshold: ' % h
                inp = input(mess)
                while inp:
                    if inp in ['QUIT', 'STOP', 'EXIT']:
                        print('\t <<< Skipping config >>>')
                        return True
                    if 'Count' in fil:
                        if inp.isdigit():
                            config.append('- %s\n' % inp)
                    else:
                        try:
                            if not 0 < float(inp) < 1:
                                print('"%s" not a [0-1] float' % inp)
                            config.append('- %s\n' % inp)
                        except ValueError:
                            print('"%s" not a float' % inp)
                    inp = input(mess)
        if config:
            self.config.extend(config)

    def create_config(self):
        none = True
        for arg, boolean in list(self.__dict__.items()):
            func = 'prep_%s' % arg
            if hasattr(self, func) and callable(getattr(self, func)):
                if boolean or self.all_configs:
                    none = False
                    self.empty = True
                    self.analysis = arg
                    self.setup()
                    getattr(self, func)()
                    self.teardown()
        if none:
            print('\nNo config written? Either:')
            print('- Activate `--all-configs` to prepare all configurations')
            print('- Activate `--<STEP>` to prepare individual configurations')


    def show_config(self):
        if not self.configs:
            sys.exit('No template configuration file to prepare')
        else:
            print('\nDone! Please check/edit your template configuration files')
            for analysis, config_fp in self.configs.items():
                print('>%s\n%s' % (analysis, config_fp))

    def write_config(self):
        with open(self.config_fp, 'w') as o:
            for line in self.config:
                o.write(line)

    def teardown(self):
        if self.config:
            self.write_config()
            self.configs[self.analysis] = self.config_fp
            print('\nWritten: %s' % self.config_fp)
            self.config = []


