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
import subprocess

import pandas as pd
import pkg_resources
import numpy as np
from datetime import datetime as dt
from os.path import abspath, dirname, exists, isfile, isdir
from microbiome_analyzer._inputs import read_yaml_file, read_meta_pd

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
        if self.rarefy:
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
        self.config = ''
        self.config_fp = ''
        self.analysis = ''
        self.time = dt.now().strftime("%d-%m-%Y_%H") + 'h'
        self.out = '%s/%s' % (self.configs_folder, self.time)
        self.check_io()
        self.get_vars()

    def check_io(self):
        if not isdir(self.analysis_folder):
            sys.exit('Input folder not found: "%s"' % self.analysis_folder)
        if not isdir('%s/metadata' % self.analysis_folder):
            sys.exit('"%s/metadata" not found' % self.analysis_folder)
        for dat in self.datasets:
            meta_fp = '%s/metadata/%s.tsv' % (self.analysis_folder, dat)
            if not isfile(meta_fp):
                sys.exit('No metadata file "%s"' % meta_fp)
            self.meta[dat] = read_meta_pd(meta_fp)
        if not isdir(dirname(self.configs_folder)):
            sys.exit('Output folder not found: "%s"' % self.configs_folder)
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
        # for idx, (i, j) in enumerate(self.vars.items()):
        #     if idx == 4:
        #         break
        #     print(j[0], sorted(j[1])[:4])

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

    def check_user_input(self, inp, digit=False, num=False):
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

    def setup(self):
        self.prep = '[Config: %s]' % self.analysis
        self.config_fp = '%s/%s.yml' % (self.out, self.analysis)
        sep = '\n%s\n' % ('=' * (len(self.config_fp) + 2))
        print('%s%s\nPreparing template:\n%s' % (
            sep, self.prep, self.config_fp))
        print('<Enter for none/exit>)%s' % sep)

    def teardown(self):
        if self.empty:
            os.remove(self.config_fp)
        else:
            self.configs[self.analysis] = self.config_fp
            print('\nWritten: %s' % self.config_fp)
            time.sleep(2)

    def prep_filtering(self):
        with open(self.config_fp, 'w') as o:
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
                    self.check_user_input(inp)
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
                    self.check_user_input(inp, True, True)
                    writer += '  features: %s\n' % inp

                inp = input('%s(%s) Sample threshold: ' % (self.prep, dat))
                if inp:
                    self.check_user_input(inp, True, True)
                    writer += '  samples: %s\n' % inp

                if self.all_datasets:
                    self.write_all_datasets(writer, o)
                    break
                else:
                    o.write(writer)

    def prep_rarefations(self):
        with open(self.config_fp, 'w') as o:
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
                        self.check_user_input(inp, True)
                        writer += '  - "%s"\n' % inp
                    inp = input(mess)

                if self.all_datasets:
                    self.write_all_datasets(writer, o)
                    break
                else:
                    o.write(writer)

    def prep_feature_subsets(self):
        with open(self.config_fp, 'w') as o:
            for dat in self.datasets:
                writer = '%s:\n' % dat
                mess = '%s(%s) Enter subset name: ' % (self.prep, dat)
                inp = input(mess)
                while inp:
                    self.check_user_input(inp)
                    writer += '  %s:\n' % inp
                    mess2 = '%s(%s) Feature name/regex: ' % (self.prep, dat)
                    inp2 = input(mess2)
                    while input(inp2):
                        self.check_user_input(inp2)
                        writer += '    - %s\n' % inp2
                        inp2 = input(mess2)
                    inp = input(mess)

                if self.all_datasets:
                    self.write_all_datasets(writer, o)
                    break
                else:
                    o.write(writer)

    def get_str_bool_vals(self, inp, vals):
        groups = []
        vals_idx = {str(x): y for x, y in enumerate(vals)}
        p = '\n'.join(['%s\t: %s' % (x, y) for x, y in vals_idx.items()])
        mess = '%s\tFactors to group for "%s":' % (self.prep, inp)
        i = input('%s\n%s\nspace-separated indices: ' % (mess, p))
        while i:
            j = set(i.strip().split())
            off = j.difference(set(vals_idx))
            if off:
                sys.exit('%s\tWrong index selection: %s' % (self.prep, off))
            groups.extend([
                '  - "%s"' % vals_idx[k] if x else '- - "%s"' % vals_idx[k]
                for x, k in enumerate(j)])
            self.empty = False
            i = input('%s\n%s\nspace-separated indices: ' % (mess, p))
        return groups

    def get_int_float_vals(self, inp, vals, typ):
        groups = []
        desc = str(pd.Series(list(vals), dtype=typ).describe()).split('dtyp')[0]
        mess = '%s\tBounding values to make group for "%s":' % (self.prep, inp)
        i = input('%s\n%s\n("<" or ">" + value): ' % (mess, desc))
        while i:
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
            groups.append('- - %s' % i)
            self.empty = False
            i = input('%s\n%s\n("<" or ">" + value): ' % (mess, desc))
        return groups

    def get_subsets(self, dat):
        mess = '%s(%s) Enter metadata variable name: ' % (self.prep, dat)
        inp = input(mess)
        subsets = {}
        while inp:
            self.check_user_input(inp)
            if inp not in self.vars:
                print('No variable "%s" in metadata files' % inp)
            else:
                typ, vals = self.vars[inp]
                if typ in [str, bool]:
                    groups = self.get_str_bool_vals(inp, vals)
                else:
                    groups = self.get_int_float_vals(inp, vals, typ)
                if groups:
                    subsets[inp] = groups
            inp = input(mess)
        return subsets

    def prep_sample_subsets(self):
        groups = {}
        with open(self.config_fp, 'w') as o:
            dat = 'All datasets'
            subsets = self.get_subsets(dat)
            if subsets:
                for inp, groups in subsets.items():
                    o.write('%s:\n%s\n' % (inp, '\n'.join(groups)))

    def write_all_datasets(self, writer, o, head=''):
        if head:
            o.write(head)
        for dat in self.datasets:
            o.write(writer.replace('All datasets', dat))

    def prep_adonis(self):
        with open(self.config_fp, 'w') as o:
            stratas = {}
            for ddx, dat in enumerate(self.datasets):
                if self.all_datasets:
                    dat = 'All datasets'
                h = '%s(%s) ' % (self.prep, dat)
                writer = '  %s:\n' % dat
                print('%s* Defining models:' % h)
                i = input('%sModel name: ' % h)
                while i:
                    self.check_user_input(i)
                    form = input('%s> Formula: ' % h)
                    if form:
                        writer += '    %s: "%s"\n' % (i, form)
                    else:
                        write += '    %s: "FORMULA_NEEDED_HERE"\n' % i
                    if dat not in stratas:
                        stratas[dat] = {}
                    stratas[dat][model] = set()
                    i = input('%sModel name: ' % h)

                if self.all_datasets:
                    write_all_datasets(self, writer, o, 'models:\n')
                    break
                else:
                    if not ddx:
                        o.write('models:\n' + writer)
                    o.write(writer)

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
                        self.check_user_input(i)
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
                self.check_user_input(i)
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
                            self.write_all_datasets(writer, o, 'strata:\n')
                            if global_strata:
                                o.write('  global: "%s"\n' % global_strata)
                            break
                        else:
                            if not ddx:
                                o.write('strata:\n')
                            o.write(writer)

            if not self.all_datasets and global_strata:
                o.write('  global: "%s"\n' % global_strata)

    def prep_nestedness(self):
        with open(self.config_fp, 'w') as o:
            h = '%s(All samples) ' % self.prep
            binary = input('%s* Path to "Nestedness" java binary: ' % h)
            if binary:
                o.write('soft: "%s"\n' % binary)
            else:
                o.write('soft: "BINARY_NEEDED_HERE"\n')

            subsets = self.get_subsets('All samples')
            if subsets:
                o.write('subsets:\n')
                for inp, groups in subsets.items():
                    o.write('  %s:\n  %s\n' % (inp, '\n  '.join(groups)))

            print('%s* Samples coloring: ' % h)
            mess = '%s  Metadata variable (categorical): ' % h
            inp = input(mess)
            cols = set()
            while inp:
                self.check_user_input(inp)
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
                o.write('sample_colors:\n- "%s"\n' % (
                    '"\n- "'.join(sorted(cols))))

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
                o.write('feature_colors:\n- %s\n' % ('\n- '.join(sorted(cols))))

            print('%s* NODF variables: ' % h)
            mess = '%s  Metadata variable (categorical): ' % h
            inp = input(mess)
            cols = set()
            while inp:
                self.check_user_input(inp)
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
                o.write('nodfs:\n- %s\n' % ('\n- '.join(sorted(cols))))

            print('%s* Comparison modes: ' % h)
            modes = set()
            for m in ["overall", "withineachtype", "betweeneachpairoftypes"]:
                inp = input('%s  Use mode "%s": <[y]/n>' % (h, m))
                if not inp or inp == 'y':
                    modes.add(m)
                elif inp != 'n':
                    print('Must be "y" or "n"')
            if modes:
                o.write('modes:\n- "%s"\n' % ('"\n- "'.join(sorted(modes))))

            print('%s* Null models: ' % h)
            nulls = set()
            for n in ["equiprobablefixed", "fixedequiprobable"]:
                inp = input('%s  Use null "%s": <[y]/n>' % (h, n))
                if not inp or inp == 'y':
                    nulls.add(n)
                elif inp != 'n':
                    print('Must be "y" or "n"')
            if nulls:
                o.write('nulls:\n- "%s"\n' % ('"\n- "'.join(sorted(nulls))))

            iterations = input('%s* Number of iterations (default=20): ' % h)
            if not iterations.isdigit():
                print('Not a numeric: "%s" (using default)' % inp)
            elif iterations:
                o.write('params:\n  iterations: %s\n' % iterations)
            else:
                o.write('params:\n  iterations: 20\n')

    def get_pairs(self):

        sams_idx = {str(x): y for x, y in enumerate(self.datasets)}
        sams = '\n'.join(['%s\t: %s' % s for s in sams_idx.items()])

        h = '%s(sample pairs) ' % self.prep
        mess = '%sSamples to pair:\n%s\n' % (h, sams)
        i = input('%s%sTwo space-separated indices: ' % (mess, h))

        pairs = {}
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
                    self.check_user_input(name)
                    pairs[name] = [sams_idx[x] for x in j]
                    self.empty = False
                    break
                name = input('%sNickname for this pair: ' % h)
            i = input('%s%sTwo space-separated indices: ' % (mess, h))
        return pairs

    def prep_procrustes(self):
        with open(self.config_fp, 'w') as o:
            pairs = self.get_pairs()
            if pairs:
                o.write('pairs:\n')
                for pair, sams in pairs.items():
                    o.write('  %s:\n' % pair)
                    for sam in sams:
                        o.write('  - "%s"\n' % sam)

    def prep_mantel(self):
        with open(self.config_fp, 'w') as o:
            pairs = self.get_pairs()
            if pairs:
                o.write('pairs:\n')
                for pair, sams in pairs.items():
                    o.write('  %s:\n' % pair)
                    for sam in sams:
                        o.write('  - "%s"\n' % sam)

    def prep_dm_decay(self):
        with open(self.config_fp, 'w') as o:
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
                        self.check_user_input(v)
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
                o.write('modes:\n')
                print(modes)
                for m, vs in modes.items():
                    o.write('  %s:\n  - "%s"\n' % (
                        m, '"\n  - "'.join(sorted(vs))))

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
                o.write('params:\n')
                for p, v in params.items():
                    o.write('  %s: %s\n' % (p, v))

    def prep_geo_decay(self):
        print('NOT AVAILABLE YET :(')

    def prep_collapse(self):
        with open(self.config_fp, 'w') as o:
            for dat in self.datasets:
                if self.all_datasets:
                    dat = 'All datasets'
                h = '%s(%s) ' % (self.prep, dat)
                mess = '%sEnter a collapse name: ' % h
                collapses = {}
                inp = input(mess)
                while inp:
                    if inp:
                        self.check_user_input(inp)
                        name = input('%s> Taxon identifier (or index): ' % h)
                        while 1:
                            if name:
                                collapses[inp] = name
                                break
                            name = input('%sNickname for this pair: ' % h)
                    inp = input(mess)
                if collapses:
                    writer = '%s:\n' % dat
                    for i, j in collapses.items():
                        writer += '  %s: "%s"\n' % (i, j)

                    if self.all_datasets:
                        self.write_all_datasets(writer, o)
                        break
                    else:
                        o.write(writer)

    def prep_train_test(self):
        with open(self.config_fp, 'w') as o:
            trains = {}
            for dat in self.datasets:
                if self.all_datasets:
                    dat = 'All datasets'
                h = '%s(%s) ' % (self.prep, dat)
                print('%s* Set training/testing splits:' % h)
                i = input('%s  Enter name of a split: ' % h)
                while i:
                    self.check_user_input(i)
                    inp = input('%s  Metadata variable: ' % h)
                    while inp:
                        self.check_user_input(inp)
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
                    writer += '  %s:\n' % dat
                    for k, vs in trains[dat].items():
                        writer += '    %s:\n' % k
                        for v in vs:
                            writer += '    - "%s"\n' % v
                    if self.all_datasets:
                        self.write_all_datasets(writer, o, 'datasets:\n')
                        break
                    else:
                        if not ddx:
                            o.write('datasets:\n')
                        o.write(writer)

            inp = input('%s* Fraction for training (default=0.7): ' % h)
            try:
                float(inp)
                o.write('train: %s\n' % inp)
            except ValueError:
                print('Training input must be a float (using 0.7)')
                o.write('train: 0.7\n')

    def prep_phate(self):
        with open(self.config_fp, 'w') as o:
            h = '%s(All samples) ' % self.prep
            print('%s* Labels for coloring: ' % h)
            inp = input('%s  - Metadata variable: ' % h)
            cols = set()
            while inp:
                self.check_user_input(inp)
                if inp not in self.vars:
                    print('No variable "%s" in metadata files' % inp)
                else:
                    cols.add(inp)
                    # typ, vals = self.vars[inp]
                    # if typ != str:
                    #     print('No variable "%s" not categorical' % inp)
                    # else:
                    #     cols.add(inp)
                inp = input('%s  - Metadata variable: ' % h)
            if cols:
                o.write('labels:\n- "%s"\n' % ('"\n- "'.join(sorted(cols))))

            print('%s* Setup PHATE key parameters: ' % h)
            o.write('params:\n')
            for p, (s, e, t) in {'t': [1, 20, 2],
                                 'decay': [5, 25, 5],
                                 'knn': [3, 30, 3]}.items():
                o.write('  %s:\n' % p)
                d = '%s (start=%s; end=%s; step=%s)' % (
                    str(list(range(s, (e + 1), t))), s, e, t)
                print('%s  - Parameters for "%s" (defaults: %s): ' % (h, p, d))
                for step, val in [('start', s), ('end', e), ('step', t)]:
                    inp = input('%s    * %s [default=%s]: ' % (h, step, val))
                    if not inp.isdigit():
                        print('Not a numeric: "%s" (using default)' % inp)
                        o.write('  - %s\n' % val)
                    elif inp:
                        o.write('  - %s\n' % inp)

    def prep_sourcetracking(self):
        pass

    def prep_doc(self):
        with open(self.config_fp, 'w') as o:
            h = '%s(All samples) ' % self.prep
            print('%s* Setup DOC parameters: ' % h)
            o.write('params:\n')
            for p, v in {'r': 2, 'subr': 50, 'mov_avg': 5, 'degree': 1,
                               'iterations': 2, 'nulls': 1, 'non_zero': 1,
                               'null': 1}.items():
                mess = '%s  - Param for "%s" (INT; default=%s): ' % (h, p, v)
                inp = input(mess)
                if not inp.isdigit():
                    if inp:
                        print('Not an acceptable numeric: %s' % inp)
                    o.write('  %s: %s\n' % (p, v))
                elif inp:
                    o.write('  %s: %s\n' % (p, inp))
                self.empty = False

            for p, v in {'span': 0.2}.items():
                mess = '%s  - Param for "%s" (FLOAT; default=%s): ' % (h, p, v)
                inp = input(mess)
                if not inp.isdigit():
                    if inp:
                        print('Not an acceptable numeric: %s' % inp)
                    o.write('  %s: %s\n' % (p, v))
                elif inp:
                    o.write('  %s: %s\n' % (p, inp))
                self.empty = False

            for p, v in {'family': 'symmetric', 'use_mp': 'false',
                         'surface': 'interpolate'}.items():
                mess = '%s  - Param for "%s" (TEXT default=%s): ' % (h, p, v)
                inp = input(mess)
                if inp:
                    o.write('  %s: %s\n' % (p, inp))
                else:
                    o.write('  %s: %s\n' % (p, v))
                self.empty = False

            print('%s  - Param for "ci" (BOUNDS): ')
            o.write('  ci:\n')
            for n, v in [('low', '0.025'), ('mid', '0.5'), ('top', '0.975')]:
                inp = input('%s    * %s [default=%s]: ' % (h, n, v))
                try:
                    float(inp)
                    o.write('  - "%s"\n' % inp)
                except ValueError:
                    print('Not a float: "%s"' % inp)
                    o.write('  - "%s"\n' % v)
                self.empty = False

            inp = input('%s* Use also the R version? <[y]/n>: ' % h)
            if not inp or inp == 'y':
                o.write('do_r: 1\n')
            else:
                o.write('do_r: 0\n')
            self.empty = False

    def prep_xmmvec(self):
        pass

    def prep_mmvec_highlights(self):
        pass

    def prep_mmvec_pairs(self):
        pass

    def prep_diff_models(self):
        with open(self.config_fp, 'w') as o:
            baselines = {}
            for dat in self.datasets:
                if self.all_datasets:
                    dat = 'All datasets'
                h = '%s(%s) ' % (self.prep, dat)
                writer = '  %s:\n' % dat
                print('%s* Defining models:' % h)
                i = input('%sModel name: ' % h)
                while i:
                    self.check_user_input(i)
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
                    write_all_datasets(self, writer, o, 'models:\n')
                    break
                else:
                    if not ddx:
                        o.write('models:\n' + writer)
                    o.write(writer)

            for ddx, dat in enumerate(baselines):
                mess = '%s\n%s* Define models baselines: '% (h, h)
                print(mess)
                for model in sorted(baselines[dat]):
                    print('%s  - Model "%s":' % (h, model))
                    name = input('%s    > Baseline name: ' % h)
                    while name:
                        if name:
                            self.check_user_input(name)
                            form = input('%s      Formula: ' % h)
                            if form:
                                baselines[dat][model][name] = form
                            else:
                                baselines[dat][model][name] = "FORMULA_HERE"
                        name = input('%s    > Baseline name: ' % h)

            if baselines:
                for ddx, dat in enumerate(baselines):
                    writer = '  %s:\n' % dat
                    for model in baselines[dat]:
                        writer += '    %s:\n' % model
                        for k, v in baselines[dat][model].items():
                            writer += '      %s: "%s"\n' % (k, v)

                        if self.all_datasets:
                            self.write_all_datasets(writer, o, 'baselines:\n')
                        else:
                            if not ddx:
                                o.write('baselines:\n')
                            o.write(writer)

            o.write('params:\n')
            for p, v in {'batches': '20', 'epochs': "2000",
                         'thresh_feats': '0', 'thresh_samples': '0'}.items():
                mess = '%s  - Param for "%s" (INT; default=%s): ' % (h, p, v)
                inp = input(mess)
                o.write('  %s:\n' % p)
                if not inp:
                    o.write('  - %s\n' % v)
                else:
                    while inp:
                        if not inp.isdigit():
                            if inp:
                                print('Not an acceptable numeric: %s' % inp)
                        elif inp:
                            o.write('  - %s\n' % inp)
                        inp = input(mess)

            for p, v in {'learns': "1e-2", 'diff_priors': "0.5"}.items():
                mess = '%s  - Param for "%s" (FLOAT; default=%s): ' % (h, p, v)
                inp = input(mess)
                o.write('  %s:\n' % p)
                if not inp:
                    o.write('  - %s\n' % v)
                else:
                    while inp:
                        try:
                            float(inp)
                            o.write('  - "%s"\n' % inp)
                        except ValueError:
                            print('Not a float: "%s"' % inp)
                        inp = input(mess)

            for p, v in {'train': '0.7'}.items():
                mess = '%s  - Param for "%s" (FLOAT/METADATA; default=%s): ' % (
                    h, p, v)
                inp = input(mess)
                o.write('  %s:\n' % p)
                if not inp:
                    o.write('  - %s\n' % v)
                else:
                    while inp:
                        try:
                            float(inp)
                            o.write('  - "%s"\n' % inp)
                        except ValueError:
                            if inp not in self.vars:
                                print('"%s" not float/metadata column' % inp)
                            else:
                                if self.vars[inp][0] in [str, bool]:
                                    o.write('  - "%s"\n' % inp)
                        inp = input(mess)

    def prep_filt3d(self):
        filts = {'prevalCount': ['0', '1', '2', '3', '5', '10', '20', '30'],
                 'prevalPercent': ['0', '0.01', '0.02', '0.03', '0.05',
                                   '0.1', '0.2', '0.3'],
                 'abundCount': ['0', '1', '2', '5', '10', '50', '100',
                                '500', '1000'],
                 'abundPercent': ['0', '0.0001', '0.0005', '0.001', '0.005',
                                  '0.01', '0.02', '0.03', '0.05', '0.1']}
        with open(self.config_fp, 'w') as o:
            h = '%s(All samples) ' % self.prep
            print('%sDefining filtering thresholds, based on:')
            for fil, threshs in filts.items():
                o.write('%s:\n' % fil)
                print('%s* %s (defaults=%s): ' % (h, fil, str(list(threshs))))
                inp = input('%s  Use default? <[y]/n>: ' % h)
                self.empty = False
                if not inp or inp == 'y':
                    o.write('- %s\n' % '\n- '.join(map(str, threshs)))
                else:
                    mess = '%s  Enter percentage threshold: ' % h
                    if 'Count' in fil:
                        mess = '%s  Enter count threshold: ' % h
                    inp = input(mess)
                    while inp:
                        if 'Count' in fil:
                            if inp.isdigit():
                                o.write('- %s\n' % inp)
                        else:
                            try:
                                if not 0 < float(inp) < 1:
                                    print('"%s" not a [0-1] float' % inp)
                                o.write('- %s\n' % inp)
                            except ValueError:
                                print('"%s" not a float' % inp)
                        inp = input(mess)

    def create_config(self):
        for arg, boolean in list(self.__dict__.items()):
            func = 'prep_%s' % arg
            if hasattr(self, func) and callable(getattr(self, func)):
                if not boolean:
                    continue
                self.empty = True
                self.analysis = arg
                self.setup()
                getattr(self, func)()
                self.teardown()

    def show_config(self):
        if not self.configs:
            sys.exit('No template configuration file to prepare')
        else:
            print('\nDone! Please check/edit your template configuration files')
            for analysis, config_fp in self.configs.items():
                print('>%s\n%s' % (analysis, config_fp))


