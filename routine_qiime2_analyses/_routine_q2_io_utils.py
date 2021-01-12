# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import sys
import yaml
import glob
import pkg_resources
import pandas as pd
import numpy as np
from biom import load_table

from pandas.util import hash_pandas_object
from os.path import basename, dirname, splitext, isfile, isdir, abspath

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_cmds import run_import, run_export, get_case, get_new_meta_pd
from routine_qiime2_analyses._routine_q2_metadata import check_metadata_cases_dict

RESOURCES = pkg_resources.resource_filename("routine_qiime2_analyses", "resources")


def summarize_songbirds(i_datasets_folder) -> pd.DataFrame:
    q2s = []
    songbird_folder = get_analysis_folder(i_datasets_folder, 'songbird')
    for root, dirs, files in os.walk(songbird_folder):
        for fil in files:
            if fil == 'tensorboard.html':
                path = root + '/' + fil
                diff = '%s/differentials.tsv' % dirname(root)
                root_split = root.split('%s/' % songbird_folder)[-1].split('/')
                if len(root_split) == 8:
                    dat, pair, dataset_filter, subset, songbird_filter, parameters, model, baseline = root_split
                else:
                    pair = 'no_pair'
                    dat, dataset_filter, subset, songbird_filter, parameters, model, baseline = root_split
                with open(path) as f:
                    for line in f:
                        if 'Pseudo Q-squared' in line:
                            q2s.append([
                                pair, dat, dataset_filter, subset, model, songbird_filter, parameters, baseline, diff,
                                float(line.split('Pseudo Q-squared:</a></strong> ')[-1].split('<')[0])
                            ])
    q2s_pd = pd.DataFrame(q2s, columns=['pair', 'dat', 'dataset_filter', 'subset', 'model',
                                        'songbird_filter', 'parameters', 'baseline',
                                        'differentials', 'Pseudo_Q_squared'])
    return q2s_pd


def get_songbird_outputs(songbird_outputs: list) -> pd.DataFrame:
    songbird_outputs_pd = pd.DataFrame(
        songbird_outputs,
        columns=[
            'songbird_dat',
            'songbird_filt',
            'songbird_parameters',
            'songbird_case',
            'songbird_fp',
            'songbird_baseline',
            'songbird_q2',
            'pair'
        ])
    songbird_outputs_pd['pair_case_omic_filt'] = songbird_outputs_pd['pair'] + '__' + songbird_outputs_pd[
        'songbird_case'] + '__' + songbird_outputs_pd['songbird_dat'] + '__' + songbird_outputs_pd['songbird_filt']
    songbird_outputs_pd['params'] = songbird_outputs_pd['songbird_parameters']
    songbird_outputs_pd['baseline'] = songbird_outputs_pd['songbird_baseline']
    songbird_outputs_drop_pd = songbird_outputs_pd.drop(
        columns=['pair', 'songbird_case', 'songbird_dat', 'songbird_filt',
                 'songbird_parameters', 'songbird_baseline'])

    songbird_outputs_pd = songbird_outputs_drop_pd[
        ['params', 'pair_case_omic_filt', 'songbird_fp']
    ].drop_duplicates().pivot(
        columns='params', index='pair_case_omic_filt'
    )

    # songbird_outputs_pd['pair_omic_filt'] = songbird_outputs_pd['pair'] + '__' + songbird_outputs_pd[
    #     'songbird_dat'] + '__' + songbird_outputs_pd['songbird_filt']
    # songbird_outputs_pd['case_params'] = songbird_outputs_pd['songbird_case'] + '__' + songbird_outputs_pd[
    #     'songbird_parameters']
    # songbird_outputs_pd['case_params_baseline'] = songbird_outputs_pd[
    #     'case_params'] + '__' + songbird_outputs_pd['songbird_baseline']
    # songbird_outputs_drop_pd = songbird_outputs_pd.drop(
    #     columns=['pair', 'songbird_dat', 'songbird_filt', 'songbird_case',
    #              'songbird_parameters', 'songbird_baseline'])
    #
    # songbird_outputs_pd = songbird_outputs_drop_pd[
    #     ['case_params', 'pair_omic_filt', 'songbird_fp']
    # ].drop_duplicates().pivot(
    #     columns='case_params', index='pair_omic_filt'
    # )
    songbird_outputs_pd.columns = songbird_outputs_pd.columns.droplevel()
    songbird_outputs_pd = songbird_outputs_pd.reset_index()
    return songbird_outputs_pd


def read_yaml_file(p_yaml_file: str) -> dict:
    """
    :param p_subsets: Subsets for alpha diversity.
    """
    if p_yaml_file:
        if not isfile(p_yaml_file):
            print('[Warning] yaml file for subsets does not exist: %s\n' % p_yaml_file)
        else:
            with open(p_yaml_file) as handle:
                try:
                    yaml_content = yaml.load(handle, Loader=yaml.FullLoader)
                except AttributeError:
                    yaml_content = yaml.load(handle)
            if not isinstance(yaml_content, dict):
                print('[Warning] %s must be a dictionary\n' % p_yaml_file)
            else:
                return yaml_content
    return {}


def update_filtering_abundance(mmvec_dict: dict, p_mmvec_pairs: str, filtering: dict) -> dict:
    if 'abundance' in mmvec_dict['filtering']:
        if not isinstance(mmvec_dict['filtering']['abundance'], list):
            print('Filtering parameter "abundance" should be a list (see %s)\n' % p_mmvec_pairs)
            sys.exit(0)
        not_int = []
        for abundance_list in mmvec_dict['filtering']['abundance']:
            if not isinstance(abundance_list, list):
                print('Filtering parameter "abundance" should be a list of lists (see %s)\nExiting\n' % p_mmvec_pairs)
                sys.exit(0)
            for abundance in abundance_list:
                try:
                    int(abundance)
                except ValueError:
                    not_int.append(abundance_list)
                    break
        if not_int:
            print('Filtering parameter "abundance" should contain two integers (see %s)\n' % p_mmvec_pairs)
            print('  Not integer(s) in\n:%s\nExiting\n' % ', '.join(["['%s']" % "', '".join(x) for x in not_int]))
            sys.exit(0)
    else:
        print('No "abundance" filter specified in %s:\nUsing defaults:' % p_mmvec_pairs)
        for k, v in filtering['abundance'].items():
            print(' -', k, ' ,'.join(v))
        return {'filtering': {'abundance': filtering['abundance']}}
    return {}


def update_filtering_prevalence(filtering_dict: dict, p_yml: str, filtering: dict) -> dict:
    if 'prevalence' in filtering_dict['filtering']:
        if not isinstance(filtering_dict['filtering']['prevalence'], list):
            print('Filtering parameter "prevalence" should be a list (see %s)\nExiting\n' % p_yml)
            sys.exit(0)
        not_int = []
        for prevalence in filtering_dict['filtering']['prevalence']:
            try:
                int(prevalence)
            except ValueError:
                not_int.append(prevalence)
        if not_int:
            print('Filtering parameter "prevalence" should contain integers (see %s)\n' % p_yml)
            print('  Not integer(s)\n:%s\nExiting\n' % ', '.join(not_int))
            sys.exit(0)
    else:
        print('No "prevalence" filter specified in %s:\nUsing defaults: %s' %
              (p_yml, ','.join(filtering['prevalence'])))
        return {'filtering': {'prevalence': filtering['prevalence']}}
    return {}


def get_dat_mb_or_not(dat: str) -> tuple:
    if dat[-1] == '*':
        return (dat[:-1], 1)
    else:
        return (dat, 0)


def get_filtering(p_yml: str, filtering_dict: dict,
                  songbird_mmvec: dict, analysis: str) -> dict:
    """
    Get the parameters for songbird passed by the user.
    :param p_mmvec_pairs: file containing the parameters.
    :param mmvec_dict: parsed content of the file containing the parameters.
    :return: parameters.
    """
    dats = []
    filtering = {}
    if analysis == 'mmvec':
        for pair, dats_pair in songbird_mmvec.items():
            if pair not in filtering:
                filtering[pair] = {'0_0': {}}
            for dat in dats_pair:
                dats.append(dat)
                filtering[pair]['0_0'][dat] = ['0', '0']
    elif analysis == 'songbird':
        filtering[''] = {'0_0': {}}
        for dat_ in songbird_mmvec.keys():
            dat = get_dat_mb_or_not(dat_)
            dats.append(dat)
            filtering['']['0_0'][dat] = ['0', '0']

    if 'filtering' not in filtering_dict:
        print('No filtering thresholds set in %s\n:' % p_yml)
    elif analysis == 'mmvec':
        if 'global' in filtering_dict['filtering']:
            for filt_name, prev_abund in filtering_dict['filtering']['global'].items():
                for pair, dats_pair in songbird_mmvec.items():
                    if pair not in filtering:
                        filtering[pair] = {}
                    if filt_name not in filtering[pair]:
                        filtering[pair][filt_name] = {}
                    for dat in dats_pair:
                        dats.append(dat)
                        filtering[pair][filt_name][dat] = prev_abund
        for pair, pair_d in filtering_dict['filtering'].items():
            if pair == 'global':
                continue
            filtering[pair] = {}
            for filt_name, dats_d in pair_d.items():
                filtering[pair][filt_name] = {}
                for dat_, prev_abund in dats_d.items():
                    dat = get_dat_mb_or_not(dat_)
                    if dat in dats:
                        filtering[pair][filt_name][dat] = prev_abund

    elif analysis == 'songbird':
        if 'global' in filtering_dict['filtering']:
            for filt_name, prev_abund in filtering_dict['filtering']['global'].items():
                for dat_ in songbird_mmvec.keys():
                    dat = get_dat_mb_or_not(dat_)
                    if filt_name not in filtering['']:
                        filtering[''][filt_name] = {}
                    filtering[''][filt_name][dat] = prev_abund

        for dat_, filts in filtering_dict['filtering'].items():
            if dat_ == 'global':
                continue
            dat = get_dat_mb_or_not(dat_)
            for filt_name, prev_abund in filts.items():
                if filt_name not in filtering['']:
                    filtering[''][filt_name] = {}
                if dat in dats:
                    filtering[''][filt_name][dat] = prev_abund

    return filtering


def get_mmvec_params(p_mmvec_pairs: str, mmvec_dict: dict) -> dict:
    """
    Get the parameters for songbird passed by the user.
    :param p_mmvec_pairs: file containing the parameters.
    :param mmvec_dict: parsed content of the file containing the parameters.
    :return: parameters.
    """
    params = {
        'train_column': ['None'],
        'n_examples': ['10'],
        'batches': ['2'],
        'learns': ['1e-4'],
        'epochs': ['5000'],
        'priors': ['0.1', '1'],
        'thresh_feats': ['0'],
        'latent_dims': ['3']
    }
    if 'params' not in mmvec_dict:
        print('No parameters set in %s:\nUsing defaults: %s' % (
            p_mmvec_pairs, ', '.join(['%s: %s' % (k, v) for k,v in params.items()])))
    else:
        for param, cur_param in mmvec_dict['params'].items():
            if not isinstance(cur_param, list):
                print('Parameter %s should be a list (correct in %s)\n' % (param, p_mmvec_pairs))
                sys.exit(0)
            params[param] = cur_param
    return params


def get_mmvec_subsets(mmvec_pairs: dict, mmvec_dict: dict) -> dict:
    mmvec_subsets = {'ALL': [[]]}
    if 'subsets' in mmvec_dict:
        mmvec_subsets.update(mmvec_dict['subsets'])
    return mmvec_subsets


def get_mmvec_pairs(p_mmvec_pairs: str, mmvec_dict: dict) -> dict:
    """
    Get the parameters for mmvec pairs to process.
    :param p_mmvec_pairs: file containing the parameters.
    :param mmvec_dict: parsed content of the file containing the parameters.
    :return: parameters.
    """
    if 'pairs' not in mmvec_dict:
        print('No datasets pairs specified in %s:\nExiting\n' % p_mmvec_pairs)
        sys.exit(0)
    mmvec_pairs = {}
    for pair, paired_datasets in mmvec_dict['pairs'].items():
        n_dats = len(paired_datasets)
        if n_dats != 2:
            print('Must be two datasets per mmvec pair (found %s in %s)\n'
                  'Exiting\n' % (n_dats, p_mmvec_pairs))
            sys.exit(0)
        paired = []
        for dat in paired_datasets:
            if dat[-1] == '*':
                paired.append((dat[:-1], 1))
            else:
                paired.append((dat, 0))
        mmvec_pairs[pair] = paired
    return mmvec_pairs


def get_mmvec_dicts(p_mmvec_pairs: str) -> (dict, dict, dict, dict):
    """
    Collect pairs of datasets from the passed yaml file:
    :param p_mmvec_pairs: Pairs of datasets for which to compute co-occurrences probabilities.
    :param datasets_filt:
    :return: datasets pairs, filtering thresholds, mmvec run parameters.
    """
    if not isfile(p_mmvec_pairs):
        print('yaml file for mmvec pairs does not exist:\n%s\nExiting...' % p_mmvec_pairs)
        sys.exit(0)
    with open(p_mmvec_pairs) as handle:
        try:
            mmvec_dict = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            mmvec_dict = yaml.load(handle)

    mmvec_pairs = get_mmvec_pairs(p_mmvec_pairs, mmvec_dict)
    mmvec_filtering = get_filtering(p_mmvec_pairs, mmvec_dict, mmvec_pairs, 'mmvec')
    mmvec_params = get_mmvec_params(p_mmvec_pairs, mmvec_dict)
    mmvec_subsets = get_mmvec_subsets(mmvec_pairs, mmvec_dict)
    return mmvec_pairs, mmvec_filtering, mmvec_params, mmvec_subsets


def get_songbird_params(p_diff_models: str, diff_dict: dict) -> dict:
    """
    Get the parameters for songbird passed by the user.
    :param p_diff_models: file containing the parameters.
    :param diff_dict: parsed content of the file containing the parameters.
    :return: parameters.
    """
    params = {
        'train': ['0.2'],
        'batches': ['2'],
        'learns': ['1e-4'],
        'epochs': ['5000'],
        'thresh_feats': ['0'],
        'thresh_samples': ['0'],
        'diff_priors': ['0.5']
    }
    if 'params' not in diff_dict:
        print('No parameters set in %s:\nUsing defaults: %s' % (
            diff_dict, ', '.join(['%s: %s' % (k,v) for k,v in params.items()])))
    else:
        for param in diff_dict['params']:
            cur_param = diff_dict['params'][param]
            if not isinstance(cur_param, list):
                print('Parameter %s should be a list (correct in %s)\n' % (param, p_diff_models))
                sys.exit(0)
            params[param] = cur_param
    return params


def get_songbird_baselines(diff_dict: dict) -> dict:
    """
    Get potential different baselines for songbird modesl.
    :param p_diff_models: file containing the parameters.
    :param diff_dict: parsed content of the file containing the parameters.
    :return: baselines per model.
    """
    baselines = {}
    if 'baselines' in diff_dict:
        return diff_dict['baselines']
    return baselines


def get_songbird_models(p_diff_models: str, diff_dict: dict) -> dict:
    """
    Get the models for songbird passed by the user.
    :param p_diff_models: file containing the models.
    :param diff_dict: parsed content of the file containing the models.
    :return: models.
    """
    if 'models' not in diff_dict:
        print('No models in %s' % p_diff_models)
        sys.exit(0)
    return diff_dict['models']


def get_highlights_mmbird(highlights_mmbird_fp: str) -> dict:
    highlights_mmbird = {}
    if isfile(highlights_mmbird_fp):
        with open(highlights_mmbird_fp) as handle:
            try:
                highlights_mmbird = yaml.load(handle, Loader=yaml.FullLoader)
            except AttributeError:
                highlights_mmbird = yaml.load(handle)

    return highlights_mmbird


def get_songbird_dicts(p_diff_models: str) -> (dict, dict, dict, dict, dict, dict):
    """
    Collect from on the passed yaml file:
    - subsets to perform songbird on
    - formulas for songbird
    - paramters for the modelling.
    :param p_perm_groups: path to the yaml file containing groups.
    :return: subset groups.
    """
    if not isfile(p_diff_models):
        print('yaml file containing groups does not exist:\n%s\nExiting...' % p_diff_models)
        sys.exit(0)
    with open(p_diff_models) as handle:
        try:
            diff_dict = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            diff_dict = yaml.load(handle)

    main_cases_dict = {'ALL': [[]]}
    if 'subsets' in diff_dict:
        main_cases_dict.update(diff_dict['subsets'])

    models = get_songbird_models(p_diff_models, diff_dict)
    params = get_songbird_params(p_diff_models, diff_dict)
    baselines = get_songbird_baselines(diff_dict)
    filtering = get_filtering(p_diff_models, diff_dict,  models, 'songbird')
    datasets = [(dat[:-1], 1) if dat[-1] == '*' else (dat, 0) for dat in models]
    return models, filtering, params, baselines, datasets, main_cases_dict


def get_phate_dicts(p_phate_config: str) -> (dict, list, dict, dict):
    if not isfile(p_phate_config):
        print('yaml file containing groups does not exist:\n%s\nExiting...' % p_phate_config)
        sys.exit(0)
    with open(p_phate_config) as handle:
        try:
            diff_dict = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            diff_dict = yaml.load(handle)

    main_cases_dict = {'ALL': [[]]}
    if 'subsets' in diff_dict:
        main_cases_dict.update(diff_dict['subsets'])

    phate_filtering = {}
    if 'filtering' in diff_dict:
        phate_filtering = diff_dict['filtering']

    phate_params = {'t': (None,), 'decay': (15,), 'knn': (5,)}
    if 'parameters' in diff_dict:
        phate_params.update(diff_dict['parameters'])

    phate_labels = []
    if 'labels' in diff_dict:
        phate_labels.extend(diff_dict['labels'])

    return phate_filtering, phate_labels, phate_params, main_cases_dict


def get_doc_config(p_doc_config: str) -> (dict, dict, dict):
    doc_filtering = {}
    doc_params = {
        'r': 50,
        'subr': 0,
        'mov_avg': 5,
        'ci': ['0.025', '0.5', '0.975'],
        'span': 0.2,
        'degree': 1,
        'family': 'symmetric',
        'iterations': 2,
        'surface': 'direct',
        'nulls': 1,
        'non_zero': 1,
        'null': 1,
        'use_mp': False
    }
    main_cases_dict = {'ALL': [[]]}
    if p_doc_config:
        if not isfile(p_doc_config):
            print('DOC config yaml file does not exist:\n%s\nExiting...' % p_doc_config)
            sys.exit(0)
        with open(p_doc_config) as handle:
            try:
                doc_config = yaml.load(handle, Loader=yaml.FullLoader)
            except AttributeError:
                doc_config = yaml.load(handle)

        if 'filtering' in doc_config:
            doc_filtering = doc_config['filtering']
        if 'subsets' in doc_config:
            main_cases_dict.update(doc_config['subsets'])
        if 'params' in doc_config:
            doc_params.update(doc_config['params'])
    return doc_filtering, doc_params, main_cases_dict


def get_sourcetracking_config(p_sourcetracking_config: str) -> (dict, dict, dict, dict):
    sourcetracking_filtering = {}
    sourcetracking_params = {
        'method': ['q2'],
        'times': 1,
        'iterations': None,
        'rarefaction': None
    }
    main_cases_dict = {'ALL': [[]]}
    if not p_sourcetracking_config or not isfile(p_sourcetracking_config):
        print('DOC config yaml file does not exist:\n%s\nExiting...' % p_sourcetracking_config)
        sys.exit(0)
    with open(p_sourcetracking_config) as handle:
        try:
            sourcetracking_config = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            sourcetracking_config = yaml.load(handle)

    if 'sourcesink' not in sourcetracking_config:
        raise IOError('At least one sink for one metadata column must be set '
                      '(no "sourcesink" in %s)' % p_sourcetracking_config)
    sourcetracking_sourcesink = sourcetracking_config['sourcesink']
    if 'filtering' in sourcetracking_config:
        sourcetracking_filtering = sourcetracking_config['filtering']
    if 'subsets' in sourcetracking_config:
        main_cases_dict.update(sourcetracking_config['subsets'])
    if 'params' in sourcetracking_config:
        sourcetracking_params.update(sourcetracking_config['params'])
    return sourcetracking_sourcesink, sourcetracking_filtering, sourcetracking_params, main_cases_dict


def get_main_cases_dict(p_perm_groups: str) -> dict:
    """
    Collect the subsets to perform based on the passed yaml file.

    :param p_perm_groups: path to the yaml file containing groups.
    :return: subset groups.
    """
    main_cases_dict = {'ALL': [[]]}
    if p_perm_groups:
        if not isfile(p_perm_groups):
            print('yaml file containing groups does not exist:\n%s\nExiting...' % p_perm_groups)
            sys.exit(0)
        with open(p_perm_groups) as handle:
            try:
                main_cases_dict.update(yaml.load(handle, Loader=yaml.FullLoader))
            except AttributeError:
                main_cases_dict.update(yaml.load(handle))

    return main_cases_dict


def get_formulas_dict(p_formulas: str) -> dict:
    """
    Collect the formulas for ADONIS from the passed yaml file.

    :param p_formulas: path to the yaml file containing fomulas.
    :return: formulas.
    """
    if not isfile(p_formulas):
        print('yaml file containing formulas does not exist:\n%s\nExiting...' % p_formulas)
        sys.exit(0)
    with open(p_formulas) as handle:
        try:
            formulas = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            formulas = yaml.load(handle)

    return formulas


def get_feature_sample_col(meta_tab: str) -> str:

    """
    Get the first column of the metadata file.

    :param meta: metadata file name
    :return: column name
    """
    n = 0
    with open(meta_tab) as f:
        for line in f:
            n += 1
            break
    if n:
        return line.split()[0]
    else:
        print('Empty now: %s (possily being written elsewhere..)' % meta_tab)
        sys.exit(0)


def get_raref_tab_meta_pds(meta: str, tsv: str) -> (pd.DataFrame, pd.DataFrame):
    """
                # --> datasets_read <--
                # path_pd : indexed with feature name
                # meta_pd : not indexed -> "sample_name" as first column

    :param meta: metadata for non-rarefied data.
    :param tsv: rarefied table.
    :return:
    """
    tsv_pd = pd.read_csv(tsv, header=0, index_col=0, sep='\t', low_memory=False)
    meta_pd = read_meta_pd(meta)
    meta_raref_pd = meta_pd.loc[meta_pd.sample_name.isin(tsv_pd.columns.tolist()), :].copy()
    meta_raref_pd.to_csv(meta, index=False, sep='\t')
    return tsv_pd, meta_raref_pd


def read_meta_pd(meta_tab: str, rep_col: str ='sample_name') -> pd.DataFrame:
    """
    Read metadata wit first column as index.
    :param meta: file path to the metadata file.
    :return: metadata table.
    """
    meta_tab_sam_col = get_feature_sample_col(meta_tab)
    meta_tab_pd = pd.read_csv(meta_tab, header=0, sep='\t', dtype={meta_tab_sam_col: str}, low_memory=False)
    meta_tab_pd.rename(columns={meta_tab_sam_col: rep_col}, inplace=True)
    return meta_tab_pd


def simple_chunks(run_pbs, job_folder2, to_chunk, analysis: str,
                  prjct_nm: str, time: str, n_nodes: str, n_procs: str,
                  mem_num: str, mem_dim: str, qiime_env: str, chmod: str,
                  noloc: bool, jobs: bool, chunkit: int, tmp: str = None) -> None:

    warning = 0
    with open(run_pbs, 'w') as main_o:

        chunks = {}
        if chunkit and len(to_chunk) > chunkit:
            for idx, keys in enumerate(np.array_split(to_chunk, chunkit)):
                head_sh = '%s/%s_chunk%s.sh' % (job_folder2, analysis, idx)
                chunks[(idx, head_sh)] = sorted(keys)
        else:
            chunks = dict(
                ((idx, '%s/%s_chunk%s.sh' % (job_folder2, analysis, idx)), [x]) for idx, x in enumerate(to_chunk))

        for (dat, out_sh), cur_shs in chunks.items():
            cur_written = False
            with open(out_sh, 'w') as sh:
                for cur_sh in cur_shs:
                    if isfile(cur_sh):
                        with open(cur_sh) as f:
                            for line in f:
                                sh.write(line)
                                cur_written = True
                        os.remove(cur_sh)
            if jobs:
                if cur_written:
                    out_pbs = '%s.pbs' % splitext(out_sh)[0]
                    run_xpbs(out_sh, out_pbs, '%s.%s' % (prjct_nm, dat), qiime_env,
                             time, n_nodes, n_procs, mem_num, mem_dim, chmod, 1,
                             '', None, noloc, jobs, tmp)
                    if os.getcwd().startswith('/panfs'):
                        out_pbs = out_pbs.replace(os.getcwd(), '')
                    main_o.write('qsub %s\n' % out_pbs)
                    warning += 1
                else:
                    os.remove(out_sh)
            else:
                if cur_written:
                    main_o.write('sh %s\n' % out_sh)
                    warning += 1


def write_main_sh(job_folder: str, analysis: str, all_sh_pbs: dict,
                  prjct_nm: str, time: str, n_nodes: str, n_procs: str,
                  mem_num: str, mem_dim: str, qiime_env: str, chmod: str,
                  noloc: bool, jobs: bool, chunkit: int, tmp: str = None) -> str:
    """
    Write the main launcher of pbs scripts, written during using multiprocessing.

    :param job_folder: folder where the main job is to be written.
    :param analysis: current qqime2 analysis (e.g. PERMANOVA).
    :param all_sh_pbs: collection of all the sh scripts transformed to pbs.
    :param prjct_nm: Nick name for your project.
    :param time: walltime in hours.
    :param n_nodes: number of nodes to use.
    :param n_procs: number of processors to use.
    :param mem_num: memory in number.
    :param mem_dim: memory dimension to the number.
    :param qiime_env: qiime2-xxxx.xx conda environment.
    :param chmod: whether to change permission of output files (defalt: 775).
    :return: either the written launcher or nothing.
    """
    main_sh = '%s/%s.sh' % (job_folder, analysis)
    out_main_sh = ''
    warning = 0
    with open(main_sh, 'w') as main_o:
        chunks = {}
        if chunkit and len(all_sh_pbs) > chunkit:
            for idx, keys in enumerate(np.array_split(list(all_sh_pbs.keys()), chunkit)):
                head_sh = '%s/chunks/%s_chunk%s_%s.sh' % (job_folder, analysis, idx, prjct_nm)
                chunks[(idx, head_sh)] = [x for key in keys for x in all_sh_pbs[tuple(key)]]
        else:
            chunks = all_sh_pbs.copy()

        for (dat, out_sh), cur_shs in chunks.items():
            cur_written = False
            with open(out_sh, 'w') as sh:
                for cur_sh in cur_shs:
                    if isfile(cur_sh):
                        with open(cur_sh) as f:
                            for line in f:
                                sh.write(line)
                                cur_written = True
                        os.remove(cur_sh)
            if jobs:
                if cur_written:
                    out_pbs = '%s.pbs' % splitext(out_sh)[0]
                    run_xpbs(out_sh, out_pbs, '%s.%s' % (prjct_nm, dat), qiime_env,
                             time, n_nodes, n_procs, mem_num, mem_dim, chmod, 1,
                             '', None, noloc, jobs, tmp)
                    if os.getcwd().startswith('/panfs'):
                        out_pbs = out_pbs.replace(os.getcwd(), '')
                    main_o.write('qsub %s\n' % out_pbs)
                    out_main_sh = main_sh
                    warning += 1
                else:
                    os.remove(out_sh)
            else:
                if cur_written:
                    main_o.write('sh %s\n' % out_sh)
                    out_main_sh = main_sh
                    warning += 1
    if warning > 40:
        print(' -> [WARNING] >40 jobs here: please check before running!')
    return out_main_sh


def get_corresponding_meta(path: str) -> str:
    """
    Automatically gets the metadata file corresponding to the tsv / biom file.

    :param path: Path of the tsv / biom file.
    :return: metadata file path.
    """
    path_split = abspath(path).split('/')
    meta_rad = '%s/metadata/meta_%s' % ('/'.join(path_split[:-2]),
                                        '_'.join(splitext(basename(path))[0].split('_')[1:]))
    meta_tsv = '%s.tsv' % meta_rad
    meta_txt = '%s.txt' % meta_rad
    if isfile(meta_tsv):
        return meta_tsv
    elif isfile(meta_txt):
        return meta_txt
    else:
        print('No metadata found for %s\n(was looking for:\n- %s\n- %s)' % (path, meta_tsv, meta_txt))
        sys.exit(0)


def get_paths(i_datasets: tuple, i_datasets_folder: str) -> dict:
    """
    Collect the data and metadata files pairs per dataset.

    :param i_datasets: Internal name identifying the datasets in the input folder.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :return: datasets.
    """
    tsvs = []
    paths = {}
    for i_dataset in i_datasets:
        tsv = '%s/data/tab_%s.tsv' % (i_datasets_folder, i_dataset)
        biom = '%s.biom' % splitext(tsv)[0]
        tsvs.append(tsv)
        if isfile(tsv):
            paths[i_dataset] = tsv
        elif isfile(biom):
            paths[i_dataset] = biom
    if not paths:
        print('None of these target files found in input folder %s:' % i_datasets_folder)
        for tsv in tsvs:
            print(' - %s (or .biom)' % tsv)
        print('Exiting...')
        sys.exit(0)
    return paths


def gID_or_DNA(dat: str, path: str, path_pd: pd.DataFrame, datasets_read: dict,
               datasets_features: dict, datasets_phylo: dict) -> None:
    """
    Check whether the features of the current dataset are or contain:
    - genome IDs: then collect the gID -> corrected feature names for Web of Life tree shearing.
    - DNA sequences (e.g. typically deblur): then have a flag for sepp/phylo placement.
    (- to be developed for non-DNA OTU IDs associated with fasta sequences for sepp/phylo placement.)

    :param dat: name of the current dataset.
    :param path: feature table file path in the ./data folder.
    :param path_pd: feature table cotaining the features names in index.
    :param datasets_read: dataset -> [tsv table, meta table] (here it updates tsv table after features correction)
    :param datasets_features: to be updated with {gID: corrected feature name (no ';' or ' ')} per dataset.
    :param datasets_phylo: to be updated with ('tree_to_use', 'corrected_or_not') per dataset.
    """
    # regex to find the fist non-DNA character
    not_dna = re.compile('[^ACGTN].*?')
    if str(path_pd.index.dtype) == 'object':
        features_names = path_pd.index.tolist()
        # check if that's genome IDs
        found_gids = {}
        dna = True
        correction_needed = False
        for features_name in features_names:
            if dna and bool(not_dna.search(features_name)):
                dna = False
            if re.search('G\d{9}', features_name):
                if ';' in features_name:
                    correction_needed = True
                    features_name_corr = features_name.replace(';', '|').replace(' ', '')
                else:
                    features_name_corr = features_name
                found_gids[re.search('G\d{9}', features_name).group(0)] = features_name_corr
        if len(found_gids) == len(features_names):
            datasets_features[dat] = found_gids
            if correction_needed:
                path_pd.index = path_pd.index.str.replace(r'[; ]+', '|')
                path_pd.reset_index().to_csv(path, index=False, sep='\t')
                datasets_read[dat][0] = path_pd
                datasets_phylo[dat] = ('wol', 1)
            else:
                datasets_phylo[dat] = ('wol', 0)
        elif dna:
            datasets_phylo[dat] = ('amplicon', 0)


def get_datasets(i_datasets: tuple, i_datasets_folder: str) -> (dict, dict, dict, dict, dict):
    """
    Collect the pairs of features tables / metadata tables, formatted as in qiime2. e.g:

        --> Feature table example:
        #Feature ID  BVC.1591.10.10  BVC.1509.10.36  BVC.1584.10.10
        G000006785              0.0             0.0           175.0
        G000006925          34614.0          5973.0         12375.0
        G000007265              0.0           903.0           619.0

        --> Metadata table example:
        SampleID        age_years  age_wk40
        BVC.1591.10.10       0.75      0.79
        BVC.1509.10.36       3.00      0.77
        BVC.1584.10.10       0.75      0.77

    :param i_datasets: Internal name identifying the datasets in the input folder.
    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :return datasets infomration.
    """
    print('# Fetching data and metadata, for:\n - %s' % '\n - '.join(i_datasets))
    paths = get_paths(i_datasets, i_datasets_folder)

    datasets = {}
    datasets_read = {}
    datasets_phylo = {}
    datasets_features = {}
    datasets_rarefs = {}
    for dat, path in paths.items():
        meta = get_corresponding_meta(path)
        if not isfile(meta):
            print(meta, 'does not exist\n Skipping', dat)
            continue
        if path.endswith('.biom'):
            path_pd = load_table(path).to_dataframe()
        else:
            path_pd = read_meta_pd(path, '#OTU ID')
            #
            # tab_feat_col = get_feature_sample_col(path)
            # path_pd = pd.read_csv(path, header=0, index_col=0, sep='\t')
        meta_sam_col = get_feature_sample_col(meta)
        meta_pd = pd.read_csv(meta, header=0, sep='\t', dtype={meta_sam_col: str}, low_memory=False)
        meta_pd.rename(columns={meta_sam_col: 'sample_name'}, inplace=True)
        datasets[dat] = [[path, meta]]
        datasets_read[dat] = [[path_pd, meta_pd]]
        datasets_features[dat] = {}
        datasets_phylo[dat] = ('', 0)
        datasets_rarefs[dat] = ['']
        gID_or_DNA(dat, path, path_pd, datasets_read, datasets_features, datasets_phylo)
    return datasets, datasets_read, datasets_features, datasets_phylo, datasets_rarefs


def get_prjct_nm(project_name: str) -> str:
    """
    Get a smaller name for printing in qstat / squeue.

    :param project_name: Nick name for your project.
    :return: Shorter name (without vows).
    """
    alpha = 'aeiouy'
    prjct_nm = ''.join(x for x in project_name if x.lower() not in alpha)
    if prjct_nm == '':
        prjct_nm = 'q2.routine'
    return prjct_nm


def get_job_folder(i_datasets_folder: str, analysis: str):
    """
    Get the job folder name.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param analysis: name of the qiime2 analysis (e.g. beta).
    :return: job folder name.
    """

    job_folder = '%s/jobs/%s' % (i_datasets_folder, analysis)
    if not isdir(job_folder):
        os.makedirs(job_folder)
    return job_folder


def get_analysis_folder(i_datasets_folder, analysis):
    """
    Get the output folder name.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param analysis: name of the qiime2 analysis (e.g. beta).
    :return: output folder name.
    """
    odir = '%s/qiime/%s' % (i_datasets_folder, analysis)
    if not isdir(odir):
        os.makedirs(odir)
    return odir


def get_metrics(file_name: str, ABs: tuple) -> list:
    """
    Collect the alpha or beta diversity metrics from a resources file.

    :param file_name: name of the *_metrics file.
    :return: alpha or beta diversity metrics.
    """
    metrics = []
    with open('%s/%s.txt' % (RESOURCES, file_name)) as f:
        for line in f:
            line_strip = line.strip()
            if len(line_strip):
                if ABs:
                    if line_strip in ABs:
                        metrics.append(line_strip)
                else:
                    metrics.append(line_strip)
    return metrics


def get_wol_tree(i_wol_tree: str) -> str:
    """
    :param i_wol_tree: passed path to a tree.
    :return: path to a verified tree .nwk file.
    """
    if i_wol_tree == 'resources/wol_tree.nwk':
        return '%s/wol_tree.nwk' % RESOURCES
    if not isfile(i_wol_tree):
        print('%s does not exist\nExiting...' % i_wol_tree)
        sys.exit(0)
    if not i_wol_tree.endswith('.nwk'):
        if i_wol_tree.endswith('qza'):
            i_wol_tree_nwk = '%s.nwk' % splitext(i_wol_tree)[0]
            if isfile(i_wol_tree_nwk):
                print('Warning: about to overwrite %s\nExiting' % i_wol_tree_nwk)
                sys.exit(0)
            run_export(i_wol_tree, i_wol_tree_nwk, 'Phylogeny')
            return i_wol_tree_nwk
        else:
            # need more formal checks (sniff in skbio / stdout in "qiime tools peek")
            print('%s is not a .nwk (tree) file or not a qiime2 Phylogeny artefact\nExiting...' % i_wol_tree)
            sys.exit(0)
    else:
        return i_wol_tree


def get_sepp_tree(i_sepp_tree: str) -> str:
    """
    Get the full path of the reference database for SEPP.

    :param i_sepp_tree: database to use.
    :return: path of the reference database for SEPP.
    """
    if not i_sepp_tree or not isfile(i_sepp_tree):
        print('%s does not exist\nExiting...' % i_sepp_tree)
        sys.exit(0)
    if not i_sepp_tree.endswith('qza'):
        print('%s is not a qiime2 Phylogeny artefact\nExiting...' % i_sepp_tree)
        sys.exit(0)
    if basename(i_sepp_tree) in ['sepp-refs-silva-128.qza',
                                 'sepp-refs-gg-13-8.qza']:
        return i_sepp_tree
    else:
        print('%s is not:\n- "sepp-refs-silva-128.qza"\n- "sepp-refs-gg-13-8.qza"\n'
              'Download: https://docs.qiime2.org/2019.10/data-resources/#sepp-reference-databases)\n'
              'Exiting...' % i_sepp_tree)
        sys.exit(0)


def get_taxonomy_classifier(i_classifier: str) -> str:
    """
    Get the full path of the reference taxonomic classifier.

    :param i_classifier: database to use.
    :return: path of the reference taxonomy classifier.
    """
    if not isfile(i_classifier):
        print('%s does not exist\nExiting...' % i_classifier)
        sys.exit(0)
    if not i_classifier.endswith('qza'):
        print('%s is not a qiime2 artefact\nExiting...' % i_classifier)
        sys.exit(0)

    allowed_classifiers = [
        'silva-132-99-nb-classifier.qza',
        'silva-132-99-515-806-nb-classifier.qza',
        'gg-13-8-99-nb-classifier.qza',
        'gg-13-8-99-515-806-nb-classifier.qza'
    ]
    if basename(i_classifier) in allowed_classifiers:
        return i_classifier
    else:
        print('%s is not:\n%s'
              'Download: https://docs.qiime2.org/2020.2/data-resources/'
              '#taxonomy-classifiers-for-use-with-q2-feature-classifier)\n'
              'Exiting...' % (i_classifier, '\n - %s' % '\n - '.join(allowed_classifiers)))
        sys.exit(0)


def parse_g2lineage() -> dict:
    """
    :return: gOTU ID-to-gOTU dict.
    """
    g2lineage_fp = '%s/g2lineage.txt' % RESOURCES
    g2lineage = {}
    for line in open(g2lineage_fp).readlines():
        line_split = line.strip().split('\t')
        g2lineage[line_split[0]] = line_split[1]
    return g2lineage


def get_run_params(p_run_params: str) -> dict:

    run_params_default_fp = '%s/run_params.yml' % RESOURCES
    with open(run_params_default_fp) as handle:
        try:
            run_params_default = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            run_params_default = yaml.load(handle)

    if p_run_params and isfile(p_run_params):
        run_params_fp = p_run_params
        with open(run_params_fp) as handle:
            try:
                run_params = yaml.load(handle, Loader=yaml.FullLoader)
            except AttributeError:
                run_params = yaml.load(handle)
        run_params_default.update(run_params)
    else:
        print('using run parameters from', run_params_default_fp)

    return run_params_default


def filter_mb_table(preval: str, abund: str,
                    tsv_pd: pd.DataFrame,
                    do_res: bool=False) -> (pd.DataFrame, list):
    preval = float(preval)
    abund = float(abund)
    if abund:
        new_cols = {}
        cur_index = tsv_pd.index.tolist()
        cur_columns = tsv_pd.columns.tolist()
        for col in cur_columns:
            cur_col = tsv_pd[col]
            min_thresh = min([x for x in cur_col if x > 0]) * abund
            new_col = [x if x > min_thresh else 0 for x in cur_col]
            new_cols[col] = new_col
        tsv_pd = pd.DataFrame(new_cols, index=cur_index, columns = cur_columns)
        tsv_pd = tsv_pd[tsv_pd.sum(1) > 1]
    if preval:
        if preval < 1:
            n_perc = tsv_pd.shape[1] * preval
        else:
            n_perc = preval
        tsv_pd = tsv_pd.loc[tsv_pd.astype(bool).sum(1) >= n_perc, :]
    tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
    res = []
    if do_res:
        res = [preval, abund, tsv_pd.shape[0], tsv_pd.shape[1]]
    return tsv_pd, res


def filter_non_mb_table(preval: str, abund: str,
                        tsv_pd: pd.DataFrame,
                        do_res: bool=False) -> (pd.DataFrame, list):
    res = []
    preval = float(preval)
    abund = float(abund)
    if preval + abund == 0:
        if do_res:
            res = [0, 0, tsv_pd.shape[0], tsv_pd.shape[1]]
        return tsv_pd, res
    tsv_filt_pd = tsv_pd.copy()
    # get the min number of samples based on prevalence percent
    if preval < 1:
        n_perc = tsv_pd.shape[1] * preval
    else:
        n_perc = preval
    # abundance filter in terms of min reads counts
    tsv_pd_perc = tsv_filt_pd.copy()
    tsv_pd_perc_sum = tsv_filt_pd.sum(1)
    if abund < 1:
        tsv_pd_perc = tsv_filt_pd / tsv_filt_pd.sum()
        tsv_pd_perc_sum = tsv_filt_pd.sum(1) / tsv_filt_pd.sum(1).sum()
    abund_mode = 'sample'
    # remove features from feature table that are not present
    # in enough samples with the minimum number/percent of reads in these samples
    if abund_mode == 'sample':
        tsv_filt_pd = tsv_filt_pd.loc[(tsv_pd_perc > abund).sum(1) > n_perc, :]
    elif abund_mode == 'dataset':
        tsv_filt_pd = tsv_filt_pd.loc[tsv_pd_perc_sum > abund, :]
    elif abund_mode == 'both':
        tsv_filt_pd = tsv_filt_pd.loc[(tsv_pd_perc > abund).sum(1) > n_perc, :]
        fil_pd_perc_sum = tsv_filt_pd.sum(1)
        if abund < 1:
            fil_pd_perc_sum = tsv_filt_pd.sum(1) / tsv_filt_pd.sum(1).sum()
        tsv_filt_pd = tsv_filt_pd.loc[fil_pd_perc_sum > abund, :]
    else:
        raise Exception('"%s" mode not recognized' % abund_mode)
    tsv_filt_pd = tsv_filt_pd.loc[tsv_filt_pd.sum(1) > 0, tsv_filt_pd.sum(0) > 0]
    if do_res:
        res = [preval, abund, tsv_filt_pd.shape[0], tsv_filt_pd.shape[1]]
    return tsv_filt_pd, res


def get_meta_alpha(raref_dir, dat_rt, raref):
    meta_rgx = '%s/meta_%s%s*_alphas_full.tsv' % (raref_dir, dat_rt, raref)
    meta = glob.glob(meta_rgx)
    if not len(meta):
        meta_rgx = '%s/meta_%s%s*_alphas.tsv' % (raref_dir, dat_rt, raref)
        meta = glob.glob(meta_rgx)
        if not len(meta):
            meta_rgx = '%s/meta_%s%s*.tsv' % (raref_dir, dat_rt, raref)
            meta = glob.glob(meta_rgx)
            if not len(meta):
                meta = ''
            else:
                meta = sorted(meta)[0]
        else:
            meta = sorted(meta)[0]
    else:
        meta = sorted(meta)[0]
    return meta


def get_raref_table(dat_rt: str, raref: str, i_datasets_folder: str,
                    analysis: str) -> (pd.DataFrame, pd.DataFrame, str):
    raref_dir = get_analysis_folder(i_datasets_folder, 'rarefy/%s' % dat_rt)
    tsv_rgx = '%s/tab_%s%s*.tsv' % (raref_dir, dat_rt, raref)
    tsv_globbed = glob.glob(tsv_rgx)
    if len(tsv_globbed) >= 1:
        tsv = sorted(tsv_globbed)[0]
    else:
        print('Must have one rarefaction for "%s" to use it further in %s...'
              '\nExiting' % (dat_rt, analysis))
        sys.exit(0)
    meta = get_meta_alpha(raref_dir, dat_rt, raref)
    if not meta:
        return pd.DataFrame(), pd.DataFrame(), ''
    tsv_pd_, meta_pd_ = get_raref_tab_meta_pds(meta, tsv)
    return tsv_pd_, meta_pd_, meta


def write_filtered_tsv(tsv_out: str, tsv_pd: pd.DataFrame) -> None:
    tsv_sams_col = tsv_pd.reset_index().columns[0]
    tsv_pd = tsv_pd.reset_index().rename(columns={tsv_sams_col: 'Feature ID'}).set_index('Feature ID')
    tsv_pd.reset_index().to_csv(tsv_out, index=False, sep='\t')


def write_filtered_meta(meta_out: str, meta_pd_: pd.DataFrame, tsv_pd: pd.DataFrame) -> pd.DataFrame:
    meta_filt_pd = meta_pd_.loc[meta_pd_.sample_name.isin(tsv_pd.columns),:].copy()
    meta_filt_pd.to_csv(meta_out, index=False, sep='\t')
    return meta_filt_pd


def get_datasets_filtered(
        i_datasets_folder: str, datasets: dict,
        datasets_read: dict, datasets_filt: dict,
        unique_datasets: list, filtering: dict, force: bool,
        analysis: str, filt_datasets_done: dict,
        input_to_filtered: dict,
        already_computed: dict, subsets: dict) -> (dict, dict, list):
    """
    Filter the datasets for use in mmvec.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of data_sets.
    :param datasets_read: dataset -> [tsv table, meta table] (here it updates tsv table after features correction)
    :param datasets: datasets names from the yaml pairs.
    :param filtering: validated filtering thersholds.
    :param force: Force the re-writing of scripts for all commands.
    :return: list of datasets from filtered threshold.
    """

    drop_keys = {}
    filt_jobs = []
    filt_datasets = {}
    for (dat_, mb) in unique_datasets:

        if dat_ in datasets_filt:
            dat = datasets_filt[dat_]
        else:
            dat = dat_
        if dat not in datasets:
            if '__raref' in dat:
                split = dat.split('__raref')
                dat = '__raref'.join(split[:-1])
                raref = '_raref%s' % '__raref'.join(split[-1:])
                if dat in datasets_filt:
                    dat = datasets_filt[dat]
                tsv_pd_, meta_pd_, meta = get_raref_table(dat, raref, i_datasets_folder, analysis)
                if not tsv_pd_.shape[0]:
                    continue
                dat = '%s_%s' % (dat, raref)
                input_to_filtered[dat_] = dat
            else:
                print('%s dataset "%s" not found...' % (analysis, dat))
                continue

        elif not isinstance(datasets_read[dat][0], pd.DataFrame) and datasets_read[dat][0] == 'raref':
            tsv, meta = datasets[dat]
            if not isfile(tsv):
                print(analysis, 'Must have run rarefaction to use it further...\nExiting')
                sys.exit(0)
            tsv_pd_, meta_pd_ = get_raref_tab_meta_pds(meta, tsv)
            datasets_read[dat] = [tsv_pd_, meta_pd_]
            input_to_filtered[dat_] = dat
        else:
            tsv_pd_, meta_pd_ = datasets_read[dat][0]
            tsv, meta = datasets[dat][0]
            meta_alphas = get_meta_alpha(dirname(meta), dat, '')
            if meta_alphas and meta_alphas != meta:
                meta_pd_ = read_meta_pd(meta_alphas)
            input_to_filtered[dat_] = dat

        tsv_pd_ = tsv_pd_.loc[tsv_pd_.sum(1) > 0, :]
        tsv_pd_ = tsv_pd_.loc[:, tsv_pd_.sum(0) > 0]
        dat_filts = {}
        cases_dict = check_metadata_cases_dict(meta, meta_pd_, dict(subsets), 'songbird')
        for case_var, case_vals_list in cases_dict.items():
            for case_vals in case_vals_list:
                case = get_case(case_vals, case_var)
                case_meta_pd = get_new_meta_pd(meta_pd_, case, case_var, case_vals)
                case_tsv_pd = tsv_pd_[case_meta_pd.sample_name.tolist()]
                dat_dir = get_analysis_folder(i_datasets_folder, '%s/datasets/%s/%s' % (analysis, dat, case))
                prevals_abunds = filtering[(dat_, mb)]
                for (preval_abund, preval, abund) in sorted(prevals_abunds):
                    # make sure there's no empty row / column
                    if mb:
                        tsv_pd, res = filter_mb_table(preval, abund, case_tsv_pd)
                    else:
                        tsv_pd, res = filter_non_mb_table(preval, abund, case_tsv_pd)
                    rad_out = '%s_%s_%ss' % (dat, preval_abund, tsv_pd.shape[1])
                    tsv_out = '%s/tab_%s.tsv' % (dat_dir, rad_out)
                    tsv_qza = '%s.qza' % splitext(tsv_out)[0]
                    meta_out = '%s/meta_%s.tsv' % (dat_dir, rad_out)
                    tsv_hash = hash_pandas_object(tsv_pd).sum()
                    if len(filt_datasets_done[(dat, mb)][(case, preval_abund)]):
                        print('\t\t\t*', '[DONE]', dat, mb, case, preval_abund)
                        dat_filts[(case, preval_abund)] = filt_datasets_done[(dat, mb)][(case, preval_abund)]
                        meta_pd = write_filtered_meta(meta_out, case_meta_pd, tsv_pd)
                    else:
                        if analysis == 'songbird':
                            meta_out_mmvec = meta_out.replace('/songbird/', '/mmvec/')
                            tsv_out_mmvec = tsv_out.replace('/songbird/', '/mmvec/')
                            tsv_qza_mmvec = tsv_qza.replace('/songbird/', '/mmvec/')
                            if isfile(meta_out_mmvec):
                                meta_out = meta_out_mmvec
                                with open(meta_out) as f:
                                    for line in f:
                                        break
                                meta_pd = pd.read_csv(meta_out, header=0, sep='\t',
                                                      dtype={line.split('\t')[0]: str},
                                                      low_memory=False)
                            else:
                                meta_pd = write_filtered_meta(meta_out, case_meta_pd, tsv_pd)

                            if tsv_hash in already_computed:
                                drop_keys.setdefault((dat_, mb), []).append((preval_abund, preval, abund))
                                already_computed[tsv_hash].append([tsv_out, tsv_qza, meta_out])
                                tsv_out_src = already_computed[tsv_hash][0][0]
                                # tsv_qza_src = already_computed[tsv_hash][0][1]
                                # meta_out_src = already_computed[tsv_hash][0][2]
                                print('\t\t\t  [%s] Skip "%s" (same input as "%s")' % (
                                    analysis, basename(tsv_out), basename(tsv_out_src)))
                                continue
                                # already_computed[tsv_hash].append([tsv_out, tsv_qza, meta_out])
                                # tsv_out_src = already_computed[tsv_hash][0][0]
                                # tsv_qza_src = already_computed[tsv_hash][0][1]
                                # meta_out_src = already_computed[tsv_hash][0][2]
                                # if isfile(tsv_out):
                                #     cmd = '\nrm %s\nln -s %s %s\n' % (tsv_out, tsv_out_src, tsv_out)
                                #     filt_jobs.append(cmd)
                                # # if isfile(tsv_qza):
                                # #     cmd = '\nrm %s\nln -s %s %s\n' % (tsv_qza, tsv_qza_src, tsv_qza)
                                # #     filt_jobs.append(cmd)
                                # if isfile(meta_out):
                                #     cmd = '\nrm %s\nln -s %s %s\n' % (meta_out, meta_out_src, meta_out)
                                #     filt_jobs.append(cmd)
                            else:
                                if isfile(tsv_out_mmvec):
                                    # print(analysis, 'is file: tsv_out_mmvec', tsv_out_mmvec)
                                    # print('USE MMVECs tsv')
                                    # print(' - - -', tsv_out_mmvec)
                                    tsv_out = tsv_out_mmvec
                                elif force or not isfile(tsv_out):
                                    # print(analysis, 'write: tsv_out', tsv_out)
                                    write_filtered_tsv(tsv_out, tsv_pd)

                                if isfile(tsv_qza_mmvec):
                                    # print('USE MMVECs qza')
                                    # print(' - - -', tsv_qza_mmvec)
                                    # print(analysis, 'is file: tsv_qza_mmvec', tsv_qza_mmvec)
                                    tsv_qza = tsv_qza_mmvec
                                elif force or not isfile(tsv_qza):
                                    cmd = run_import(tsv_out, tsv_qza, 'FeatureTable[Frequency]')
                                    filt_jobs.append(cmd)
                                    # print(analysis, 'write (job): tsv_qza', tsv_qza)
                                already_computed[tsv_hash] = [[tsv_out, tsv_qza, meta_out]]
                        else:
                            meta_pd = write_filtered_meta(meta_out, case_meta_pd, tsv_pd)
                            if tsv_hash in already_computed:
                                drop_keys.setdefault((dat_, mb), []).append((preval_abund, preval, abund))
                                already_computed[tsv_hash].append([tsv_out, tsv_qza, meta_out])
                                tsv_out_src = already_computed[tsv_hash][0][0]
                                # tsv_qza_src = already_computed[tsv_hash][0][1]
                                # meta_out_src = already_computed[tsv_hash][0][2]
                                print('\t\t\t  [%s] Skip "%s" (same input as "%s")' % (
                                    analysis, basename(tsv_out), basename(tsv_out_src)))
                                continue
                                # if isfile(tsv_out):
                                #     cmd = '\nrm %s\nln -s %s %s\n' % (tsv_out, tsv_out_src, tsv_out)
                                #     filt_jobs.append(cmd)
                                # # if isfile(tsv_qza):
                                # #     cmd = '\nrm %s\nln -s %s %s\n' % (tsv_qza, tsv_qza_src, tsv_qza)
                                # #     filt_jobs.append(cmd)
                                # if isfile(meta_out):
                                #     cmd = '\nrm %s\nln -s %s %s\n' % (meta_out, meta_out_src, meta_out)
                                #     filt_jobs.append(cmd)
                            else:
                                if force or not isfile(tsv_out):
                                    write_filtered_tsv(tsv_out, tsv_pd)
                                if force or not isfile(tsv_qza):
                                    cmd = run_import(tsv_out, tsv_qza, 'FeatureTable[Frequency]')
                                    filt_jobs.append(cmd)
                                already_computed[tsv_hash] = [[tsv_out, tsv_qza, meta_out]]
                        print('\t\t\t* [TODO]', dat, mb, case, preval_abund, ':', tsv_pd.shape)
                        dat_filts[(case, preval_abund)] = [
                            tsv_out, tsv_qza, meta_out, meta_pd, tsv_pd.columns.tolist()]
        filt_datasets[(dat, mb)] = dat_filts
    return filt_datasets, filt_jobs


def check_datasets_filtered(
        i_datasets_folder: str, datasets: dict, datasets_filt: dict,
        unique_datasets: list, unique_filterings: dict,
        analysis: str, input_to_filtered: dict, subsets: dict) -> (dict, dict, list):
    """
    Filter the datasets for use in mmvec.

    :param i_datasets_folder: Path to the folder containing the data/metadata subfolders.
    :param datasets: list of data_sets.
    :param datasets_read: dataset -> [tsv table, meta table] (here it updates tsv table after features correction)
    :param datasets: datasets names from the yaml pairs.
    :param filtering: validated filtering thersholds.
    :param force: Force the re-writing of scripts for all commands.
    :return: list of datasets from filtered threshold.
    """

    filt_datasets_pass = {}
    for (dat_, mb) in unique_datasets:
        if dat_ in datasets_filt:
            dat = datasets_filt[dat_]
        else:
            dat = dat_
        if dat not in datasets:
            if '__raref' in dat:
                split = dat_.split('__raref')
                dat = '__raref'.join(split[:-1])
                raref = '_raref%s' % '__raref'.join(split[-1:])
                if dat in datasets_filt:
                    dat = datasets_filt[dat]
                dat = '%s_%s' % (dat, raref)
            else:
                print('%s dataset "%s" not found...' % (analysis, dat))
                continue
        input_to_filtered[dat_] = dat

        dat_filts_pass = {}
        for case_var, case_vals_list in subsets.items():
            for case_vals in case_vals_list:
                case = get_case(case_vals, case_var)
                dat_dir = get_analysis_folder(i_datasets_folder, '%s/datasets/%s/%s' % (analysis, dat, case))

                prevals_abunds = unique_filterings[(dat_, mb)]
                for (preval_abund, preval, abund) in sorted(prevals_abunds, key=lambda x: (x[1], x[2])):
                    # make sure there's no empty row / column
                    rad_out = '%s_%s_*s' % (dat, preval_abund)
                    tsv_out = '%s/tab_%s.tsv' % (dat_dir, rad_out)
                    tsv_qza = '%s.qza' % splitext(tsv_out)[0]
                    meta_out = '%s/meta_%s.tsv' % (dat_dir, rad_out)

                    if analysis == 'songbird':
                        meta_outs = glob.glob(meta_out.replace('/songbird/', '/mmvec/'))
                        tsv_outs = glob.glob(tsv_out.replace('/songbird/', '/mmvec/'))
                        tsv_qzas = glob.glob(tsv_qza.replace('/songbird/', '/mmvec/'))
                    else:
                        meta_outs = glob.glob(meta_out)
                        tsv_outs = glob.glob(tsv_out)
                        tsv_qzas = glob.glob(tsv_qza)

                    if len(meta_outs) == 1 and len(tsv_outs) == 1 and len(tsv_qzas) == 1:
                        meta_out = meta_outs[0]
                        tsv_out = tsv_outs[0]
                        tsv_qza = tsv_qzas[0]
                        with open(meta_out) as f:
                            for line in f:
                                sample_name = line.split('\t')[0]
                                break
                        meta_pd = pd.read_csv(meta_out, header=0, sep='\t',
                                              dtype={sample_name: str},
                                              low_memory=False)
                        dat_filts_pass[(case, preval_abund)] = [
                        # dat_filts_pass[(preval, abund)] = [
                            tsv_out, tsv_qza, meta_out,
                            meta_pd, meta_pd[sample_name].tolist()]
                    else:
                        dat_filts_pass[(case, preval_abund)] = []
                        # dat_filts_pass[(preval, abund)] = []
        filt_datasets_pass[(dat, mb)] = dat_filts_pass
    return filt_datasets_pass


def get_procrustes_mantel_dicts(p_procrustes_mantel):
    if not isfile(p_procrustes_mantel):
        print('yaml file for procrustes pairs does not exist:\n%s\nExiting...' % p_procrustes_mantel)
        sys.exit(0)
    with open(p_procrustes_mantel) as handle:
        try:
            procrustes_mantel_dict = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            procrustes_mantel_dict = yaml.load(handle)

    if 'pairs' not in procrustes_mantel_dict:
        print('No datasets pairs specified in %s:\nExiting\n' % p_procrustes_mantel)
        sys.exit(0)

    procrustes_mantel_pairs = {}
    for pair, paired_datasets in procrustes_mantel_dict['pairs'].items():
        n_dats = len(paired_datasets)
        if n_dats != 2:
            print('Must be two datasets per mmvec pair (found %s in %s)\n'
                  'Exiting\n' % (n_dats, p_procrustes_mantel))
            sys.exit(0)
        procrustes_mantel_pairs[pair] = paired_datasets

    procrustes_mantel_subsets = {'ALL': [[]]}
    if 'subsets' in procrustes_mantel_dict:
        procrustes_mantel_subsets.update(procrustes_mantel_dict['subsets'])

    return procrustes_mantel_pairs, procrustes_mantel_subsets


def get_collapse_taxo(p_collapse_taxo):
    if not isfile(p_collapse_taxo):
        print('yaml file for taxonomic collapse does not exist:\n%s\nExiting...' % p_collapse_taxo)
        sys.exit(0)
    with open(p_collapse_taxo) as handle:
        try:
            collapse_taxo = yaml.load(handle, Loader=yaml.FullLoader)
        except AttributeError:
            collapse_taxo = yaml.load(handle)

    return collapse_taxo
