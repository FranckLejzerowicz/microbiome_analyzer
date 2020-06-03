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
from os.path import basename, splitext, isfile, isdir, abspath

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_cmds import run_import, run_export
from routine_qiime2_analyses._routine_q2_metadata import rename_duplicate_columns

RESOURCES = pkg_resources.resource_filename("routine_qiime2_analyses", "resources")


def get_subsets(p_subsets: str) -> dict:
    """
    :param p_subsets: Subsets for alpha diversity.
    """
    if p_subsets:
        if not isfile(p_subsets):
            print('[Warning] yaml file for subsets does not exist: %s\n' % p_subsets)
        else:
            with open(p_subsets) as handle:
                subsets = yaml.load(handle, Loader=yaml.FullLoader)
            if not isinstance(subsets, dict):
                print('[Warning] %s must be a dictionary\n' % p_subsets)
            else:
                return subsets
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


def get_filtering(p_yml: str, filtering_dict: dict) -> dict:
    """
    Get the parameters for songbird passed by the user.
    :param p_mmvec_pairs: file containing the parameters.
    :param mmvec_dict: parsed content of the file containing the parameters.
    :return: parameters.
    """
    filtering = {
        'prevalence': [
            '0'
        ],
        'abundance': [
            ['0', '0']
        ]
    }
    if 'filtering' not in filtering_dict:
        print('No filtering thresholds set in %s:\nUsing defaults:' % p_yml)
        for k, v in filtering.items():
            print(k, v)
        return filtering
    else:
        filtering_dict.update(update_filtering_prevalence(filtering_dict, p_yml, filtering))
        filtering_dict.update(update_filtering_abundance(filtering_dict, p_yml, filtering))
        return filtering_dict['filtering']


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
            mmvec_dict, ', '.join(['%s: %s' % (k,v) for k,v in params.items()])))
    else:
        for param, cur_param in mmvec_dict['params'].items():
            if not isinstance(cur_param, list):
                print('Parameter %s should be a list (correct in %s)\n' % (param, p_mmvec_pairs))
                sys.exit(0)
            params[param] = cur_param
    return params


def get_mmvec_subsets(p_mmvec_pairs: str, mmvec_dict: dict) -> dict:
    mmvec_subsets = {}
    if 'subsets' in mmvec_dict:
        for pair, paired_datasets in mmvec_dict['subsets'].items():
            n_dats = len(paired_datasets)
            if n_dats != 2:
                print('Must be two datasets per mmvec pair (found %s in %s)\n'
                      'Exiting\n' % (n_dats, p_mmvec_pairs))
                sys.exit(0)
            mmvec_subsets[pair] = [(dat[:-1], 1) if dat[-1] == '*' else (dat, 0) for dat in paired_datasets]
    return mmvec_subsets


def get_mmvec_pairs(p_mmvec_pairs: str, mmvec_dict: dict, datasets_filt: dict) -> dict:
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
        for dat_ in paired_datasets:
            if dat_[-1] == '*':
                dat_, mb = dat_[:-1], 1
            else:
                mb = 0
            if dat_ in datasets_filt:
                dat = datasets_filt[dat_]
            else:
                dat = dat_
            paired.append((dat, mb))
        mmvec_pairs[pair] = paired
    return mmvec_pairs


def get_mmvec_dicts(p_mmvec_pairs: str, datasets_filt: dict) -> (dict, dict, dict):
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
        mmvec_dict = yaml.load(handle, Loader=yaml.FullLoader)

    mmvec_pairs = get_mmvec_pairs(p_mmvec_pairs, mmvec_dict, datasets_filt)
    mmvec_filtering = get_filtering(p_mmvec_pairs, mmvec_dict)
    mmvec_params = get_mmvec_params(p_mmvec_pairs, mmvec_dict)
    mmvec_subsets = get_mmvec_subsets(p_mmvec_pairs, mmvec_dict)
    return mmvec_pairs, mmvec_filtering, mmvec_params


def get_songbird_params(p_diff_models: str, diff_dict: dict) -> dict:
    """
    Get the parameters for songbird passed by the user.
    :param p_diff_models: file containing the parameters.
    :param diff_dict: parsed content of the file containing the parameters.
    :return: parameters.
    """
    params = {
        'batches': ['2'],
        'learns': ['1e-4'],
        'epochs': ['5000'],
        'thresh_feats': ['0'],
        'thresh_samples': ['0'],
        'diff_priors': ['0.1', '1']
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


def get_songbird_dicts(p_diff_models: str) -> (dict, dict, dict):
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
        diff_dict = yaml.load(handle, Loader=yaml.FullLoader)

    main_cases_dict = {'ALL': [[]]}
    if 'subsets' in diff_dict:
        main_cases_dict.update(diff_dict['subsets'])

    models = get_songbird_models(p_diff_models, diff_dict)
    params = get_songbird_params(p_diff_models, diff_dict)
    filtering = get_filtering(p_diff_models, diff_dict)
    datasets = [(dat[:-1], 1) if dat[-1] == '*' else (dat, 0) for dat in models]
    return models, filtering, params, datasets, main_cases_dict


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
            main_cases_dict.update(yaml.load(handle, Loader=yaml.FullLoader))
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
        formulas = yaml.load(handle, Loader=yaml.FullLoader)
    return formulas


def get_sample_col(meta: str) -> str:
    """
    Get the first column of the metadata file.

    :param meta: metadata file name
    :return: column name
    """
    n = 0
    with open(meta) as f:
        for line in f:
            n += 1
            break
    if n:
        return line.split()[0]
    else:
        print('Empty now: %s (possily being written elsewhere..)' % meta)
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
    meta_raref_pd = meta_pd.loc[meta_pd.sample_name.isin(tsv_pd.columns.tolist()),:].copy()
    meta_raref_pd.to_csv(meta, index=False, sep='\t')
    return tsv_pd, meta_raref_pd


def read_meta_pd(meta: str) -> pd.DataFrame:
    """
    Read metadata wit first column as index.
    :param meta: file path to the metadata file.
    :return: metadata table.
    """
    meta_sam_col = get_sample_col(meta)
    meta_pd = pd.read_csv(meta, header=0, sep='\t', dtype={meta_sam_col: str}, low_memory=False)
    meta_pd.rename(columns={meta_sam_col: 'sample_name'}, inplace=True)
    # meta_pd.set_index('sample_name', inplace=True)
    return meta_pd


def write_main_sh(job_folder: str, analysis: str, all_sh_pbs: dict,
                  prjct_nm: str, time: str, n_nodes: str, n_procs: str,
                  mem_num: str, mem_dim: str, qiime_env: str, chmod: str, noloc: bool) -> str:
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
        for (dat, out_sh), cur_shs in all_sh_pbs.items():
            cur_written = False
            with open(out_sh, 'w') as sh:
                for cur_sh in cur_shs:
                    if isfile(cur_sh):
                        with open(cur_sh) as f:
                            for line in f:
                                sh.write(line)
                                cur_written = True
                        os.remove(cur_sh)
            if cur_written:
                out_pbs = '%s.pbs' % splitext(out_sh)[0]
                run_xpbs(out_sh, out_pbs, '%s.%s' % (prjct_nm, dat), qiime_env,
                         time, n_nodes, n_procs, mem_num, mem_dim, chmod, 1, '', None, noloc)
                if os.getcwd().startswith('/panfs'):
                    out_pbs = out_pbs.replace(os.getcwd(), '')
                main_o.write('qsub %s\n' % out_pbs)
                out_main_sh = main_sh
                warning += 1
            else:
                os.remove(out_sh)
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
        path_pd = pd.read_csv(path, header=0, index_col=0, sep='\t')
        meta_sam_col = get_sample_col(meta)
        meta_pd = pd.read_csv(meta, header=0, sep='\t', dtype={meta_sam_col: str}, low_memory=False)
        meta_pd.rename(columns={meta_sam_col: 'sample_name'}, inplace=True)
        datasets[dat] = [path, meta]
        # path_pd : indexed with feature name
        # meta_pd : not indexed -> "sample_name" as first column
        datasets_read[dat] = [path_pd, meta_pd]
        datasets_features[dat] = {}
        datasets_phylo[dat] = ('', 0)
        datasets_rarefs[dat] = ''
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
    for line in open(g2lineage_fp).readlines()[1:]:
        line_split = line.strip().split('\t')
        g2lineage[line_split[0]] = line_split[1]
    return g2lineage


def get_run_params(p_run_params: str) -> dict:

    run_params_default_fp = '%s/run_params.yml' % RESOURCES
    with open(run_params_default_fp) as handle:
        run_params_default = yaml.load(handle, Loader=yaml.FullLoader)

    if p_run_params and isfile(p_run_params):
        run_params_fp = p_run_params
        with open(run_params_fp) as handle:
            run_params = yaml.load(handle, Loader=yaml.FullLoader)
        run_params_default.update(run_params)
    else:
        print('using run parameters from', run_params_default_fp)

    return run_params_default


def filter_mb_table(preval_filt: int, abund_filt: int,
                    tsv_pd: pd.DataFrame) -> pd.DataFrame:
    if preval_filt or abund_filt:
        new_cols = []
        cur_index = tsv_pd.index.tolist()
        cur_columns = tsv_pd.columns.tolist()
        for r, row in tsv_pd.iterrows():
            if sum(row == 0):
                min_thresh = min([x for x in row if x > 0]) * abund_filt
                cur_row = [x if x >= min_thresh else 0 for x in row]
                new_cols.append(cur_row)
            else:
                new_cols.append(row.tolist())
        tsv_pd = pd.DataFrame(
            new_cols,
            index=cur_index,
            columns=cur_columns
        )
        tsv_pd = tsv_pd[tsv_pd.sum(1) > 1]
        n_perc = (preval_filt / tsv_pd.shape[1]) * 100
        if preval_filt and abund_filt:
            tsv_pd = tsv_pd.loc[tsv_pd.astype(bool).sum(1) >= n_perc, :]
        tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
    return tsv_pd


def filter_non_mb_table(preval_filt: int, abund_filt: int,
                        tsv_pd: pd.DataFrame) -> pd.DataFrame:
    n_perc = (preval_filt / tsv_pd.shape[1]) * 100
    if preval_filt and abund_filt:
        tsv_filt_pd = tsv_pd.loc[(tsv_pd.values >= abund_filt).sum(1) >= n_perc, :].copy()
    elif preval_filt:
        tsv_filt_pd = tsv_pd.loc[tsv_pd.values.astype(bool).sum(1) >= n_perc, :].copy()
    elif abund_filt:
        tsv_filt_pd = tsv_pd.loc[tsv_pd.sum(1) > abund_filt, :].copy()
    else:
        tsv_filt_pd = tsv_pd.copy()
    tsv_filt_pd = tsv_filt_pd.loc[:, tsv_filt_pd.sum(0) > 0]
    return tsv_filt_pd


def get_raref_table(dat: str, i_datasets_folder: str,
                    analysis: str) -> (pd.DataFrame, pd.DataFrame):
    dat_rt = dat.split('__raref')[0]
    raref_dir = get_analysis_folder(i_datasets_folder, 'rarefy/%s' % dat_rt)
    tsv_globbed = glob.glob('%s/tab_%s_raref*.tsv' % (raref_dir, dat_rt))
    if len(tsv_globbed) != 1:
        print('Must have one rarefaction for "%s" to use it further in %s ("%s")...'
              '\nExiting' % (dat_rt, analysis, dat))
        sys.exit(0)
    tsv = tsv_globbed[0]
    meta = glob.glob('%s/meta_%s_raref*_alphas.tsv' % (raref_dir, dat_rt))
    if not len(meta):
        meta = glob.glob('%s/meta_%s_raref*.tsv' % (raref_dir, dat_rt))
        if not len(meta):
            return pd.DataFrame(), pd.DataFrame()
        else:
            meta = meta[0]
    else:
        meta = meta[0]
    tsv_pd_, meta_pd_ = get_raref_tab_meta_pds(meta, tsv)
    return tsv_pd_, meta_pd_


def write_filtered_tsv(tsv_out: str, tsv_pd: pd.DataFrame) -> None:
    tsv_sams_col = tsv_pd.reset_index().columns[0]
    tsv_pd = tsv_pd.reset_index().rename(columns={tsv_sams_col: 'Feature ID'}).set_index('Feature ID')
    tsv_pd.reset_index().to_csv(tsv_out, index=False, sep='\t')


def write_filtered_meta(meta_out: str, meta_pd_: pd.DataFrame, tsv_pd: pd.DataFrame) -> pd.DataFrame:
    meta_filt_pd = meta_pd_.loc[meta_pd_.sample_name.isin(tsv_pd.columns),:].copy()
    meta_filt_pd.to_csv(meta_out, index=False, sep='\t')
    return meta_filt_pd


def get_datasets_filtered(i_datasets_folder: str, datasets: dict,
                          datasets_read: dict, unique_datasets: list,
                          filtering: dict, force: bool, analysis: str) -> (dict, list):
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
    filt_jobs = []
    filt_datasets = {}
    for (dat, mb) in unique_datasets:
        if dat not in datasets:
            if dat.endswith('__raref'):
                tsv_pd_, meta_pd_ = get_raref_table(dat, i_datasets_folder, analysis)
                if not tsv_pd_.shape[0]:
                    continue
            else:
                print('%s dataset "%s" not found...' % (analysis, dat))
                continue
        elif datasets_read[dat] == 'raref':
            tsv, meta = datasets[dat]
            if not isfile(tsv):
                print('Must have run rarefaction to use it further...\nExiting')
                sys.exit(0)
            tsv_pd_, meta_pd_ = get_raref_tab_meta_pds(meta, tsv)
            datasets_read[dat] = [tsv_pd_, meta_pd_]
        else:
            tsv_pd_, meta_pd_ = datasets_read[dat]

        dat_filts = {}
        dat_dir = get_analysis_folder(i_datasets_folder, '%s/datasets/%s' % (analysis, dat))
        for preval_filt in filtering['prevalence']:
            for abund_filt in filtering['abundance']:
                # make sure there's no empty row / column
                tsv_pd = tsv_pd_.loc[tsv_pd_.sum(1) > 0, :].copy()
                tsv_pd = tsv_pd.loc[:, tsv_pd.sum(0) > 0]
                if mb:
                    abund_filter = int(abund_filt[1])
                    tsv_pd = filter_mb_table(int(preval_filt), abund_filter, tsv_pd)
                else:
                    abund_filter = int(abund_filt[0])
                    tsv_pd = filter_non_mb_table(int(preval_filt), abund_filter, tsv_pd)
                rad_out = '%s_%s_%s_%ss' % (dat, preval_filt, abund_filter, tsv_pd.shape[1])
                tsv_out = '%s/tab_%s.tsv' % (dat_dir, rad_out)
                tsv_qza = '%s.qza' % splitext(tsv_out)[0]
                meta_out = '%s/meta_%s.tsv' % (dat_dir, rad_out)

                if analysis == 'songbird':
                    meta_out_mmvec = meta_out.replace('/songbird/', '/mmvec/')
                    tsv_out_mmvec = tsv_out.replace('/songbird/', '/mmvec/')
                    tsv_qza_mmvec = tsv_qza.replace('/songbird/', '/mmvec/')
                    if isfile(meta_out_mmvec):
                        meta_out = meta_out_mmvec
                        meta_pd = pd.read_csv(meta_out, header=0, sep='\t',
                                              dtype={'sample_name': str},
                                              low_memory=False)
                    else:
                        meta_pd = write_filtered_meta(meta_out, meta_pd_, tsv_pd)

                    if isfile(tsv_out_mmvec):
                        tsv_out = tsv_out_mmvec
                    elif force or not isfile(tsv_out):
                        write_filtered_tsv(tsv_out, tsv_pd)

                    if isfile(tsv_qza_mmvec):
                        tsv_qza = tsv_qza_mmvec
                    elif force or not isfile(tsv_qza):
                        cmd = run_import(tsv_out, tsv_qza, 'FeatureTable[Frequency]')
                        filt_jobs.append(cmd)
                else:
                    meta_pd = write_filtered_meta(meta_out, meta_pd_, tsv_pd)
                    if force or not isfile(tsv_out):
                        write_filtered_tsv(tsv_out, tsv_pd)
                    if force or not isfile(tsv_qza):
                        cmd = run_import(tsv_out, tsv_qza, 'FeatureTable[Frequency]')
                        filt_jobs.append(cmd)
                dat_filts[(preval_filt, str(abund_filter))] = [
                    tsv_out, tsv_qza, meta_out, meta_pd, tsv_pd.columns.tolist()
                ]
        filt_datasets[dat] = dat_filts
    return filt_datasets, filt_jobs


def get_procrustes_dicts(p_procrustes):
    if not isfile(p_procrustes):
        print('yaml file for procrustes pairs does not exist:\n%s\nExiting...' % p_mmvec_pairs)
        sys.exit(0)
    with open(p_procrustes) as handle:
        procrustes_dict = yaml.load(handle, Loader=yaml.FullLoader)

    if 'pairs' not in procrustes_dict:
        print('No datasets pairs specified in %s:\nExiting\n' % p_procrustes)
        sys.exit(0)

    procrustes_pairs = {}
    for pair, paired_datasets in procrustes_dict['pairs'].items():
        n_dats = len(paired_datasets)
        if n_dats != 2:
            print('Must be two datasets per mmvec pair (found %s in %s)\n'
                  'Exiting\n' % (n_dats, p_procrustes))
            sys.exit(0)
        procrustes_pairs[pair] = paired_datasets

    procrustes_subsets = {'ALL': [[]]}
    if 'subsets' in procrustes_dict:
        procrustes_subsets.update(procrustes_dict['subsets'])

    return procrustes_pairs, procrustes_subsets
