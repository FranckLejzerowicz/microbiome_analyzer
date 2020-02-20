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
import pkg_resources
import pandas as pd
from os.path import basename, splitext, isfile, isdir, abspath

from routine_qiime2_analyses._routine_q2_xpbs import run_xpbs
from routine_qiime2_analyses._routine_q2_cmds import run_export

RESOURCES = pkg_resources.resource_filename("routine_qiime2_analyses", "resources")


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
            sys.exit(1)
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
        sys.exit(1)
    with open(p_formulas) as handle:
        formulas = yaml.load(handle, Loader=yaml.FullLoader)
    return formulas



def write_main_sh(job_folder: str, analysis: str, all_sh_pbs: dict,
                  prjct_nm: str, time: str, n_nodes: str, n_procs: str,
                  mem_num: str, mem_dim: str, qiime_env: str) -> str:
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
    :return: either the written launcher or nothing.
    """
    main_sh = '%s/%s.sh' % (job_folder, analysis)
    out_main_sh = ''
    with open(main_sh, 'w') as main_o:
        for (dat, out_sh), cur_shs in all_sh_pbs.items():
            cur_written = False
            with open(out_sh, 'w') as sh:
                for cur_sh in cur_shs:
                    if isfile(cur_sh):
                        sh.write('sh %s\n' % cur_sh)
                        cur_written = True
            if cur_written:
                out_pbs = '%s.pbs' % splitext(out_sh)[0]
                run_xpbs(out_sh, out_pbs, '%s.%s' % (prjct_nm, dat), qiime_env,
                         time, n_nodes, n_procs, mem_num, mem_dim, 1, '', None)
                main_o.write('qsub %s\n' % out_pbs)
                out_main_sh = main_sh
            else:
                os.remove(out_sh)
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
        sys.exit(1)


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
        sys.exit(1)
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


def get_datasets(i_datasets: tuple, i_datasets_folder: str) -> (dict, dict, dict, dict):
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
    print('# Fetching data and metadata (in %s)' % ', '.join(i_datasets))
    paths = get_paths(i_datasets, i_datasets_folder)

    datasets = {}
    datasets_read = {}
    datasets_phylo = {}
    datasets_features = {}
    for dat, path in paths.items():
        meta = get_corresponding_meta(path)
        if not isfile(meta):
            print(meta)
        path_pd = pd.read_csv(path, header=0, index_col=0, sep='\t')
        meta_pd = pd.read_csv(meta, header=0, sep='\t')
        datasets[dat] = [path, meta]
        datasets_read[dat] = [path_pd, meta_pd]
        datasets_features[dat] = {}
        datasets_phylo[dat] = ('', 0)
        gID_or_DNA(dat, path, path_pd, datasets_read, datasets_features, datasets_phylo)
    return datasets, datasets_read, datasets_features, datasets_phylo


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


def get_metrics(file_name: str) -> list:
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
        sys.exit(1)
    if not i_wol_tree.endswith('.nwk'):
        if i_wol_tree.endswith('qza'):
            i_wol_tree_nwk = '%s.nwk' % splitext(i_wol_tree)[0]
            if isfile(i_wol_tree_nwk):
                print('Warning: about to overwrite %s\nExiting' % i_wol_tree_nwk)
                sys.exit(1)
            run_export(i_wol_tree, i_wol_tree_nwk, 'Phylogeny')
            return i_wol_tree_nwk
        else:
            # need more formal checks (sniff in skbio / stdout in "qiime tools peek")
            print('%s is not a .nwk (tree) file or not a qiime2 Phylogeny artefact\nExiting...' % i_wol_tree)
            sys.exit(1)
    else:
        return i_wol_tree


def get_sepp_tree(i_sepp_tree: str) -> str:
    """
    Get the full path of the reference database for SEPP.

    :param i_sepp_tree: database to use.
    :return: path of the reference database for SEPP.
    """
    if not isfile(i_sepp_tree):
        print('%s does not exist\nExiting...' % i_sepp_tree)
        sys.exit(1)
    if not i_sepp_tree.endswith('qza'):
        print('%s is not a qiime2 Phylogeny artefact\nExiting...' % i_sepp_tree)
        sys.exit(1)
    if basename(i_sepp_tree) in ['sepp-refs-silva-128.qza',
                                 'sepp-refs-gg-13-8.qza']:
        return i_sepp_tree
    else:
        print('%s is not:\n- "sepp-refs-silva-128.qza"\n- "sepp-refs-gg-13-8.qza"\n'
              'Download: https://docs.qiime2.org/2019.10/data-resources/#sepp-reference-databases)\n'
              'Exiting...' % i_sepp_tree)
        sys.exit(1)
