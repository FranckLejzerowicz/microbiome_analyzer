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
import pandas as pd
from glob import glob
from os.path import splitext, basename, isfile, isdir


def run_import(input_path: str, output_path: str, typ: str) -> str:
    """
    Return the import qiime2 command.

    :param input_path: input file path.
    :param output_path: output file path.
    :param typ: qiime2 type.
    :return: command to qiime2.
    """
    cmd = ''
    if typ.startswith("FeatureTable"):
        if not input_path.endswith('biom'):
            cur_biom = '%s.biom' % splitext(input_path)[0]
            cmd += 'biom convert'
            cmd += ' -i %s -o %s' % (input_path, cur_biom)
            cmd += ' --table-type="OTU table" --to-hdf5\n'
            cmd += 'qiime tools import'
            cmd += ' --input-path %s --output-path %s' % (cur_biom, output_path)
            cmd += ' --type "FeatureTable[Frequency]"\n'
        else:
            cmd += 'qiime tools import'
            cmd += ' --input-path %s --output-path %s' % (input_path, output_path)
            cmd += ' --type "FeatureTable[Frequency]"\n'
    else:
        cmd += 'qiime tools import'
        cmd += ' --input-path %s --output-path %s' % (input_path, output_path)
        cmd += ' --type "%s"\n' % typ
    return cmd


def run_export(input_path: str, output_path: str, typ: str) -> str:
    """
    Return the export qiime2 command.

    :param input_path: input file path.
    :param output_path: output file path.
    :param typ: qiime2 type.
    :return: command to qiime2.
    """
    cmd = ''
    if typ.startswith("FeatureTable"):
        if not output_path.endswith('biom'):
            cur_biom = '%s.biom' % splitext(output_path)[0]
            cmd += 'qiime tools export --input-path %s --output-path %s\n' % (input_path, splitext(output_path)[0])
            cmd += 'mv %s/*.biom %s\n' % (splitext(output_path)[0], cur_biom)
            cmd += 'biom convert -i %s -o %s.tmp --to-tsv\n' % (cur_biom, output_path)
            cmd += 'tail -n +2 %s.tmp > %s\n' % (output_path, output_path)
            cmd += 'rm -rf %s %s.tmp\n' % (splitext(output_path)[0], output_path)
        else:
            cmd += 'qiime tools export --input-path %s --output-path %s\n' % (input_path, splitext(output_path)[0])
            cmd += 'mv %s/*.biom %s\n' % (splitext(input_path)[0], output_path)
            cmd += 'rm -rf %s\n' % splitext(input_path)[0]
    else:
        cmd += 'qiime tools export --input-path %s --output-path %s\n' % (input_path, splitext(output_path)[0])
        cmd += 'mv %s/*.tsv %s\n' % (splitext(output_path)[0], output_path)
        cmd += 'rm -rf %s\n' % splitext(output_path)[0]
    return cmd


def get_corresponding_meta(path):
    """
    Automatically gets the metadata file corresponding to the tsv / biom file.

    :param path: Path of the tsv / biom file.
    :return:
    """
    meta_rad = splitext(path.replace('/data/', '/metadata/').replace('/tab_', '/meta_'))[0]
    meta_files = glob('%s.*' % meta_rad)
    if len(meta_files) != 1:
        print('No metadata found for %s\n(was looking for %s)' % (path, meta_rad))
        sys.exit(1)
    meta = meta_files[0]
    return meta


def get_paths(i_datasets: tuple, i_folder: str) -> dict:
    """

    :param i_datasets: Internal name identifying the datasets in the input folder.
    :param i_folder: Path to the folder containing the .tsv datasets.
    :return:
    """
    tsvs = []
    paths = {}
    for i_dataset in i_datasets:
        tsv = '%s/data/tab_%s.tsv' % (i_folder, i_dataset)
        biom = '%s.biom' % splitext(tsv)[0]
        tsvs.append(tsv)
        if isfile(tsv):
            paths[i_dataset] = tsv
        elif isfile(biom):
            paths[i_dataset] = biom
    if not paths:
        print('None of these target files found in input folder %s:' % i_folder)
        for tsv in tsvs:
            print(' - %s (or .biom)' % tsv)
        print('Exiting...')
        sys.exit(1)
    return paths


def get_datasets(i_datasets: tuple, i_folder: str, gid: bool, biom: bool) -> (dict, dict, dict):
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
    :param i_folder: Path to the folder containing the .tsv datasets.
    :param gid: If feature names have the genome ID (to use the Web of Life tree).
    :param biom: Use biom files in the input folder
    :return
    """
    print('# Fetching data and metadata (in %s)' % ', '.join(i_datasets))
    paths = get_paths(i_datasets, i_folder)

    datasets = {}
    datasets_read = {}
    datasets_features = {}
    for dat, path in paths.items():
        meta = get_corresponding_meta(path)
        if not isfile(meta):
            print(meta)
        path_pd = pd.read_csv(path, header=0, index_col=0, sep='\t')
        meta_pd = pd.read_csv(meta, header=0, sep='\t')
        datasets.setdefault(dat, []).append([path, meta])
        datasets_read.setdefault(dat, []).append([path_pd, meta_pd])
        if gid and str(tab_filt_pd.index.dtype) == 'object':
            features_names = path_pd.index.tolist()
            found_gids = {}
            for features_name in features_names:
                if re.search('G\d{9}', features_name):
                    found_gids[re.search('G\d{9}', features_name).group(0)] = features_name
            if found_gids:
                datasets_features[dat] = found_gids

    return datasets, datasets_read, datasets_features


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


def get_job_folder(i_folder: str, analysis: str):
    """
    Get the job folder name.

    :param i_folder: Path to the folder containing the .tsv datasets.
    :param analysis: name of the qiime2 analysis (e.g. beta)
    :return: job folder name
    """
    job_folder = '%s/jobs/%s' % (i_folder, analysis)
    if not isdir(job_folder):
        os.makedirs(job_folder)
    return job_folder


def get_analysis_folder(i_folder, analysis):
    """
    Get the output folder name.

    :param i_folder: Path to the folder containing the .tsv datasets.
    :param analysis: name of the qiime2 analysis (e.g. beta)
    :return: output folder name
    """
    odir = '%s/qiime/%s' % (i_folder, analysis)
    if not isdir(odir):
        os.makedirs(odir)
    return odir

