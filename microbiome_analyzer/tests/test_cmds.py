# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest
import numpy as np
import pkg_resources
from biom.table import Table
from os.path import isfile
from microbiome_analyzer.core.commands import (
    run_import, run_export, write_rarefy, write_fasta, write_taxonomy_sklearn,
)

OUT = pkg_resources.resource_filename("microbiome_analyzer", "tests")


class RunImports(unittest.TestCase):

    def setUp(self):

        self.tsv = 'tab.tsv'
        self.biom = 'tab.biom'
        self.qza = 'tab.qza'
        self.fas = 'tab.fasta'
        self.taxa_tsv = 'taxa.tsv'
        self.taxa_qza = 'taxa.qza'
        self.biplot_tsv = 'biplot.tsv'
        self.biplot_qza = 'biplot.qza'
        self.diff_tsv = 'diff.tsv'
        self.diff_qza = 'diff.qza'

        self.tab_cmd = 'qiime tools import \\\n'
        self.tab_cmd += '  --input-path %s \\\n' % self.biom
        self.tab_cmd += '  --output-path %s \\\n' % self.qza

        self.fas_cmd = 'qiime tools import \\\n'
        self.fas_cmd += '  --input-path %s \\\n' % self.fas
        self.fas_cmd += '  --output-path %s \\\n' % self.qza

        self.taxa_cmd = 'qiime tools import \\\n'
        self.taxa_cmd += '  --input-path %s \\\n' % self.taxa_tsv
        self.taxa_cmd += '  --output-path %s \\\n' % self.taxa_qza

        self.biplot_cmd = 'qiime tools import \\\n'
        self.biplot_cmd += '  --input-path %s \\\n' % self.biplot_tsv
        self.biplot_cmd += '  --output-path %s \\\n' % self.biplot_qza

        self.diff_cmd = 'qiime tools import \\\n'
        self.diff_cmd += '  --input-path %s \\\n' % self.diff_tsv
        self.diff_cmd += '  --output-path %s \\\n' % self.diff_qza

        self.biom_cmd = 'biom convert \\\n'
        self.biom_cmd += '  -i %s \\\n' % self.tsv
        self.biom_cmd += '  -o %s \\\n' % self.biom
        self.biom_cmd += '  --table-type="OTU table" \\\n'
        self.biom_cmd += '  --to-hdf5\n\n'

    def test_run_import(self):

        typ = 'FeatureTable[Frequency]'
        cmd = run_import(self.tsv, self.qza, typ)
        ref_cmd = self.biom_cmd + self.tab_cmd
        ref_cmd += '  --type "%s"\n' % typ
        self.assertEqual(ref_cmd, cmd)

        typ = 'FeatureData[Sequence]'
        cmd = run_import(self.fas, self.qza, typ)
        ref_cmd = self.fas_cmd + '  --type "%s"\n' % typ
        self.assertEqual(ref_cmd, cmd)

        typ = "FeatureData[Taxonomy]"
        cmd = run_import(self.taxa_tsv, self.taxa_qza, typ)
        ref_cmd = self.taxa_cmd + '  --type "%s"\n' % typ
        self.assertEqual(ref_cmd, cmd)

        typ = "FeatureData[Differential]"
        cmd = run_import(self.diff_tsv, self.diff_qza, typ)
        ref_cmd = self.diff_cmd + '  --type "%s"\n' % typ
        self.assertEqual(ref_cmd, cmd)

        typ = "PCoAResults % Properties('biplot')"
        cmd = run_import(self.biplot_tsv, self.biplot_qza, typ)
        ref_cmd = self.biplot_cmd + '  --type "%s"\n' % typ
        self.assertEqual(ref_cmd, cmd)


class RunExports(unittest.TestCase):

    def setUp(self):

        self.tsv = 'file.tsv'
        self.biom = 'file.biom'
        self.qza = 'file.qza'
        self.txt = 'file.txt'
        self.nwk = 'file.nwk'
        self.html = 'file.html'

        self.exp_cmd = 'qiime tools export \\\n'
        self.exp_cmd += '  --input-path %s \\\n' % self.qza
        self.exp_cmd += '  --output-path file\n'

        self.biom_cmd = 'mv file/*.biom file.biom\n'
        self.biom_cmd += 'rm -rf file\n'

        self.conv_cmd = 'mv file/*.biom file.biom\n'
        self.conv_cmd += 'biom convert'
        self.conv_cmd += '  -i file.biom \\\n'
        self.conv_cmd += '  -o file.tsv.tmp \\\n'
        self.conv_cmd += '  --to-tsv\n\n'
        self.conv_cmd += 'tail -n +2 file.tsv.tmp > file.tsv\n\n'
        self.conv_cmd += 'rm -rf file file.tsv.tmp\n'

        self.txt_cmd = 'mv file/*.txt file.txt\n'
        self.nwk_cmd = 'mv file/*.nwk file.nwk\n'
        self.tsv_cmd = 'mv file/*.tsv file.tsv\n'
        self.html_cmd = 'mv file/index.html file.html\n'

        self.rm = 'rm -rf file\n'

    def test_run_export(self):
        typ = 'FeatureTable'
        cmd = run_export(self.qza, self.tsv, typ)
        ref_cmd = self.exp_cmd + self.conv_cmd
        self.assertEqual(ref_cmd, cmd)
        cmd = run_export(self.qza, self.biom, typ)
        ref_cmd = self.exp_cmd + self.biom_cmd
        self.assertEqual(ref_cmd, cmd)

        for typ in ['pcoa', 'biplot', 'mmvec']:
            cmd = run_export(self.qza, self.txt, typ)
            ref_cmd = self.exp_cmd + self.txt_cmd + self.rm
            self.assertEqual(ref_cmd, cmd)

        for typ in ['perms', 'mantel', 'songbird', 'mmvec_summary']:
            cmd = run_export(self.qza, self.html, typ)
            ref_cmd = self.exp_cmd + self.html_cmd + self.rm
            self.assertEqual(ref_cmd, cmd)

        typ = 'phylogeny'
        cmd = run_export(self.qza, self.nwk, typ)
        ref_cmd = self.exp_cmd + self.nwk_cmd + self.rm
        self.assertEqual(ref_cmd, cmd)


class Writers(unittest.TestCase):

    def setUp(self):

        self.tsv = 'file.tsv'
        self.qza = 'in.qza'
        self.qza_out = 'out.qza'

        data = np.arange(15).reshape(5, 3)
        sample_ids = ['S%d' % i for i in range(3)]
        observ_ids = ['O%d' % i for i in range(5)]
        self.biom = Table(data, observ_ids, sample_ids, table_id='unittest')
        self.seqs_fas = '%s/seqs.fasta' % OUT
        self.seqs_qza = '%s/seqs.qza' % OUT

        self.cmd = 'qiime feature-table rarefy \\\n'
        self.cmd += '--i-table %s \\\n' % self.qza
        self.cmd += '--p-sampling-depth RAREF \\\n'
        self.cmd += '--o-rarefied-table %s\n' % self.qza_out

        self.rarefs = {1: self.cmd.replace('RAREF', '1'),
                       2: self.cmd.replace('RAREF', '2'),
                       3: self.cmd.replace('RAREF', '3')}

        self.cmd = '# Write features as fasta file:\n'
        self.cmd += '#  - Features from: file.tsv\n'
        self.cmd += '# Snippet:\n'
        self.cmd += '# ```:\n'
        self.cmd += "# with open(fasta_out, 'w') as o:\n"
        self.cmd += "#     for seq in tsv_pd.index:\n"
        self.cmd += "#         o.write('>%s\\n%s\\n' % (seq.strip(), seq.strip()))\n"
        self.cmd += '# ```:\n'
        self.cmd += run_import(self.seqs_fas, self.seqs_qza,
                               'FeatureData[Sequence]')

        self.sklearn_cmd = 'qiime feature-classifier classify-sklearn \\\n'
        self.sklearn_cmd += '--i-reads %s \\\n' % self.qza
        self.sklearn_cmd += '--i-classifier classifier.qza \\\n'
        self.sklearn_cmd += '--p-n-jobs 4 \\\n'
        self.sklearn_cmd += '--o-classification %s\n' % self.qza_out

    def test_write_rarefy(self):
        for raref, ref_cmd in self.rarefs.items():
            cmd = write_rarefy(self.qza, self.qza_out, raref)
            self.assertEqual(ref_cmd, cmd)

    def test_write_fasta(self):
        cmd = write_fasta(self.seqs_fas, self.seqs_qza, self.biom, self.tsv)
        self.assertEqual(True, isfile(self.seqs_fas))
        self.assertEqual(cmd, self.cmd)
        os.remove(self.seqs_fas)

    def test_write_taxonomy_sklearn(self):
        cmd = write_taxonomy_sklearn(self.qza_out, self.qza, 'classifier.qza')
        self.assertEqual(cmd, self.sklearn_cmd)


if __name__ == '__main__':
    unittest.main()
