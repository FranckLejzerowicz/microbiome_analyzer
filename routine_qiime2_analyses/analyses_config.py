# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from routine_qiime2_analyses._routine_q2_io_utils import (
    check_input, get_prjct_anlss_nm, get_conda_envs, get_run_params,
    get_filt_raref_suffix, read_yaml_file
)


class AnalysesConfig(object):
    """Collect the data associated with each dataset passed but the user
    """
    # def __init__(self, args, kwargs) -> None:
    #
    #     """Initialize the class instance with the dataset name"""
    #     self.i_datasets = args[0]
    #     self.i_datasets_folder = args[1]
    #     self._check_input()
    #
    #     self.prjct_nm = get_prjct_nm(args[2])
    #     self.qiime_env = args[3]
    #     self.conda_envs = self._get_conda_envs()
    #     self.run_params = self._get_run_params(kwargs['p_run_params'])
    #
    #     self.raref = kwargs['raref']
    #     self.p_filt_threshs = kwargs['p_filt_threshs']
    #     self.filt_raref = self._get_filt_raref_suffix()
    #
    #     self.jobs = kwargs['jobs']
    #
    #     self.p_longi_column = kwargs['p_longi_column']
    #     self.p_raref_depths = kwargs['p_raref_depths']
    #     self.eval_rarefs = kwargs['eval_rarefs']
    #     self.p_alpha_subsets = kwargs['p_alpha_subsets']
    #     self.p_beta_subsets = kwargs['p_beta_subsets']
    #     self.p_perm_tests = kwargs['p_perm_tests']
    #     self.p_perm_tests_min = kwargs['p_perm_tests_min']
    #     self.p_beta_groups = kwargs['p_beta_groups']
    #     self.p_nestedness_groups = kwargs['p_nestedness_groups']
    #     self.p_beta_type = kwargs['p_beta_type']
    #     self.p_procrustes = kwargs['p_procrustes']
    #     self.p_mantel = kwargs['p_mantel']
    #     self.p_distance_decay = kwargs['p_distance_decay']
    #     self.p_collapse_taxo = kwargs['p_collapse_taxo']
    #     self.p_train_test = kwargs['p_train_test']
    #     self.p_adonis_formulas = kwargs['p_adonis_formulas']
    #     self.p_doc_config = kwargs['p_doc_config']
    #     self.p_sourcetracking_config = kwargs['p_sourcetracking_config']
    #     self.p_phate_config = kwargs['p_phate_config']
    #     self.do_biplots = kwargs['do_biplots']
    #     self.force = kwargs['force']
    #     self.i_classifier = kwargs['i_classifier']
    #     self.i_wol_tree = kwargs['i_wol_tree']
    #     self.i_sepp_tree = kwargs['i_sepp_tree']
    #     self.i_qemistree = kwargs['i_qemistree']
    #     self.p_diff_models = kwargs['p_diff_models']
    #     self.p_mmvec_pairs = kwargs['p_mmvec_pairs']
    #     self.p_mmvec_highlights = kwargs['p_mmvec_highlights']
    #     self.p_xmmvec = kwargs['p_xmmvec']
    #     self.p_chmod = kwargs['p_chmod']
    #     self.skip = kwargs['skip']
    #     self.gpu = kwargs['gpu']
    #     self.standalone = kwargs['standalone']
    #     self.loc = kwargs['loc']
    #     self.p_alphas = kwargs['p_alphas']
    #     self.p_betas = kwargs['p_betas']
    #     self.split = kwargs['split']
    #     self.dropout = kwargs['dropout']
    #     self.doc_phate = kwargs['doc_phate']
    #     self.filt3d = kwargs['filt3d']
    #     self.p_filt3d_config = kwargs['p_filt3d_config']
    #     self.filt_only = kwargs['filt_only']
    #     self.p_chunkit = kwargs['p_chunkit']

    def __init__(self, *args) -> None:

        """Initialize the class instance with the dataset name"""
        self.i_datasets = args[0]
        self.i_datasets_folder = args[1]
        self.project_name = args[2]
        self.qiime_env = args[3]

        self.raref = args[4]
        self.filt_threshs = args[5]
        self.longi_column = args[6]
        self.raref_depths = args[7]
        self.eval_rarefs = args[8]
        self.alpha_subsets = args[9]
        self.beta_subsets = args[10]
        self.perm_tests = args[11]
        self.perm_tests_min = args[12]
        self.beta_groups = args[13]
        self.nestedness_groups = args[14]
        self.beta_type = args[15]
        self.procrustes = args[16]
        self.mantel = args[17]
        self.distance_decay = args[18]
        self.collapse_taxo = args[19]
        self.train_test = args[20]
        self.train_test_dict = {}
        self.adonis_formulas = args[21]
        self.doc_config = args[22]
        self.sourcetracking_config = args[23]
        self.phate_config = args[24]
        self.do_biplots = args[25]
        self.force = args[26]
        self.i_classifier = args[27]
        self.i_wol_tree = args[28]
        self.i_sepp_tree = args[29]
        self.i_qemistree = args[30]
        self.diff_models = args[31]
        self.mmvec_pairs = args[32]
        self.mmvec_highlights = args[33]
        self.xmmvec = args[34]
        self.run_params = args[35]
        self.chmod = args[36]
        self.skip = args[37]
        self.gpu = args[38]
        self.standalone = args[39]
        self.loc = args[40]
        self.lphas = args[41]
        self.betas = args[42]
        self.split = args[43]
        self.dropout = args[44]
        self.doc_phate = args[45]
        self.filt3d = args[46]
        self.filt3d_config = args[47]
        self.filt_only = args[48]
        self.jobs = args[49]
        self.chunkit = args[50]

    def get_prjct_anlss_nm(self):
        return get_prjct_anlss_nm(self.project_name)

    def get_conda_envs(self):
        return get_conda_envs(self.qiime_env)

    def get_run_params(self):
        conda_envs = get_conda_envs(self.qiime_env)
        return get_run_params(self.run_params, conda_envs)

    def get_filt_raref_suffix(self):
        return get_filt_raref_suffix(self.filt_threshs, self.raref)

    def get_train_test_dict(self):
        self.train_test_dict = read_yaml_file(self.train_test)
        if 'train' not in self.train_test_dict:
            self.train_test_dict['train'] = 0.7
        elif float(self.train_test_dict['train']) < 0:
            self.train_test_dict['train'] = 0.7
        elif float(self.train_test_dict['train']) > 1:
            self.train_test_dict['train'] = 0.7

    def get_workflow(self):
        workflow = [('Import',),]
        if self.filt_threshs:
            workflow.append(('Filter',))
        if self.raref:
            workflow.append(('Rarefy',))
        if self.filt3d:
            workflow.append(('Filter 3D',))
            return workflow

        if self.i_qemistree and 'qemistree' not in self.skip:
            workflow.append(('Qemistree',))
        if 'taxonomy' not in self.skip:
            workflow.append(('Taxonomy',))
            if 'barplot' not in self.skip:
                workflow.append(('edit taxonomy',))
            workflow.append(('Edit taxonomy',))
        if 'wol' not in self.skip:
            workflow.append(('Shear WOL tree',))
        if self.i_sepp_tree and 'sepp' not in self.skip:
            workflow.append(('SEPP reads placement',))
        if self.filt_only:
            workflow.append(('Delete non-filtered',))
        if 'do_pies' in self.skip:
            workflow.append(('Make taxonomy pie charts',))
        if self.collapse_taxo and 'collapse' not in self.skip:
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
            if self.do_biplots and 'biplot' not in self.skip:
                workflow.append(('Principal coordinate analysis (biplot)',))
                if 'emperor_biplot' not in self.skip:
                    workflow.append(('EMPeror biplot',))
                if 'empress_biplot' not in self.skip:
                    workflow.append(('EMPress biplot',))
            if 'deicode' not in self.skip:
                workflow.append(('DEICODE',))
            if self.perm_tests and 'permanova' not in self.skip:
                workflow.append(('PERMANOVA',))
            if self.adonis_formulas and 'adonis' not in self.skip:
                workflow.append(('Adonis',))
            if self.procrustes and 'procrustes' not in self.skip:
                workflow.append(('Procrustes analysis',))
            if self.mantel and 'mantel' not in self.skip:
                workflow.append(('Mantel tests',))
            if self.nestedness_groups and 'nestedness' not in self.skip:
                workflow.append(('Nestedness analysis',))
            if self.distance_decay and 'decay' not in self.skip:
                workflow.append(('Beta-distance decay',))

        if self.phate_config and 'phate' not in self.skip:
            workflow.append(('phate',))
        if 'doc' not in self.skip and self.doc_config:
            workflow.append(('Dissimilarity-Overlap Curves',))
        if self.mmvec_pairs and 'mmvec' not in self.skip:
            workflow.append(('MMVEC',))
            if 'mmbird' not in self.skip:
                workflow.append(('MMBird',))
        if self.diff_models and 'songbird' not in self.skip:
            workflow.append(('Songbird',))
        return workflow
