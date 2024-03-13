import pandas as pd
from qiime2 import Artifact
from biom.table import Table

ranks_fp = 'RANKS_FP'
ranks_fp_tmp = 'RANKS_FP_TMP'
with open(ranks_fp_tmp, 'w') as o, open(ranks_fp) as f:
    for ldx, line in enumerate(f):
        ls = line.strip().split('\t')
        if ldx:
            if ls[0].isdigit():
                o.write('s%s\t%s\n' % (ls[0], '\t'.join(ls[1:])))
            else:
                o.write(line)
        else:
            ls_edit = ['s%s' % x if x.isdigit() else x for x in ls[1:]]
            o.write('%s\t%s\n' % (ls[0], '\t'.join(ls_edit)))

ranks_tmp_pd = pd.read_csv(ranks_fp_tmp, index_col=0, sep='\t')
ranks_tmp_qza = Artifact.import_data(
    "FeatureData[Conditional]", ranks_tmp_pd)
ranks_qza_tmp = 'RANKS_QZA_TMP'
ranks_tmp_qza.save(ranks_qza_tmp)

omic1_features = set(ranks_tmp_pd.columns)
omic2_features = set(ranks_tmp_pd.index)

omic1_common_fp = 'OMIC1_COMMON_FP'
omic1_common_fp_tmp = 'OMIC1_COMMON_FP_TMP'
with open(omic1_common_fp_tmp, 'w') as o, open(omic1_common_fp) as f:
    for ldx, line in enumerate(f):
        ls = line.strip().split('\t')
        if ldx:
            if ls[0] not in omic1_features:
                continue
            if ls[0].isdigit():
                o.write('s%s\t%s\n' % (ls[0], '\t'.join(ls[1:])))
            else:
                o.write(line)
        else:
            o.write('#OTU ID\t%s\n' % '\t'.join(['s%s' % x for x in ls[1:]]))

omic1_common_qza_tmp = 'OMIC1_COMMON_QZA_TMP'
omic1_common_fp_tmp_pd = pd.read_csv(omic1_common_fp_tmp, index_col=0, sep='\t')
omic1_common_fp_tmp_biom = Table(
    omic1_common_fp_tmp_pd.values,
    observation_ids=omic1_common_fp_tmp_pd.index,
    sample_ids=omic1_common_fp_tmp_pd.columns
)
omic1_common_fp_tmp_qza = Artifact.import_data(
    "FeatureTable[Frequency]", omic1_common_fp_tmp_biom)
omic1_common_fp_tmp_qza.save(omic1_common_qza_tmp)


omic2_common_fp = 'OMIC2_COMMON_FP'
omic2_common_fp_tmp = 'OMIC2_COMMON_FP_TMP'
with open(omic2_common_fp_tmp, 'w') as o, open(omic2_common_fp) as f:
    for ldx, line in enumerate(f):
        ls = line.strip().split('\t')
        if ldx:
            if ls[0] not in omic2_features:
                continue
            if ls[0].isdigit():
                o.write('s%s\t%s\n' % (ls[0], '\t'.join(ls[1:])))
            else:
                o.write(line)
        else:
            o.write('#OTU ID\t%s\n' % '\t'.join(['s%s' % x for x in ls[1:]]))

omic2_common_qza_tmp = 'OMIC2_COMMON_QZA_TMP'
omic2_common_fp_tmp_pd = pd.read_csv(omic2_common_fp_tmp, index_col=0, sep='\t')
omic2_common_fp_tmp_biom = Table(
    omic2_common_fp_tmp_pd.values,
    observation_ids=omic2_common_fp_tmp_pd.index,
    sample_ids=omic2_common_fp_tmp_pd.columns
)
omic2_common_fp_tmp_qza = Artifact.import_data(
    "FeatureTable[Frequency]", omic2_common_fp_tmp_biom)
omic2_common_fp_tmp_qza.save(omic2_common_qza_tmp)


taxonomy_tsv = 'TAXONOMY_TSV'
taxonomy_tsv_tmp = 'TAXONOMY_TSV_TMP'
taxonomy_tsv_pd = pd.read_csv(taxonomy_tsv, index_col=0, sep='\t')
taxonomy_tsv_pd.index = [
    's%s' % x if x.isdigit() or isinstance(x, float)
    else x for x in taxonomy_tsv_pd.index
]
taxonomy_tsv_pd = taxonomy_tsv_pd.reset_index()
taxonomy_tsv_pd = taxonomy_tsv_pd.rename(
    columns={taxonomy_tsv_pd.columns[0]: 'Feature ID'})
taxonomy_tsv_pd.to_csv(taxonomy_tsv_tmp, index=False, sep='\t')
