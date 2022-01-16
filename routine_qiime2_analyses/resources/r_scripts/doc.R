# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

library(DOC)
library(ggplot2)

otu <- read.table('<TAB_TSV>', header=T, sep='\t', comment.char='', check.names=F, nrows=2)
index_name <- colnames(otu)[1]
otu <- read.table('<TAB_TSV>', header=T, sep='\t', comment.char='', check.names=F, row.names=index_name)

if (file.exists("<DOC_TSV>")) {

    if (dim(otu)[1] > 100) {
        res <- DOC(otu)
        res.null <- DOC.null(otu)
        write.table(x=res$DO, file='<R_DIR>/DO.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res$LME, file='<R_DIR>/LME.tsv', sep='\t', quote=F, row.names=F)
        colnames(res$NEG) <- c('Neg_Slope', 'Data')
        write.table(x=res$NEG, file='<R_DIR>/NEG.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res$FNS, file='<R_DIR>/FNS.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res$BOOT, file='<R_DIR>/BOOT.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res$CI, file='<R_DIR>/CI.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res.null$DO, file='<R_DIR>/null_DO.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res.null$LME, file='<R_DIR>/null_LME.tsv', sep='\t', quote=F, row.names=F)
        colnames(res.null$NEG) <- c('Neg_Slope', 'Data')
        write.table(x=res.null$NEG, file='<R_DIR>/null_NEG.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res.null$FNS, file='<R_DIR>/null_FNS.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res.null$BOOT, file='<R_DIR>/null_BOOT.tsv', sep='\t', quote=F, row.names=F)
        write.table(x=res.null$CI, file='<R_DIR>/null_CI.tsv', sep='\t', quote=F, row.names=F)
    }

}
res = list(
    BOOT=read.table('<R_DIR>/BOOT.tsv', h=T, sep='\t'),
    CI=read.table('<R_DIR>/CI.tsv', h=T, sep='\t'),
    DO=read.table('<R_DIR>/DO.tsv', h=T, sep='\t'),
    LME=read.table('<R_DIR>/LME.tsv', h=T, sep='\t'),
    FNS=read.table('<R_DIR>/FNS.tsv', h=T, sep='\t'),
    NEG=read.table('<R_DIR>/NEG.tsv', h=T, sep='\t')
)
res.null = list(
    BOOT=read.table('<R_DIR>/null_BOOT.tsv', h=T, sep='\t'),
    CI=read.table('<R_DIR>/null_CI.tsv', h=T, sep='\t'),
    DO=read.table('<R_DIR>/null_DO.tsv', h=T, sep='\t'),
    LME=read.table('<R_DIR>/null_LME.tsv', h=T, sep='\t'),
    FNS=read.table('<R_DIR>/null_FNS.tsv', h=T, sep='\t'),
    NEG=read.table('<R_DIR>/null_NEG.tsv', h=T, sep='\t'))

colnames(res$NEG) <- c('Neg.Slope', 'Data')
colnames(res.null$NEG) <- c('Neg.Slope', 'Data')
res$DO <- res$DO[which(res$DO$Overlap <= 1),]
res.null$DO <- res.null$DO[which(res.null$DO$Overlap <= 1),]
pdf('<PDF>')
merged <- DOC.merge(list(s = res, s = res.null))
plot(merged)
dev.off()
