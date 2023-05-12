# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

library(vegan)
library(readr)
library(data.table)

meta_fp <- 'META_FP'
meta_pd <- read.table(meta_fp, header = TRUE, check.names = FALSE, sep='\t', colClasses = c("sample_name"="character"))
row.names(meta_pd) <- meta_pd[,'sample_name']
meta_pd <- meta_pd[,-1]

METRIC_fp <- 'DM_FP'
METRIC <- read.table(METRIC_fp, check.names = FALSE)

res <- list()

perm <- how(nperm = NPERM)
NAME_meta_pd <- meta_pd[, c("VARS")]
NAME_meta_pd <- na.omit(NAME_meta_pd)
setBlocks(perm) <- with(NAME_meta_pd, BLOCK)

NAME_METRIC <- METRIC[row.names(NAME_meta_pd), row.names(NAME_meta_pd)]
res[['NAME_METRIC']] <- as.data.frame(adonis2(NAME_METRIC ~ FORMULA, data=NAME_meta_pdPERMUTATIONS))
res[['NAME_METRIC']]$strata <- "STRATA"
res[['NAME_METRIC']]$metric <- "METRIC"
res[['NAME_METRIC']]$factors <- row.names(res[['NAME_METRIC']])
row.names(res[['NAME_METRIC']]) <- NULL

write.table(x=rbindlist(res), file='OUT', quote=FALSE, sep='\t', row.names=FALSE)