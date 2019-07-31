rm(list=ls())

dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-amplification-oncogene.csv", stringsAsFactors = FALSE)
head(dat, n=sum(dat$p_adj < 0.1, na.rm = T))$X

nNonNAs = sum(!is.na(dat$pvalue))
head(dat, n=sum(dat$pvalue < 1/nNonNAs, na.rm = T))$X



dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-mutation-oncogene.csv", stringsAsFactors = FALSE)
head(dat, n=sum(dat$p_adj < 0.1, na.rm = T))$X

nNonNAs = sum(!is.na(dat$pvalue))
head(dat, n=sum(dat$pvalue < 1/nNonNAs, na.rm = T))$X


dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-mutation-tumor_suppressor.csv", stringsAsFactors = FALSE)
head(dat, n=sum(dat$p_adj < 0.1, na.rm = T))$X

nNonNAs = sum(!is.na(dat$pvalue))
head(dat, n=sum(dat$pvalue < 1/nNonNAs, na.rm = T))$X


