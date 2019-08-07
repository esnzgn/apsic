rm(list=ls())

minNrOfCelllines = 383/100

dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-amplification-oncogene.csv", stringsAsFactors = FALSE)
n = sum(dat$X.mut > minNrOfCelllines)
head(dat, n=sum(dat$pvalue < min(1/n, 0.05), na.rm = T))$X


dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-mutation-oncogene.csv", stringsAsFactors = FALSE)
n = sum(dat$X.mut > minNrOfCelllines)
head(dat, n=sum(dat$pvalue < min(1/n, 0.05), na.rm = T))$X


dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-mutation-tumor_suppressor.csv", stringsAsFactors = FALSE)
n = sum(dat$X.mut > minNrOfCelllines)
head(dat, n=sum(dat$pvalue < min(1/n, 0.05), na.rm = T))$X
