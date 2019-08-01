rm(list=ls())

printNamesOnPlot <- function(geneNames, startX, startY, gapX, gapY, col, genesPerRow = 3, cex=0.7) {
   for(i in 1:length(geneNames) ) { 
     x = startX + ((i-1) %% genesPerRow) * gapX
     y = startY - (floor((i-1)/genesPerRow ) ) * gapY
     text(x, y  , geneNames[i], cex=cex, col=col, font=2)
   }
 }

datMis = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-mutation-oncogene.csv", row.names = 1)
datAmp = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-amplification-oncogene.csv", row.names = 1)
datTru = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-mutation-tumor_suppressor.csv", row.names = 1)

maxThr = 8

pMis = -log10(datMis$pvalue)
names(pMis) = rownames(datMis)
pMis=pMis[!is.na(pMis)]
pMis[pMis > maxThr] = maxThr

pAmp = -log10(datAmp$pvalue)
names(pAmp) = rownames(datAmp)
pAmp=pAmp[!is.na(pAmp)]
pAmp[pAmp > maxThr] = maxThr

pTru = -log10(datTru$pvalue)
names(pTru) = rownames(datTru)
pTru=pTru[!is.na(pTru)]
pTru[pTru > maxThr] = maxThr

allPvalues = c(pMis, pAmp, pTru)
allPvalues = allPvalues[which(!is.na(allPvalues))]

pdf("figures/pan-cancer-genes.pdf", width=8, height=6)
den =density(allPvalues) 
plot(den, axes=F, xlab="-log10(p)", main="")
axis(side=1, at=c(0,1, 3, 5, 7))
axis(side=2)
polygon(den, col="gray")
abline(v=3, lty=2)
abline(v=5, lty=2)
# abline(v=7, lty=2)

cex = 0.6
gapX = 0.6

# plot missense genes
genes = pMis[pMis>=3 & pMis < 5]
genes = genes[order(names(genes))]

printNamesOnPlot(names(genes), startX=3.4, startY=1.7, gapX=gapX, gapY=0.1, col="#e41a1c", genesPerRow = 3, cex=cex)

genes = pMis[pMis>=5 ]
genes = genes[order(names(genes))]
printNamesOnPlot(names(genes), startX=5.4, startY=1.7, gapX=gapX, gapY=0.1, col="#e41a1c", genesPerRow = 3, cex=cex)


# plot amp genes
genes = pAmp[pAmp>=3 & pAmp < 5]
genes = genes[order(names(genes))]

printNamesOnPlot(names(genes), startX=3.4, startY=0.9, gapX=gapX, gapY=0.1, col="#377eb8", genesPerRow = 3, cex=cex)

genes = pAmp[pAmp>=5 ]
genes = genes[order(names(genes))]
printNamesOnPlot(names(genes), startX=5.4, startY=0.9, gapX=gapX, gapY=0.1, col="#377eb8", genesPerRow = 3, cex=cex)

# plot trun genes
genes = pTru[pTru>=3 & pTru < 5]
genes = genes[order(names(genes))]

printNamesOnPlot(names(genes), startX=3.4, startY=0.3, gapX=gapX, gapY=0.1, col="#4daf4a", genesPerRow = 3, cex=cex)

genes = pTru[pTru>=5 ]
genes = genes[order(names(genes))]
printNamesOnPlot(names(genes), startX=5.4, startY=0.3, gapX=gapX, gapY=0.1, col="#4daf4a", genesPerRow = 3, cex=0.7)

dev.off()

