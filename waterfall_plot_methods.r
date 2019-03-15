getWaterfallSettings <- function(type = "mut", cols=NULL) {
  if(type == "mut") {
    if(is.null(cols)) {
      cols =  c("darkgray","#2c7bb6", "#fdae61", "#d7191c")
    }
    res =   list(colors = cols, legendText = c("WT", "Others", "Missense", "Truncating"))
  } else if(type == "cna") {
    if(is.null(cols)) {
      cols = c("darkgray","darkblue", "lightblue",  "orange", "darkred")
    }
    res = list(colors = cols, 
               legendText = c("no change", "2 copy deletion", "1 copy deletion", "amplification", "high-amplification"))
  }
  res
}

plotRandomRankBands <- function(n, handle, sig_alpha) {
  if(is.na(sig_alpha)) {
    sig_alpha = min(0.05, 1/n)
  }
  
  alpha = 1:n
  beta = n:1
  l = qbeta(sig_alpha/2, alpha, beta)
  h = qbeta(1-sig_alpha/2, alpha, beta)
  
  polygon(c(handle, rev(handle)), c( rev(l)-0.5,  rev(rev(h)-0.5) ), lty = 2 , lwd=2, col =rgb(1, 0, 0,0.35))
}


waterfallForGene <- function(panCancerData, gene, title, rank, legenedPos="bottomleft", 
                                 cols=NULL, type="all", sig_alpha = NA, cex.axis=1.1) {
  mut_tmp = rep(1, ncol(panCancerData$mutations_all))
  mut_tmp[ which(panCancerData$silentMutations[gene, ] == 1) ] = 2
  mut_tmp[ which(panCancerData$missenseMutations[gene, ] == 1) ] = 3
  mut_tmp[ which(panCancerData$truncatingMutations[gene, ] == 1) ] = 4

    
  if(rank == TRUE){
    allRanks = apply(panCancerData$viabilities, 2, rank) / nrow(panCancerData$viabilities) - 0.5
    via_tmp = as.numeric(allRanks[gene,])
  } else {
    via_tmp = as.numeric(panCancerData$viabilities[gene,])
  }
  
  via = via_tmp[order(via_tmp, decreasing=TRUE)]
  mut = mut_tmp[order(via_tmp, decreasing=TRUE)]
  
  if(type == "only_wt") {
    indexes = which(mut==1)  # 1 means wt
    mut = mut[indexes]
    via = via[indexes]
    if(length(indexes) == 0) {
      warnings("no observation!")
      return(NULL)
    }
    
  } else if(type == "only_mut") {
    indexes = which(mut!=1) # !=1 (i.e. 2, 3, 4) means mutated cases
    if(length(indexes) == 0) {
      warning("no observation!")
      return(NULL)
    }
    mut = mut[indexes]
    via = via[indexes]
  } else if(type == "only_truncating") {
    indexes = which(mut==4) 
    if(length(indexes) == 0) {
      warning("no observation!")
      return(NULL)
    } 
    mut = mut[indexes]
    via = via[indexes]
  } else if(type == "only_missense") {
    indexes = which(mut==3) 
    if(length(indexes) == 0) {
      warning("no observation!")
      return(NULL)
    } 
    mut = mut[indexes]
    via = via[indexes]
  }
  
  wfSettings = getWaterfallSettings(cols=cols)
  col = wfSettings$colors[mut]
  
  legText = wfSettings$legendText[as.numeric(names(table(mut)> 0))]
  legColors = wfSettings$colors[as.numeric(names(table(mut)> 0))]
  
  if(rank == TRUE){
    ylim = c(-0.6, 0.6)
    ylab = "Ranked Gene Score"
    yaxt = "n"
  } else {
    ylim=c(-max(abs(via)), max(abs(via)))
    ylab = "Gene Score"
    yaxt = NULL
  }
  handle = barplot(via, col=col, border=col, space=0.5, ylim=ylim,
                   main =title, ylab=ylab, yaxt=yaxt,
                   cex.axis=1.2, cex.lab=1.4, legend.text=legText,
                   args.legend=list(title="", fill=legColors, border=NA, cex=0.6, x = legenedPos,  bty = "n",xpd=FALSE ))

  n =  length(via)
  if(rank == TRUE){
    # add rank labels    
    axis_y = seq(-0.5, 0.5, length.out=5) 
    axis(2, at=axis_y, labels=c("0", "0.25", "0.5", "0.75", "1"), cex.axis=cex.axis)
    
    plotRandomRankBands(n, handle, sig_alpha)
  }
}




waterfallForGene_CNA <- function(panCancerData, gene, title, rank, legenedPos="bottomleft", cols=NULL, type="all",
                                 sig_alpha = NA, cex.axis=1.1) {
  cna_tmp = rep(1, ncol(panCancerData$copyNumbers))
  cna_tmp[ panCancerData$copyNumbers[gene,] == -2 ] = 2
  cna_tmp[ panCancerData$copyNumbers[gene,] == -1 ] = 3
  cna_tmp[ panCancerData$copyNumbers[gene,] == 1 ] = 4
  cna_tmp[ panCancerData$copyNumbers[gene,] == 2 ] = 5
  
  
  if(rank == TRUE){
    allRanks = apply(panCancerData$viabilities, 2, rank) / nrow(panCancerData$viabilities) - 0.5
    via_tmp = as.numeric(allRanks[gene,])
  } else {
    via_tmp = as.numeric(panCancerData$viabilities[gene,])
  }
  
  
  via = via_tmp[order(via_tmp, decreasing=TRUE)]
  cna = cna_tmp[order(via_tmp, decreasing=TRUE)]
  
  wfSettings = getWaterfallSettings("cna",cols=cols)
  col = wfSettings$colors[cna]
  
  
  legText = wfSettings$legendText[as.numeric(names(table(cna)> 0))]
  legColors = wfSettings$colors[as.numeric(names(table(cna)> 0))]
  
  
  if(rank == TRUE){
    ylim = c(-0.6, 0.6)
    ylab = "Ranked Gene Score"
    yaxt = "n"
  } else {
    ylim=c(-max(abs(via)), max(abs(via)))
    ylab = "Gene Score"
    yaxt = NULL
  }

  handle = barplot(via, col=col, border=col, space=0.5, ylim=ylim,
                   main =title, ylab=ylab, yaxt=yaxt,
                   cex.axis=1.2, cex.lab=1.4, legend.text=legText,
                   args.legend=list(title="", fill=legColors, border=NA, cex=0.6, x = legenedPos,  bty = "n",xpd=FALSE ))
  
  n =  length(via)
  
  # compute bands according to the beta distribution  
  if(rank == TRUE){
    axis_y = seq(-0.5, 0.5, length.out=5) 
    axis(2, at=axis_y, labels=c("0", "0.25", "0.5", "0.75", "1"), cex.axis=cex.axis)
    
    plotRandomRankBands(n, handle, sig_alpha)
  }
}
