#### function to load
ssMwwGst <- function(geData, geneSet, ncore,  minLenGeneSet = 15 , regulon = FALSE, filename, standardize = TRUE){
  library(yaGST)
  if (standardize){
    means <- rowMeans(geData)
    sds <- apply(geData, 1, sd)
  }
  library(doMC)
  registerDoMC(ncore)
  ans <- foreach(ss = 1:ncol(geData)) %dopar% {
    # for(ss in 1:nSamples){
    if (standardize){
      currentSample <- (geData[, ss] - means)/sds
    }
    else{
      currentSample <- geData[, ss]
    }
    rankedList <- sort(currentSample, decreasing = T)
    if(regulon == FALSE){
      aMwwGST <- lapply(geneSet, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = minLenGeneSet, alternative = "two.sided"))
    }else{
      aMwwGST <- lapply(geneSet, function(x) mwwExtGST(rankedList, geneSetUp = x$pos, geneSetDown = x$neg, minLenGeneSet = minLenGeneSet))
    }
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
    ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
    return(ans)
  }
  NES <- sapply(ans, function(x) x$tmp_NES)
  pValue <- sapply(ans, function(x) x$tmp_pValue)
  colnames(NES) <- colnames(pValue) <- colnames(geData)
  FDR <- apply(pValue, 2, function(x) p.adjust(x, method = "fdr"))
  save(NES, pValue, FDR, file = paste0(filename, "_MWW.RData"))
}


NESclassification <- function(NES, pValue, FDR, pval_filter, fdr_filter, pval_cutoff, nes_cutoff, nNES){
  if(pval_filter == TRUE){
    NES[pValue >= pval_cutoff] <- 0
  }
  if(fdr_filter == TRUE){
    NES[FDR >= pval_cutoff] <- 0
  }
  if(nes_cutoff != 0){
    NES[NES < nes_cutoff] <- 0
  }
  if(nNES == 1){
    Class <- apply(NES, 2, function(x) names(which.max(x)))
  }
  
  print(sort(table(Class), decreasing = T))
  return(Class)
}