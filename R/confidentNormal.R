library(yaGST)

#' Get at most top 30 confident normal cells 
#'
#' @param NES 
#' @param pValue 
#' @param FDR 
#' @param pval_filter 
#' @param fdr_filter 
#' @param pval_cutoff 
#' @param nes_cutoff 
#' @param nNES 
#'
#' @return
#' @export
#'
#' @examples
top30classification <- function(NES, pValue, FDR, pval_filter, fdr_filter, pval_cutoff, nes_cutoff, nNES){
  
  if(fdr_filter == TRUE){
    FDR[FDR >= pval_cutoff] <- 1
  }
  if(nes_cutoff != 0){
    NES[NES < nes_cutoff] <- 0
  }
  
  comb <- NES * -log10(FDR)
  score_cells <- apply(comb, 2, max)
  conf_norm_cells <- score_cells[score_cells>0]
  
  Class <- c()
  if(length(conf_norm_cells)>0){
    conf_norm_cells <- sort(conf_norm_cells, decreasing = TRUE)[1:min(30, length(conf_norm_cells))]
    if(length(conf_norm_cells)==1){
      Class <- names(which.max(comb[,names(conf_norm_cells)]))
      names(Class) <- names(conf_norm_cells)
    }else{
      Class <- apply(comb[,names(conf_norm_cells)], 2, function(x) names(which.max(x)))
    }
    
  }
  return(Class)
}



#' getConfidentNormalCells Get at most top 30 confident normal cells from count matrix.
#'
#' @param mtx count matrix
#' @param sample sample name (optional)
#' @param par_cores = 20,
#' @param par_cores number of cores (default 20)
#' @param AdditionalGeneSets list of additional signatures of normal cell types (optional)
#' @param SCEVANsignatures FALSE if you only want to use only the signatures specified in AdditionalGeneSets (default TRUE)
#'
#' @return
#'
#' @examples
#' 
#' @export
getConfidentNormalCells <- function(mtx, sample = "", par_cores = 20, AdditionalGeneSets = NULL, SCEVANsignatures = TRUE, organism = "human", output_dir = "./output"){
  
  
  if(organism == "human"){
    geneSet <- geneSet 
  }else{
    geneSet <- geneSetMouse 
  }
  
  if(length(AdditionalGeneSets)>0 & SCEVANsignatures){
    geneSet <- append(geneSet, AdditionalGeneSets)
  }else if(length(AdditionalGeneSets)>0){
    geneSet <- AdditionalGeneSets
  }
  
  system.time(ssMwwGst(geData = mtx, geneSet = geneSet , ncore = par_cores, minLenGeneSet = 5, filename = file.path(output_dir, paste0(sample, "_confidentNormal")), standardize = TRUE))
  load(file.path(output_dir, paste0(sample, "_confidentNormal_MWW.RData")))
  norm.cell.names <- top30classification(NES = NES, FDR = FDR, pValue = pValue, fdr_filter = TRUE, pval_filter = TRUE, pval_cutoff = 1*10^(-10), nes_cutoff = 1.0, nNES = 1)
  file.remove(file.path(output_dir, paste0(sample, "_confidentNormal_MWW.RData")))
  print(paste("found", length(norm.cell.names), "confident non malignant cells", sep=" "))
  
  return(norm.cell.names)
  
}
