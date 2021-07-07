library(yaGST)

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
    }else{
      Class <- apply(comb[,names(conf_norm_cells)], 2, function(x) names(which.max(x)))
    }
    
  }
  return(Class)
}



#' func
#'
#' @param 
#' @param 
#' @param 
#'
#' @return
#'
#' @examples
#' 
#' @export
getConfidentNormalCells <- function(mtx, par_cores = 20){
  
  system.time(ssMwwGst(geData = mtx, geneSet = geneSet , ncore = par_cores, minLenGeneSet = 5, filename = "./output/confidentNormal", standardize = FALSE))
  load("./output/confidentNormal_MWW.RData")
  norm.cell.names <- top30classification(NES = NES, FDR = FDR, pValue = pValue, fdr_filter = TRUE, pval_filter = TRUE, pval_cutoff = 0.0005, nes_cutoff = 0.5, nNES = 1)
  
  print(paste("found", length(norm.cell.names), "confident non malignant cells", sep=" "))
  
  return(norm.cell.names)
  
}
