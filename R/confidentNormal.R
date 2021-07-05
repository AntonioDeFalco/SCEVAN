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


createGeneSetNormal <- function(){
  
  #ESTIMATE
  load("/storage/qnap_home/adefalco/singleCell/AllData/ClassTumorCells/SI_geneset.RData") #from https://sourceforge.net/projects/estimateproject/
  geneSet <- c()
  geneSet$stromal <- as.character(unlist(SI_geneset["StromalSignature",-1]))
  geneSet$immune <- as.character(unlist(SI_geneset["ImmuneSignature",-1]))
  
  #FRANCESCA & CANCER CELLS
  load("/storage/qnap_home/caruso/Analisi2021/signatures/scTHI_c8_signatures_968.RData")
  load("/storage/qnap_home/caruso/Analisi2021/CPTAC/FgesSignature/fges_signature.RData")
  
  findSign <- function(Phenotype){
    all_sign <- rownames(signature_Colors[signature_Colors$ALLPhenotypeFinal==Phenotype,])
    ind <- which(names(signature) %in% all_sign)
    subset_sign <- signature[ind]
    num_sin <- length(subset_sign)/2
    subset_sign <- unlist(subset_sign)
    
    i <- 1
    while( (length(unique(subset_sign)) > 400) & (i < num_sin)){
      subset_sign <- subset_sign[duplicated(subset_sign)]
      i <- i + 1
    }
    
    return(unique(subset_sign))
  }
  
  all_sign <- rownames(signature_Colors[signature_Colors$ALLPhenotypeFinal=="Oligodendrocytes",])
  ind <- which(names(signature) %in% all_sign)
  subset_sign <- signature[ind]
  
  geneSet$olig_Myelinating <- union(subset_sign$Anna_PreMyelinatingOligo, subset_sign$Anna_MyelinatingOligo)
  geneSet$olig_Myelinating <- union(geneSet$olig_Myelinating, subset_sign$CNS_Myelinating.Oligodendrocytes) 
  
  geneSet$olig <- union(subset_sign$Anna_OligoLineage, subset_sign$CNS_Oligodendrocytes)
  geneSet$olig <- union(geneSet$olig, subset_sign$CNS_Newly.Formed.Oligodendrocyte)
  
  #geneSet$olig <- findSign("Oligodendrocytes")
  geneSet$Tcell <- findSign("Tcell")
  #geneSet$astro <- findSign("Astrocytes")
  geneSet$macro <-findSign("Macrophages")
  geneSet$micro <-findSign("Microglia")
  geneSet$neuro <-findSign("Neurons")
  
  return(geneSet)
  
}



getConfidentNormalCells <- function(mtx, geneSet, n.cores = 20){
  
  system.time(ssMwwGst(geData = mtx, geneSet = geneSet , ncore = n.cores, minLenGeneSet = 5, filename = "confidentNormal", standardize = FALSE))
  load("confidentNormal_MWW.RData")
  norm.cell.names <- top30classification(NES = NES, FDR = FDR, pValue = pValue, fdr_filter = TRUE, pval_filter = TRUE, pval_cutoff = 0.0005, nes_cutoff = 0.5, nNES = 1)
  
  print(paste("found", length(norm.cell.names), "confident non malignant cells", sep=" "))
  
  return(norm.cell.names)
  
}
