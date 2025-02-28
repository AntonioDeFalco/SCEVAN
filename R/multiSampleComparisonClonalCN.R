
#' plotAllSubclonalCN Plot the copy number of each subclone of a sample.
#'
#' @param sample Name of the sample.
#' @param pathOutput Path to the output folder containing the output of pipelineCNA.
#'
#' @return
#' @export
#'
#' @examples 
#' 
plotAllSubclonalCN <- function(sample, pathOutput = "./output/"){
  
  allFile <- list.files(pathOutput, pattern = paste0(sample,"_subclone[1-9]_CN.seg") )

  CNVtot <- lapply(allFile, function(i) read.table(paste0(pathOutput,i), sep="\t", header=TRUE, as.is=TRUE))
  
  png(paste0(pathOutput,sample,"_compareSubclonalCN.png",sep=""), height=2250, width=1350, res=100)
  
  par(mfrow=c(length(allFile),1),cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  
  for(i in 1:length(allFile)){
    
    CNV <- CNVtot[[i]]
    
    plotSegmentation(CNV) 
    title(gsub("_CN.seg","",allFile[i]))
    
  }
  dev.off()
}


#' Title plotAllClonalCN
#'
#' @param samples Vector with sample names to be plotted
#' @param name Analysis name
#'
#' @return
#' @export
#'
#' @examples
plotAllClonalCN <- function(samples, name){
  
  #CNVtot <- lapply(samples, function(i) read.table(paste0("./output/"," ",i," _  _CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  CNVtot <- lapply(samples, function(i) read.table(paste0("./output/",i,"_Clonal_CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  
  png(paste0("./output/",name,"_compareClonalCN.png",sep=""), height=1550, width=2350, res=100)
  
  if(length(samples)>3){
    par(mfrow=c(3,ceiling(length(samples)/3)), cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  }else{
    par(mfrow=c(length(samples),1),cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  }
  
  for(i in 1:length(samples)){
    
    CNV <- CNVtot[[i]]
    
    plotSegmentation(CNV) 
    title(samples[i])
    
  }
  dev.off()
}

#' multiSampleComparisonClonalCN Compare the clonal Copy Number of multiple samples.
#'
#' @param listCountMtx Named list of raw count matrix of samples
#' @param analysisName Name of the analysis (default "all")
#' @param organism Organism to be analysed (optional - "mouse" or "human" - default "human")
#' @param par_cores number of cores (default 20)
#'
#' @return
#' @export
#'
#' @examples 
#' 
multiSampleComparisonClonalCN <- function(listCountMtx, listNormCells = NULL, analysisName = "all", organism = "human" , par_cores = 20, plotTree = TRUE, output_dir = "./output"){
  #TODO add the output_dir var here as well, apply it to the plotting func
  # plotCloneTreeNew
  
  resList <- lapply(names(listCountMtx), function(x) {
    pipelineCNA(listCountMtx[[x]], norm_cell = listNormCells[[x]], sample = x, SUBCLONES = FALSE, ClonalCN = TRUE, par_cores = par_cores, organism=organism)
  })
  names(resList) <- names(listCountMtx)
  
  sampleAlterList <- lapply(names(listCountMtx), function(x) {
    analyzeSegm(x, nSub = 0)
  })
  names(sampleAlterList) <- paste0(names(listCountMtx),"_subclone", 1:length(names(listCountMtx)))
  
  names(sampleAlterList) <- paste0(analysisName,"_subclone", 1:length(names(listCountMtx)))
  
  diffList <- diffSubclones(sampleAlterList, analysisName, nSub = length(names(listCountMtx)))
  
  diffList <- testSpecificAlteration(diffList, length(names(listCountMtx)), analysisName)
  
  genesMtx <- lapply(listCountMtx, function(x) rownames(x))
  
  genesMtx <- sort(unique(unlist(genesMtx)))
  
  genesMtx <- data.frame(x = rep(0, length(genesMtx)), row.names = genesMtx)
  
  annot_mtx <- annotateGenes(genesMtx)
  
  oncoHeat <- annoteBandOncoHeat(annot_mtx, diffList, length(names(listCountMtx)), organism)
  
  rownames(oncoHeat) <- names(listCountMtx)
  
  plotOncoHeatSubclones(oncoHeat, length(names(listCountMtx)), analysisName, NULL, organism)
  
  plotAllClonalCN(names(listCountMtx), name = analysisName)
  
  if(length(names(listCountMtx))>2 & plotTree) plotCloneTreeNew(names(listCountMtx), CLONAL_MULTI = TRUE, analysisName = analysisName)
  
  for(i in 1:length(names(listCountMtx))){
    names(diffList) <- gsub(paste0("subclone",i), names(listCountMtx)[i], names(diffList))
  }
  
  names(diffList) <- gsub("clone", "shared", names(diffList))

  outputAnalysis <- list(resList, diffList)
  
  
  save(outputAnalysis, file = file.path(output_dir, paste0(analysisName, "_outputAnalysis.RData")))
  
  outputAnalysis
  
}


