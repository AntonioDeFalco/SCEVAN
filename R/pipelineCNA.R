#' SCEVAN: R package that automatically classifies the cells in the scRNA data by segregating non-malignant cells of tumor microenviroment from the malignant cells. It also infers the copy number profile of malignant cells, identifies subclonal structures and analyses the specific and shared alterations of each subpopulation.
#'
#' The SCEVAN package
#' 
#' @section SCEVAN functions:
#' The SCEVAN functions ...
#'
#' @docType package
#' @name SCEVAN
#' @useDynLib SCEVAN, .registration=TRUE
NULL
#> NULL


#' pipelineCNA Executes the entire SCEVAN pipeline that classifies tumour and normal cells from the raw count matrix, infer the clonal profile of cancer cells and looks for possible sub-clones in the tumour cell matrix automatically analysing the specific and shared alterations of each subclone and a differential analysis of pathways and genes expressed in each subclone.
#'
#' @param count_mtx Raw count matrix with genes on rows (both Gene Symbol or Ensembl ID are allowed) and cells on columns.
#' @param sample Sample name to save results (optional)
#' @param par_cores Number of cores to run the pipeline (optional - default 20)
#' @param norm_cell Vector of possible known normal cells to be used as confident normal cells (optional)
#' @param SUBCLONES Boolean value TRUE if you are interested in analysing the clonal structure and FALSE if you are only interested in the classification of malignant and non-malignant cells (optional - default TRUE)
#' @param beta_vega Specifies beta parameter for segmentation, higher beta for more coarse-grained segmentation. (optional - default 0.5)
#' @param ClonalCN Get clonal CN profile inference from all tumour cells (optional)
#' @param plotTree Plot Phylogenetic tree (optional - default FALSE)
#' @param AdditionalGeneSets list of additional signatures of normal cell types (optional)
#' @param SCEVANsignatures FALSE if you only want to use only the signatures specified in AdditionalGeneSets (default TRUE)
#' @param organism Organism to be analysed (optional - "mouse" or "human" - default "human")
#' @param ngenes_chr Minimum number of genes expressed on chromosome (optional - default 5)
#' @param perc_genes Minimum percentage gene expressed in each cell (optional - default 10)
#' @param FIXED_NORMAL_CELLS TRUE if norm_cell vector to be used as fixed reference, if you are only interested in clonal structure and not normal/tumor classification (default FALSE)
#'
#' @return
#' @export
#'
#' @examples res_pip <- pipelineCNA(count_mtx)

pipelineCNA <- function(count_mtx, sample="", par_cores = 20, norm_cell = NULL, SUBCLONES = TRUE, beta_vega = 0.5, ClonalCN = TRUE, plotTree = TRUE, AdditionalGeneSets = NULL, SCEVANsignatures = TRUE, organism = "human", ngenes_chr = 5, perc_genes = 10, FIXED_NORMAL_CELLS = FALSE, output_dir = "./output"){
  
  dir.create(file.path(output_dir), showWarnings = FALSE, recursive = TRUE)
  
  start_time <- Sys.time()
  
  normalNotKnown <- length(norm_cell)==0
  
  res_proc <- preprocessingMtx(count_mtx,sample, par_cores=par_cores, findConfident = normalNotKnown, AdditionalGeneSets = AdditionalGeneSets, SCEVANsignatures = SCEVANsignatures, organism = organism, ngenes_chr = ngenes_chr, perc_genes = (perc_genes/100), output_dir = output_dir)
  
  if(normalNotKnown) norm_cell <- names(res_proc$norm_cell)
  
  res_class <- classifyTumorCells(res_proc$count_mtx_norm, res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = NULL,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE, beta_vega = beta_vega, FIXED_NORMAL_CELLS = FIXED_NORMAL_CELLS, output_dir = output_dir)
  
  print(paste("found", length(res_class$tum_cells), "tumor cells"))
  classDf <- data.frame(class = rep("filtered", length(colnames(res_proc$count_mtx))), row.names = colnames(res_proc$count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  classDf[res_class$confidentNormal, "confidentNormal"] <- "yes"
  
  end_time<- Sys.time()
  print(paste("time classify tumor cells: ", end_time -start_time))
  
  if(ClonalCN) getClonalCNProfile(res_class, res_proc, sample, par_cores, organism = organism, output_dir = output_dir)
  
  mtx_vega <- segmTumorMatrix(res_proc, res_class, sample, par_cores, beta_vega)
  
  if (SUBCLONES) {
    res_subclones <- subcloneAnalysisPipeline(res_proc$count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf, beta_vega, plotTree, organism, output_dir)
    FOUND_SUBCLONES <- res_subclones$FOUND_SUBCLONES
    classDf <- res_subclones$classDf
  }else{
    FOUND_SUBCLONES <- FALSE
  }
  
  if(!FOUND_SUBCLONES) plotCNclonal(sample,ClonalCN, organism)
  
  #save CNA matrix
  #CNAmtx <- res_class$CNAmat[,-c(1,2,3)]
  #save(CNAmtx, file = paste0("./output/",sample,"_CNAmtx.RData"))
  
  #save annotated matrix
  count_mtx_annot <- res_proc$count_mtx_annot
  
  save(count_mtx_annot, file = file.path(output_dir, paste0(sample, "_count_mtx_annot.RData")))
  
  
  #remove intermediate files
  mtx_vega_files <- list.files(path = output_dir, pattern = "_mtx_vega")
  sapply(mtx_vega_files, function(x) file.remove(file.path(output_dir, x)))
  
  return(classDf)
}


getClonalCNProfile <- function(res_class, res_proc, sample, par_cores, beta_vega = 3, organism = "human", output_dir = "./output"){
  
  mtx <- res_class$CNAmat[,res_class$tum_cells]
  # hcc <- hclust(parallelDist::parDist(t(mtx),threads =par_cores, method = "euclidean"), method = "ward.D")
  # hcc2 <- cutree(hcc,2)
  # clonalClust <- as.integer(names(which.max(table(hcc2))))
  # mtx <- mtx[,names(hcc2[hcc2==clonalClust])]
  
  mtx_vega <- cbind(res_class$CNAmat[,1:3], mtx)
  colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
  breaks_tumor <- getBreaksVegaMC(mtx_vega, res_class$CNAmat[,3], paste0(sample,"ClonalCNProfile"), beta_vega = beta_vega, output_dir = output_dir)
  
  #mtx_CNA3 <- computeCNAmtx(mtx, breaks_tumor, par_cores, rep(TRUE, length(breaks_tumor)))
  
  #colnames(mtx_CNA3) <- res_class$tum_cells
  
  #save(mtx_CNA3, file = paste0("./output/",sample,"_mtx_CNA3.RData"))
  
  CNV <- getCNcall(mtx, res_proc$count_mtx_annot, breaks_tumor, sample = sample, CLONAL = TRUE, organism = organism)
  
  segm.mean <- getScevanCNV(sample, beta = "ClonalCNProfile")$Mean
  CNV <- cbind(CNV,segm.mean)
  
  write.table(CNV, file = file.path(output_dir, paste0(sample, "_Clonal_CN.seg")) , sep = "\t", quote = FALSE)
  
  #TODO I need to trace back where this file is created before this is changed,
  # as the removal might be important for size considerations. 

  # file.path(output_dir, paste0(sample, "ClonalCNProfile"))
  file.remove(paste0("./output/ ",paste0(sample,"ClonalCNProfile")," vega_output"))
  
  CNV
}

segmTumorMatrix <- function(res_proc, res_class, sample, par_cores, beta_vega = 0.5){
  
  mtx_vega <- cbind(res_class$CNAmat[,1:3], res_class$CNAmat[,res_class$tum_cells])
  colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
  breaks_tumor <- getBreaksVegaMC(mtx_vega, res_proc$count_mtx_annot[,3], paste0(sample,"onlytumor"), beta_vega = beta_vega, output_dir = output_dir)
  
  subSegm <- read.csv(paste0("./output/ ",paste0(sample,"onlytumor")," vega_output"), sep = "\t")
  #segmAlt <- abs(subSegm$Mean)>=0.10 | (subSegm$G.pv<=0.5 | subSegm$L.pv<=0.5)
  segmAlt <- (subSegm$G.pv<=0.5 | subSegm$L.pv<=0.5)
  mtx_vega <- computeCNAmtx(res_class$CNAmat[,res_class$tum_cells], breaks_tumor, par_cores,segmAlt ) #rep(TRUE, length(breaks_tumor))
  
  #mtx_vega <- computeCNAmtx(res_class$CNAmat[,res_class$tum_cells], breaks_tumor, par_cores, rep(TRUE, length(breaks_tumor)))
  colnames(mtx_vega) <- colnames(res_class$CNAmat[,res_class$tum_cells])
  rownames(mtx_vega) <- rownames(res_class$CNAmat[,res_class$tum_cells])
  hcc <- hclust(parallelDist::parDist(t(mtx_vega),threads = par_cores, method = "euclidean"), method = "ward.D")
  plotCNA(res_proc$count_mtx_annot$seqnames, mtx_vega, hcc, paste0(sample,"onlytumor"))
  
  return(mtx_vega)
}


subcloneAnalysisPipeline <- function(count_mtx, res_class, res_proc, mtx_vega,  sample, par_cores, classDf, beta_vega, plotTree = FALSE, organism  = "human", output_dir = "./output"){
  
  #save(count_mtx, res_class, res_proc, mtx_vega,  sample, par_cores, classDf, beta_vega, file = "debug_subclonesTumorCells.RData")
  
  start_time <- Sys.time()
  
  FOUND_SUBCLONES <- FALSE
  
  res_subclones <- subclonesTumorCells(res_class$tum_cells, res_class$CNAmat, sample, par_cores, beta_vega, res_proc, NULL, mtx_vega, organism = organism)
  
  if(length(setdiff(res_class$tum_cells,names(res_subclones$clustersSub)))>0){
    classDf[setdiff(res_class$tum_cells,names(res_subclones$clustersSub)),]$class <- "normal"
    res_class$tum_cells <- names(res_subclones$clustersSub)
  }
  
  tum_cells <- res_class$tum_cells
  clustersSub <- res_subclones$clustersSub
  #save(tum_cells,clustersSub, file = paste0(sample,"subcl.RData"))
  
  if(res_subclones$n_subclones>1){
    sampleAlter <- analyzeSegm(sample, nSub = res_subclones$n_subclones)
    
    if(length(sampleAlter)>1){
      
      diffSubcl <- diffSubclones(sampleAlter, sample, nSub = res_subclones$n_subclones)
      
      diffSubcl <- testSpecificAlteration(diffSubcl, res_subclones$n_subclones, sample)
      
      print(diffSubcl)
      
      ## new aggregation subclone
      oncoHeat <- annoteBandOncoHeat(res_proc$count_mtx_annot, diffSubcl, res_subclones$n_subclones, organism)
      
      res <- list()
      for (sub in 1:nrow(oncoHeat)){
        res[[sub]] <- apply(oncoHeat[-sub,], 1, function(x) sum(oncoHeat[sub,]==x) == ncol(oncoHeat))
      }
      if(any(unlist(lapply(res, function(x) any(x))))){
        
        shInd <- unlist(lapply(res, function(x) any(x)))
        removInd <- c()
        for(ind in which(shInd)){
          shNam <- names(res[[ind]][res[[ind]]>0]) 
          indSh <- as.numeric(substr(shNam,nchar(shNam[1]),nchar(shNam[1])))
          
          for(ind2 in indSh){          
            if(ind2>ind){
              res_subclones$clustersSub[res_subclones$clustersSub==ind2] <- ind
              removInd <- append(removInd,ind2)
            }
          }
          
        }
        unique(res_subclones$clustersSub)
        for(i in 1:length(removInd)){
          res_subclones$clustersSub[res_subclones$clustersSub>(removInd[i]-(i-1))] <- res_subclones$clustersSub[res_subclones$clustersSub>(removInd[i]-(i-1))] - 1
        }
        res_subclones$n_subclones <- length(unique(res_subclones$clustersSub))
        
        #remove previous segm file
        
        mtx_vega_files <- list.files(path = output_dir, pattern = "vega")
        mtx_vega_files <- mtx_vega_files[grep(sample,mtx_vega_files)]
        mtx_vega_files <- mtx_vega_files[grep("subclone", mtx_vega_files)]
        sapply(mtx_vega_files, function(x) file.remove(file.path(output_dir, x)))
        
        if(res_subclones$n_subclones>1){
          #res_subclones <- ReSegmSubclones(res_class$tum_cells, res_class$CNAmat, sample, res_subclones$clustersSub, par_cores, beta_vega)
          
          res_subclones <- subclonesTumorCells(res_class$tum_cells, res_class$CNAmat, sample, par_cores, beta_vega, res_proc, res_subclones$clustersSub, organism = organism)
          
          sampleAlter <- analyzeSegm(sample, nSub = res_subclones$n_subclones)
          
          if(length(sampleAlter)>1){
            diffSubcl <- diffSubclones(sampleAlter, sample, nSub = res_subclones$n_subclones)
            diffSubcl <- testSpecificAlteration(diffSubcl, res_subclones$n_subclones, sample)
          } 
        }
      }
      
      if(res_subclones$n_subclones>1){
        #segmList <- lapply(1:res_subclones$n_subclones, function(x) read.table(paste0("./output/ ",sample,"_subclone",x," vega_output"), sep="\t", header=TRUE, as.is=TRUE))
        #names(segmList) <- paste0("subclone",1:res_subclones$n_subclones)
        
        #save(res_proc, res_subclones, segmList,diffSubcl,sample, file = "plotcnaline.RData")
        #plotCNAline(segmList, diffSubcl, sample, res_subclones$n_subclones)
        
        #diffSubcl[[grep("_clone",names(diffSubcl))]] <- diffSubcl[[grep("_clone",names(diffSubcl))]][1:min(10,nrow(diffSubcl[[grep("_clone",names(diffSubcl))]])),]
        
        plotConsensusCNA(samp = sample, nSub = res_subclones$n_subclones, pathOutput = output_dir, organism = organism)
        
        perc_cells_subclones <- table(res_subclones$clustersSub)/length(res_subclones$clustersSub)
        
        oncoHeat <- annoteBandOncoHeat(res_proc$count_mtx_annot, diffSubcl, res_subclones$n_subclones, organism)
        #save(oncoHeat, file = paste0(sample,"_oncoheat.RData"))
        plotOncoHeatSubclones(oncoHeat, res_subclones$n_subclones, sample, perc_cells_subclones, organism)
        
        #TODO change to output_dir var
        plotTSNE(count_mtx, res_class$CNAmat, rownames(res_proc$count_mtx_norm), res_class$tum_cells, res_subclones$clustersSub, sample)
        classDf[names(res_subclones$clustersSub), "subclone"] <- res_subclones$clustersSub
        if(res_subclones$n_subclones>2 & plotTree) plotCloneTree(sample, res_subclones)
        
        if (length(grep("subclone",names(diffSubcl)))>0) genesDE(res_proc$count_mtx_norm, res_proc$count_mtx_annot, res_subclones$clustersSub, sample, diffSubcl[grep("subclone",names(diffSubcl))])
        pathwayAnalysis(res_proc$count_mtx_norm, res_proc$count_mtx_annot, res_subclones$clustersSub, sample, organism = organism)
        
        save(diffSubcl, file = file.path(output_dir, paste0(sample, "_SubcloneDiffAnalysis.RData")))
        
        FOUND_SUBCLONES <- TRUE
      }else{
        print("no significant subclones")
      }
      
    }else{
      print("no significant subclones")
    }
    
  }
  
  end_time<- Sys.time()
  print(paste("time subclones: ", end_time -start_time))
  
  res <- list(FOUND_SUBCLONES, classDf)
  names(res) <- c("FOUND_SUBCLONES","classDf")
  
  return(res)
}
