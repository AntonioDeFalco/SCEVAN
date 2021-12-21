#' SCEVAN: A package for
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




#' Run pipeline runs the pipeline that classifies tumour and normal cells from the raw count matrix and looks for possible sub-clones in the tumour cell matrix
#'
#' @param count_mtx raw count matrix
#' @param sample sample name (optional)
#' @param par_cores number of cores (optional)
#' @param norm_cell vector normal cells if known (optional)
#' @param gr_truth ground truth of classification (optional)
#' @param SUBCLONES find subclones (optional)
#'
#' @return
#' @export
#'
#' @examples res_pip <- pipelineCNA(count_mtx, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)

pipelineCNA <- function(count_mtx, sample="", par_cores = 20, norm_cell = NULL, SUBCLONES = TRUE){
  
  dir.create(file.path("./output"), showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  res_proc <- preprocessingMtx(count_mtx, par_cores=par_cores)
  
  if(length(norm_cell)==0) norm_cell <- names(res_proc$norm_cell)

  res_class <- classifyTumorCells(res_proc$count_mtx_norm, res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = NULL,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)
  
  print(paste("found", length(res_class$tum_cells), "tumor cells"))
  classDf <- data.frame(class = rep("filtered", length(colnames(count_mtx))), row.names = colnames(count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  classDf[res_class$confidentNormal, "confidentNormal"] <- "yes"
  
  end_time<- Sys.time()
  print(paste("time classify tumor cells: ", end_time -start_time))

  mtx_vega <- segmTumorMatrix(res_proc, res_class, sample, par_cores)

  if (SUBCLONES) {
    res_subclones <- subcloneAnalysisPipeline(count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf)
    FOUND_SUBCLONES <- res_subclones$FOUND_SUBCLONES
    classDf <- res_subclones$classDf
  }else{
    FOUND_SUBCLONES <- FALSE
  }
  
  if(!FOUND_SUBCLONES) plotCNAlineOnlyTumor(sample)
  
  #save CNA matrix
  CNAmtx <- res_class$CNAmat[,-c(1,2,3)]
  save(CNAmtx, file = paste0("./output/",sample,"_CNAmtx.RData"))
  
  #remove intermediate files
  mtx_vega_files <- list.files(path = "./output/", pattern = "_mtx_vega")
  sapply(mtx_vega_files, function(x) file.remove(paste0("./output/",x)))
  
  return(classDf)
}


segmTumorMatrix <- function(res_proc, res_class, sample, par_cores){
  
  mtx_vega <- cbind(res_class$CNAmat[,1:3], res_class$CNAmat[,res_class$tum_cells])
  colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
  breaks_tumor <- getBreaksVegaMC(mtx_vega, res_proc$count_mtx_annot[,3], paste0(sample,"onlytumor"))
  
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


subcloneAnalysisPipeline <- function(count_mtx, res_class, res_proc, mtx_vega,  sample, par_cores, classDf){
  
  start_time <- Sys.time()
  
  FOUND_SUBCLONES <- FALSE
    
  res_subclones <- subclonesTumorCells(res_class$tum_cells, res_class$CNAmat,mtx_vega, sample, par_cores)
  
  tum_cells <- res_class$tum_cells
  clustersSub <- res_subclones$clustersSub
  #save(tum_cells,clustersSub, file = paste0(sample,"subcl.RData"))
  
  if(res_subclones$n_subclones>1){
    sampleAlter <- analyzeSegm(sample, nSub = res_subclones$n_subclones)
    
    if(length(sampleAlter)>1){
      
      diffSubcl <- diffSubclones(sampleAlter, sample, nSub = res_subclones$n_subclones)
      
      diffSubcl <- testSpecificAlteration(res_class$CNAmat, res_proc$count_mtx_annot, diffSubcl, res_subclones$clustersSub, res_subclones$n_subclones, sample)
      
      print(diffSubcl)
      
      ## new aggregation subclone
      oncoHeat <- annoteBandOncoHeat(res_proc$count_mtx_annot, diffSubcl, res_subclones$n_subclones)
      
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
        mtx_vega_files <- list.files(path = "./output/", pattern = "vega")
        mtx_vega_files <- mtx_vega_files[grep(sample,mtx_vega_files)]
        mtx_vega_files <- mtx_vega_files[grep("subclone",mtx_vega_files)]
        sapply(mtx_vega_files, function(x) file.remove(paste0("./output/",x)))
        
        if(res_subclones$n_subclones>1){
        res_subclones <- ReSegmSubclones(res_class$tum_cells, res_class$CNAmat, sample, res_subclones$clustersSub, par_cores)
        sampleAlter <- analyzeSegm(sample, nSub = res_subclones$n_subclones)
        
          if(length(sampleAlter)>1){
            diffSubcl <- diffSubclones(sampleAlter, sample, nSub = res_subclones$n_subclones)
            diffSubcl <- testSpecificAlteration(res_class$CNAmat, res_proc$count_mtx_annot, diffSubcl, res_subclones$clustersSub, res_subclones$n_subclones, sample)
          } 
        }
      }
      
      if(res_subclones$n_subclones>1){
      segmList <- lapply(1:res_subclones$n_subclones, function(x) read.table(paste0("./output/ ",sample,"_subclone",x," vega_output"), sep="\t", header=TRUE, as.is=TRUE))
      names(segmList) <- paste0("subclone",1:res_subclones$n_subclones)
      
      #save(res_proc, res_subclones, segmList,diffSubcl,sample, file = "plotcnaline.RData")
      plotCNAline(segmList, diffSubcl, sample, res_subclones$n_subclones)
      
      diffSubcl[[grep("_clone",names(diffSubcl))]] <- diffSubcl[[grep("_clone",names(diffSubcl))]][1:min(10,nrow(diffSubcl[[grep("_clone",names(diffSubcl))]])),]
      
      perc_cells_subclones <- table(res_subclones$clustersSub)/length(res_subclones$clustersSub)
      
      oncoHeat <- annoteBandOncoHeat(res_proc$count_mtx_annot, diffSubcl, res_subclones$n_subclones)
      #save(oncoHeat, file = paste0(sample,"_oncoheat.RData"))
      plotOncoHeatSubclones(oncoHeat, res_subclones$n_subclones, sample, perc_cells_subclones)
      
      plotTSNE(count_mtx, res_class$CNAmat, rownames(res_proc$count_mtx_norm), res_class$tum_cells, res_subclones$clustersSub, sample)
      classDf[names(res_subclones$clustersSub), "subclone"] <- res_subclones$clustersSub
      if(res_subclones$n_subclones>2) plotCloneTree(sample, res_subclones)
      
      if (length(grep("subclone",names(diffSubcl)))>0) genesDE(res_proc$count_mtx_norm, res_proc$count_mtx_annot, res_subclones$clustersSub, sample, diffSubcl[grep("subclone",names(diffSubcl))])
      pathwayAnalysis(res_proc$count_mtx_norm, res_proc$count_mtx_annot, res_subclones$clustersSub, sample)
      
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

#' compareClonalStructure
#'
#' @param count_mtx raw count matrix
#' @param sample sample name (optional)
#' @param par_cores number of cores (optional)
#' @param SUBCLONES find subclones (optional)
#'
#' @return
#' @export
#'
#' @examples 
#' 
compareClonalStructure <- function(count_mtx1, count_mtx2 , samp_1="", samp_2="", par_cores = 20){

  nSub <- 2
  
  count_mtx <- merge(count_mtx1, count_mtx2, by="row.names" ,all = TRUE)
  count_mtx[is.na(count_mtx)] <- 0
  rownames(count_mtx) <- count_mtx$Row.names
  count_mtx<- count_mtx[,-1]
  
  print(samp_1)
  res_proc_1 <- preprocessingMtx(count_mtx1, par_cores=par_cores)
  norm_cell <- names(res_proc_1$norm_cell)
  res_class_1 <- classifyTumorCells(res_proc_1$count_mtx_norm,res_proc_1$count_mtx_annot, samp_1, par_cores=par_cores,  ground_truth = NULL,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)
  print(paste("found", length(res_class_1$tum_cells), "tumor cells"))
  
  print(samp_2)
  res_proc_2 <- preprocessingMtx(count_mtx2, par_cores=par_cores)
  norm_cell <- names(res_proc_2$norm_cell)
  res_class_2 <- classifyTumorCells(res_proc_2$count_mtx_norm,res_proc_2$count_mtx_annot, samp_2, par_cores=par_cores,  ground_truth = NULL,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)
  print(paste("found", length(res_class_2$tum_cells), "tumor cells"))
  
  sampl <- paste0(samp_1,"-vs-", samp_2)
  
  mtx_vega <- segmTumorMatrix(res_proc_1, res_class_1, samp_1, par_cores)
  mtx_vega <- segmTumorMatrix(res_proc_2, res_class_2, samp_2, par_cores)
  
  all_segm <- c()
  segmList <- list()
  segmList$subclone1 <- read.csv(paste0("./output/ ",samp_1,"onlytumor vega_output"), sep = "\t")
  all_segm[[paste0(sampl,"_subclone", 1)]] <- getPossibleSpecAltFromSeg(segmList$subclone1)
  segmList$subclone2 <- read.csv(paste0("./output/ ",samp_2,"onlytumor vega_output"), sep = "\t")
  all_segm[[paste0(sampl,"_subclone", 2)]] <- getPossibleSpecAltFromSeg(segmList$subclone2)
  
  all_segm <- diffSubclones(all_segm, sampl) 
  
  res_class <- list()
  res_class$CNAmat <- merge(res_class_1$CNAmat, res_class_2$CNAmat, by = "gene_id", all = TRUE)
  res_class$CNAmat[is.na(res_class$CNAmat)] <- 0
  
  res_proc <- list()
  res_proc$count_mtx_annot <- rbind(res_proc_1$count_mtx_annot, res_proc_2$count_mtx_annot)
  res_proc$count_mtx_annot <- res_proc$count_mtx_annot[!duplicated(res_proc$count_mtx_annot$gene_id),]
  res_proc$count_mtx_annot <- res_proc$count_mtx_annot[order(res_proc$count_mtx_annot$seqnames,res_proc$count_mtx_annot$start),]
  
  res_proc$count_mtx_norm <- merge(res_proc_1$count_mtx_norm, res_proc_2$count_mtx_norm, by="row.names" , all = TRUE)
  res_proc$count_mtx_norm[is.na(res_proc$count_mtx_norm)] <- 0
  
  rownames(res_class$CNAmat) <- res_class$CNAmat$gene_id
  res_class$CNAmat <- res_class$CNAmat[res_proc$count_mtx_annot$gene_id,]
  
  rownames(res_proc$count_mtx_norm) <- res_proc$count_mtx_norm$Row.names
  res_proc$count_mtx_norm <- res_proc$count_mtx_norm[,-1]
  res_proc$count_mtx_norm <- res_proc$count_mtx_norm[res_proc$count_mtx_annot$gene_name,]
  
  clust_subclones <- append(rep(1,length(res_class_1$tum_cells)),rep(2,length(res_class_2$tum_cells)))
  names(clust_subclones) <- append(res_class_1$tum_cells,res_class_2$tum_cells)
  
  diffSubcl <- testSpecificAlteration(res_class$CNAmat, res_proc$count_mtx_annot, all_segm, clust_subclones, nSub, sampl)
  
  colors_samp <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Set2")[6:7])
  plotCNAline(segmList, diffSubcl, sampl, nSub, colors_samp)
  
  diffSubcl[[grep("_clone",names(diffSubcl))]] <- diffSubcl[[grep("_clone",names(diffSubcl))]][1:min(10,nrow(diffSubcl[[grep("_clone",names(diffSubcl))]])),]
  
  oncoHeat <- annoteBandOncoHeat(res_proc$count_mtx_annot, diffSubcl, nSub)
  
  annotdf <- data.frame(row.names = rownames(oncoHeat), 
                        Sample = c(samp_1,samp_2) )  
  
  
  subclones <- colors_samp(nSub)
  names(subclones) <- unique(annotdf$Sample)
  mycolors <- list(Sample = subclones)
  
  ppOncoHeat <- plotOncoHeat(oncoHeat, nSub, sampl, annotdf, mycolors)
  print(sampl)
  
  classDf1 <- data.frame(class = rep("filtered", length(colnames(count_mtx1))), row.names = colnames(count_mtx1))
  classDf1[colnames(res_class_1$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf1[res_class_1$tum_cells, "class"] <- "tumor"
  classDf1[res_class_1$confidentNormal, "confidentNormal"] <- "yes"
  
  classDf2 <- data.frame(class = rep("filtered", length(colnames(count_mtx2))), row.names = colnames(count_mtx2))
  classDf2[colnames(res_class_2$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf2[res_class_2$tum_cells, "class"] <- "tumor"
  classDf2[res_class_2$confidentNormal, "confidentNormal"] <- "yes"
  
  return(rbind(classDf1,classDf2))
  
}






#' Run pipeline runs the pipeline that classifies tumour and normal cells from the raw count matrix and looks for possible sub-clones in the tumour cell matrix
#'
#' @param count_mtx raw count matrix
#' @param sample sample name (optional)
#' @param par_cores number of cores (optional)
#' @param norm_cell vector normal cells if known (optional)
#' @param gr_truth ground truth of classification (optional)
#' @param SUBCLONES find subclones (optional)
#'
#' @return
#' @export
#'
#' @examples res_pip <- pipelineCNA(count_mtx, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)

DEBUGpipelineCNA <- function(count_mtx, sample="", par_cores = 20, norm_cell = NULL,  gr_truth = NULL, SUBCLONES = TRUE){   

  dir.create(file.path("./output"), showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  res_proc <- preprocessingMtx(count_mtx, par_cores=par_cores)
  
  if(length(norm_cell)==0) norm_cell <- names(res_proc$norm_cell)
  
  print(table(gr_truth[norm_cell]))
  
  res_class <- classifyTumorCells(res_proc$count_mtx_norm,res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = gr_truth,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)

  print(paste("found", length(res_class$tum_cells), "tumor cells"))
  classDf <- data.frame(class = rep("filtered", length(colnames(count_mtx))), row.names = colnames(count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  classDf[res_class$confidentNormal, "confidentNormal"] <- "yes"
  
  end_time<- Sys.time()
  print(paste("time classify tumor cells: ", end_time -start_time))

  # DEBUG
  if(length(gr_truth)>0){
    ground_truth_mal <- names(gr_truth[gr_truth == "malignant"])
    pred_mal <- res_class$tum_cells
    F1_Score <- computeF1score(pred_mal, ground_truth_mal)
    print(paste("F1_Score: ", F1_Score))
  }
  
  mtx_vega <- segmTumorMatrix(res_proc, res_class, sample, par_cores)
  
  #save(res_class, res_proc, mtx_vega, sample, par_cores, classDf, file = paste0("beforeSubclone", sample,".RData"))
  
  if (SUBCLONES) {
    res_subclones <- subcloneAnalysisPipeline(count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf)
    FOUND_SUBCLONES <- res_subclones$FOUND_SUBCLONES
    classDf <- res_subclones$classDf
  }else{
    FOUND_SUBCLONES <- FALSE
  }
  
  if(!FOUND_SUBCLONES) plotCNAlineOnlyTumor(sample)
  
  return(classDf)
  
}


DEBUGsubclonesAnalysis <- function(count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf){   
  
  res_subclones <- subcloneAnalysisPipeline(count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf)
  FOUND_SUBCLONES <- res_subclones$FOUND_SUBCLONES
  classDf <- res_subclones$classDf
  
  return(classDf)
  
}




