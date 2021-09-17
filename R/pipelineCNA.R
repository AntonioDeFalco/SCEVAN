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
#' @param gr_truth ground truth of classification (optional)
#' @param SUBCLONES find subclones (optional)
#'
#' @return
#' @export
#'
#' @examples res_pip <- pipelineCNA(count_mtx, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)
pipelineCNA <- function(count_mtx, sample="", par_cores = 20,  gr_truth = NULL, SUBCLONES = TRUE){
  
  dir.create(file.path("./output"), showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  res_proc <- preprocessingMtx(count_mtx, par_cores=par_cores)
  #norm.cell <- getConfidentNormalCells(res_proc$count_mtx_norm, par_cores = par_cores)
  norm.cell <- names(res_proc$norm.cell)
  print(table(gr_truth[norm.cell]))
  
  res_class <- classifyTumorCells(res_proc$count_mtx_norm,res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = gr_truth,  norm.cell.names = norm.cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)
  
  #mtx_vega <- cbind(annot_mtx[,c(4,1,3)], count_mtx_relat)
  #colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
  #breaks <- getBreaksVegaMC(mtx_vega, annot_mtx[, 3], paste0(sample,"onlytumor"))
  
  mtx_vega <- cbind(res_class$CNAmat[,1:3], res_class$CNAmat[,res_class$tum_cells])
  colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
  breaks_tumor <- getBreaksVegaMC(mtx_vega, res_proc$count_mtx_annot[,3], paste0(sample,"onlytumor"))
  
  classDf <- data.frame(class = rep("filtered", length(colnames(count_mtx))), row.names = colnames(count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  
  end_time<- Sys.time()
  print(paste("time classify tumor cells: ", end_time -start_time))
  
  res_final <- list(res_class$confidentNormal, res_class$tum_cells)
  names(res_final) <- c("confidentNormalCells", "predTumorCells")
  
  if(length(gr_truth)>0){
    ground_truth_mal <- names(gr_truth[gr_truth == "malignant"])
    pred_mal <- res_class$tum_cells
    F1_Score <- computeF1score(pred_mal, ground_truth_mal)
    print(paste("F1_Score: ", F1_Score))
    res_final <- append(res_final, F1_Score)
    names(res_final)[3] <- "F1_Score"
  }
  
  start_time <- Sys.time()
  
  if (SUBCLONES){
    res_subclones <- subclonesTumorCells(res_class$tum_cells, res_class$CNAmat, sample)
    res_final <- append(res_final, list(res_subclones$n_subclones,res_subclones$breaks_subclones,res_subclones$clustersSub))
    names(res_final)[(length(names(res_final))-2):length(names(res_final))] <- c("n_subclones", "breaks_subclones", "clusters_subclones")
    
    tum_cells <- res_class$tum_cells
    clustersSub <- res_subclones$clustersSub
    save(tum_cells,clustersSub, file = paste0(sample,"subcl.RData"))
    
    if(res_subclones$n_subclones>1){
    sampleAlter <- analyzeSegm(sample, nSub = res_subclones$n_subclones)
    
      if(length(sampleAlter)>1){
      
        diffSubcl <- diffSubclones(sampleAlter, sample, nSub = res_subclones$n_subclones)
         
        diffSubcl <- testSpecificAlteration(res_class$CNAmat, res_proc$count_mtx_annot, diffSubcl, res_subclones$clustersSub, res_subclones$n_subclones, sample)
        
        print(diffSubcl)
        
        segmList <- list()
        segmList$subclone1 <- read.table(paste0("./output/ ",sample,"_subclone1 vega_output"), sep="\t", header=TRUE, as.is=TRUE)
        segmList$subclone2 <- read.table(paste0("./output/ ",sample,"_subclone2 vega_output"), sep="\t", header=TRUE, as.is=TRUE)

        save(res_proc, res_subclones, segmList,diffSubcl,sample, file = "mgh125_plotcnaline2.RData")
        plotCNAline(segmList, diffSubcl, sample)
        
        diffSubcl[[grep("_clone",names(diffSubcl))]] <- diffSubcl[[grep("_clone",names(diffSubcl))]][1:min(10,nrow(diffSubcl[[grep("_clone",names(diffSubcl))]])),]
        
        vectAlt <- annoteBand(res_proc$count_mtx_annot,diffSubcl)
        
        if(length(vectAlt)>0){
            perc_cells_subclones <- table(res_subclones$clustersSub)/length(res_subclones$clustersSub)
            plotSubclonesFish(as.integer(perc_cells_subclones[1]*100),as.integer(perc_cells_subclones[2]*100), vectAlt[[1]], vectAlt[[2]], vectAlt[[3]], sample)
            plotUMAP(count_mtx,res_class$CNAmat, rownames(res_proc$count_mtx_norm), res_final$predTumorCells, res_final$clusters_subclones, sample)
            classDf[names(res_subclones$clustersSub), "subclone"] <- res_subclones$clustersSub
      
            genesDE(res_proc$count_mtx_norm, res_proc$count_mtx_annot, res_subclones$clustersSub, sample, diffSubcl[grep("subclone",names(diffSubcl))])
            pathwayAnalysis(res_proc$count_mtx_norm, res_proc$count_mtx_annot, res_subclones$clustersSub, sample)
        }
      }else{
      print("no significant subclones")
    }
    
    }
    
   }
  
  end_time<- Sys.time()
  print(paste("time subclones: ", end_time -start_time))

  return(classDf)
}


#' compareClonalStructure
#'
#' @param count_mtx raw count matrix
#' @param sample sample name (optional)
#' @param par_cores number of cores (optional)
#' @param gr_truth ground truth of classification (optional)
#' @param SUBCLONES find subclones (optional)
#'
#' @return
#' @export
#'
#' @examples 
#' 
compareClonalStructure <- function(count_mtx1, count_mtx2 , samp_1="", samp_2="", par_cores = 20){

  
  count_mtx <- merge(count_mtx1, count_mtx2, by="row.names" ,all = TRUE)
  count_mtx[is.na(count_mtx)] <- 0
  rownames(count_mtx) <- count_mtx$Row.names
  count_mtx<- count_mtx[,-1]
  
  res_proc_1 <- preprocessingMtx(count_mtx1, par_cores=par_cores)
  norm.cell <- names(res_proc_1$norm.cell)
  res_class_1 <- classifyTumorCells(res_proc_1$count_mtx_norm,res_proc_1$count_mtx_annot, samp_1, par_cores=par_cores,  norm.cell.names = norm.cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)
  
  res_proc_2 <- preprocessingMtx(count_mtx2, par_cores=par_cores)
  norm.cell <- names(res_proc_2$norm.cell)
  res_class_2 <- classifyTumorCells(res_proc_2$count_mtx_norm,res_proc_2$count_mtx_annot, samp_2, par_cores=par_cores,  norm.cell.names = norm.cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE)
  
  sampl <- paste0(samp_1,"-vs-", samp_2)
  
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
  
  diffSubcl <- testSpecificAlteration(res_class$CNAmat, res_proc$count_mtx_annot, all_segm, clust_subclones, nSub = 2, sampl)
  
  plotCNAline(segmList, diffSubcl, sampl)
  
  diffSubcl[[grep("_clone",names(diffSubcl))]] <- diffSubcl[[grep("_clone",names(diffSubcl))]][1:min(10,nrow(diffSubcl[[grep("_clone",names(diffSubcl))]])),]
  
  vectAlt <- annoteBand(res_proc$count_mtx_annot,diffSubcl)
  
  if(length(vectAlt)>0){
    perc_cells_subclones <- table(clust_subclones)/length(clust_subclones)
    plotSubclonesFish(as.integer(perc_cells_subclones[1]*100),as.integer(perc_cells_subclones[2]*100), vectAlt[[1]], vectAlt[[2]], vectAlt[[3]], sampl)
  
    res_final <- list()
    res_final$predTumorCells <- append(res_class_1$tum_cells,res_class_2$tum_cells)
    plotUMAP(count_mtx,res_class$CNAmat, rownames(res_proc$count_mtx_norm), res_final$predTumorCells, clust_subclones, sampl)

    genesDE(res_proc$count_mtx_norm, res_proc$count_mtx_annot, clust_subclones, sampl, diffSubcl[grep("subclone",names(diffSubcl))])
    pathwayAnalysis(res_proc$count_mtx_norm, res_proc$count_mtx_annot, clust_subclones, sampl)
  }
  
}