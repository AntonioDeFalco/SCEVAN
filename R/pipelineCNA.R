#' CNAvegaMC: A package for
#'
#' The CNAvegaMC package
#' 
#' @section CNAvegaMC functions:
#' The CNAvegaMC functions ...
#'
#' @docType package
#' @name CNAvegaMC
#' @useDynLib CNAvegaMC, .registration=TRUE
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
  
  res_proc <- preprocessingMtx(count_mtx, par_cores=par_cores, SMOOTH = TRUE)
  norm.cell <- getConfidentNormalCells(res_proc$count_mtx_smooth, par_cores = par_cores)
  norm.cell <- names(norm.cell)
  print(table(gr_truth[norm.cell]))
  
  res_class <- classifyTumorCells(res_proc$count_mtx_smooth,res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = gr_truth,  norm.cell.names = norm.cell, SEGMENTATION_CLASS = TRUE)
  
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
      
        diffSubcl <- testSpecificAlteration(res_proc$count_mtx_smooth, res_proc$count_mtx_annot, diffSubcl, res_subclones$clustersSub, res_subclones$n_subclones, sample)
      
        vectAlt <- annoteBand(res_proc$count_mtx_annot,diffSubcl)
        
        if(length(vectAlt)>0){
        perc_cells_subclones <- table(res_subclones$clustersSub)/length(res_subclones$clustersSub)
        plotSubclonesFish(as.integer(perc_cells_subclones[1]*100),as.integer(perc_cells_subclones[2]*100), vectAlt$vectAlt1, vectAlt$vectAlt2, vectAlt$vectAltsh, sample)
        #plotUMAP(count_mtx,res_class$CNAmat, rownames(res_proc$count_mtx_smooth), res_final$predTumorCells, res_final$clusters_subclones, sample)
        classDf[names(res_subclones$clustersSub), "subclone"] <- res_subclones$clustersSub
  
        #mtx_tum <- count_mtx[rownames(res_proc$count_mtx_smooth),res_class$tum_cells]
        #mtx_tum_log <- log2(mtx_tum + 1)
        
        #genesDE(mtx_tum_log, res_subclones$clustersSub, sample)
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

