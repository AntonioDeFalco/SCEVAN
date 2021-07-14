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
  
  if (SUBCLONES){
    res_subclones <- subclonesTumorCells(res_class$tum_cells, res_class$CNAmat, sample)
    res_final <- append(res_final, list(res_subclones$n_subclones,res_subclones$breaks_subclones,res_subclones$clustersSub))
    names(res_final)[(length(names(res_final))-2):length(names(res_final))] <- c("n_subclones", "breaks_subclones", "clusters_subclones")
    
    if(res_subclones$n_subclones>1){
    sampleAlter <- analyzeSegm(sample, nSub = res_subclones$n_subclones)
    diffSubclones <- diffSubclones(sampleAlter, sample, nSub = res_subclones$n_subclones)
    print(diffSubclones)
    
    perc_cells_subclones <- table(res_subclones$clustersSub)/length(res_subclones$clustersSub)
    
    vectAlt1 <- apply(diffSubclones[[1]], 1, function(x){ gsub(" ","",x); paste0(x[1],x[4])})
    vectAlt1 <- do.call(paste, c(as.list(unique(vectAlt1)), sep = "\n"))
    
    vectAlt2 <- apply(diffSubclones[[2]], 1, function(x){ gsub(" ","",x); paste0(x[1],x[4])})
    vectAlt2 <- do.call(paste, c(as.list(unique(vectAlt2)), sep = "\n"))
    
    vectAltsh <- apply(diffSubclones[[3]], 1, function(x){ gsub(" ","",x); paste0(x[1],x[4])})
    vectAltsh <- do.call(paste, c(as.list(unique(vectAltsh)), sep = "\n"))
    
    plotSubclonesFish(as.integer(perc_cells_subclones[1]*100),as.integer(perc_cells_subclones[2]*100),vectAlt1, vectAlt2, vectAltsh, sample)
    
    }
  }
  
  end_time<- Sys.time()
  print(end_time -start_time)
  
  return(res_final)
}
