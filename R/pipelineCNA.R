#' Run pipeline
#'
#' @param count_mtx raw count matrix
#' @param par_cores number of cores
#' @param gr_truth ground truth of classification
#' @param SUBCLONES find subclones
#'
#' @return
#' @export
#'
#' @examples res_pip <- pipelineCNA(count_mtx, par_cores = 20, gr_truth = gr_truth, SUBCLONES = TRUE)
pipelineCNA <- function(count_mtx, sample="", par_cores = 20,  gr_truth = NULL, SUBCLONES = TRUE){
  
  start_time <- Sys.time()
  
  res_proc <- preprocessingMtx(count_mtx, par_cores=par_cores, SMOOTH = TRUE)
  norm.cell <- getConfidentNormalCells(res_proc$count_mtx_smooth, par_cores = par_cores)
  norm.cell <- names(norm.cell)
  print(table(gr_truth[norm.cell]))
  
  res_class <- classifyTumorCells(res_proc$count_mtx_smooth, res_proc$count_mtx_proc,res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = gr_truth,  norm.cell.names = norm.cell, SEGMENTATION_CLASS = TRUE)
  
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
    res_final <- append(res_final, list(res_subclones$n_subclones,res_subclones$breaks_subclones))
    names(res_final)[4:5] <- c("n_subclones", "breaks_subclones")
  }
  
  end_time<- Sys.time()
  print(end_time -start_time)
  
  return(res_final)
}
