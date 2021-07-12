
#' estimate basline copy numbers by synthetic normal cell
#'
#' @param norm.mat smoothed data matrix; genes in rows; cell names in columns.
#' @param min.cells minimal number of cells per cluster.
#' @param par_cores number of cores for parallel computing.
#'
#' @return 1) relative gene expression; 2) synthetic baseline profiles; 3) clustering results.
#'
#' @examples
#'
#' count_mtx_relat <- test.relt$expr.relat
#' @export
removeSyntheticBaseline <- function(count_mtx, par_cores = 20){ 
  
  d <- parallelDist::parDist(t(count_mtx), threads = par_cores) 
  
  hcc <- hclust(d, method="ward.D2")
  
  sCalinsky <- calinsky(hcc, d, gMax = 10)
  k <- which.max(sCalinsky)
  hcc_k <- cutree(hcc, k=k)
  
  expr.relat <- NULL
  for(i in 1:k){
    cellsClu <- count_mtx[, which(hcc_k==i)]
    sdCells <- apply(cellsClu,1,sd)
    set.seed(123)
    syn.norm <- sapply(sdCells,function(x)(x<- rnorm(1,mean = 0,sd=x)))
    cellsCluRelat <- cellsClu - syn.norm
    expr.relat <- rbind(expr.relat, t(cellsCluRelat))
  }
  
  return(t(expr.relat))
}



computeCNAmtx <- function(mtx, breaks, par_cores){
  
  n <- nrow(mtx)
  
  seg.test <- parallel::mclapply(1:ncol(mtx), function(z){
    x<-numeric(n)
    for (i in 1:(length(breaks)-1)){
      x[breaks[i]:breaks[i+1]] <- mean(mtx[breaks[i]:breaks[i+1],z])
    }
    return(x)
  }, mc.cores = par_cores)
  
  CNA <- matrix(unlist(seg.test), ncol = ncol(mtx), byrow = FALSE)
  
  return(CNA)
}

classifyCluster <- function(hcc2, norm.cell.names){
  perc_norm <- length(intersect(names(hcc2[hcc2==1]), norm.cell.names))/length(hcc2[hcc2==1])
  perc_norm <- c(perc_norm,length(intersect(names(hcc2[hcc2==2]), norm.cell.names))/length(hcc2[hcc2==2]))
  clust_norm <- which(perc_norm==max(perc_norm))
  
  cellType_pred <- names(hcc2)
  cellType_pred[hcc2 == clust_norm] <- "non malignant"
  cellType_pred[!hcc2 == clust_norm] <- "malignant"
  names(cellType_pred) <- names(hcc2)
  
  return(cellType_pred)
}



#' classifyTumorCells classify tumour and normal cells from the raw count matrix
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
classifyTumorCells <- function(count_mtx_smooth, count_mtx_annot, sample, distance="euclidean", par_cores=20, ground_truth = NULL, norm.cell.names = NULL, SEGMENTATION_CLASS = TRUE){
  
  set.seed(1)
  
  if (length(norm.cell.names) < 1){
    print("8): measuring baselines (pure tumor - synthetic normal cells)")
    count_mtx_relat <- removeSyntheticBaseline(count_mtx_smooth, par_cores=par_cores)
    #count_mtx_relat <- relt$expr.relat
    #colnames(count_mtx_relat) <- colnames(relt$expr.relat)
    
  } else {
    print("8): measuring baselines (confident normal cells)")
    
    if(length(norm.cell.names) == 1){
      basel <- count_mtx_smooth[, which(colnames(count_mtx_smooth) %in% norm.cell.names)]
    }else{
      basel <- apply(count_mtx_smooth[, which(colnames(count_mtx_smooth) %in% norm.cell.names)],1,median)
    }
    
    ##relative expression using pred.normal cells
    count_mtx_relat <- count_mtx_smooth-basel
    
  }

  ##### Segmentation with VegaMC #####
  
  if(SEGMENTATION_CLASS & length(norm.cell.names) > 0){
    
    print("9) Segmentation (VegaMC)")
    
    mtx_vega <- cbind(count_mtx_annot[,c(4,1,3)], count_mtx_relat)
    colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
    
    breaks <- getBreaksVegaMC(mtx_vega, count_mtx_annot[, 3], sample)
    
    CNA_mtx <- computeCNAmtx(count_mtx_relat, breaks, par_cores)
    SEGM <- TRUE
    
  }else{
    SEGM <- FALSE
    CNA_mtx <- count_mtx_relat
  }
  
  colnames(CNA_mtx) <- colnames(count_mtx_relat)
  CNA_mtx <- apply(CNA_mtx,2, function(x)(x <- x-mean(x)))
  
  print("10) Adjust baseline")
  
  if(length(norm.cell.names) < 1){
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx, method = distance)), method = "ward.D")
    }
    
    #plot heatmap
    print("11) plot heatmap")
    
    plotCNA(count_mtx_annot$seqnames, CNA_mtx, hcc, sample)
    

  } else {
    
    ################removed baseline adjustment
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx, method = distance)), method = "ward.D")
    }
    
    hcc2 <- cutree(hcc,2)
    names(hcc2) <- colnames(CNA_mtx)
    
    cellType_pred <- classifyCluster(hcc2, norm.cell.names)
    
    ################removed baseline adjustment
    CNA_mtx_relat <- CNA_mtx-apply(CNA_mtx[,which(cellType_pred=="non malignant")], 1, mean)
    CNA_mtx_relat <- apply(CNA_mtx_relat,2,function(x)(x <- x-mean(x)))
    CNA_mtx_relat <- CNA_mtx_relat/(0.5*(max(CNA_mtx_relat)-min(CNA_mtx_relat)))
    
    if(SEGM){
    count_mtx_relat <- count_mtx_relat-apply(count_mtx_relat[,which(cellType_pred=="non malignant")], 1, mean)
    count_mtx_relat <- apply(count_mtx_relat,2,function(x)(x <- x-mean(x)))
    count_mtx_relat <- count_mtx_relat/(0.5*(max(count_mtx_relat)-min(count_mtx_relat)))
    }
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx_relat),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx_relat, method = distance)), method = "ward.D")
    }
    
    hcc2 <- cutree(hcc,2)
    names(hcc2) <- colnames(CNA_mtx_relat)
    
    cellType_pred <- classifyCluster(hcc2, norm.cell.names)
    
    res <- cbind(names(cellType_pred), cellType_pred)
    colnames(res) <- c("cell.names", "pred")
    
    print("11) plot heatmap")
    
    plotCNA(count_mtx_annot$seqnames, CNA_mtx_relat, hcc, sample, cellType_pred, ground_truth)
    
  }
  if(length(norm.cell.names) < 1){
    tum_cells <- colnames(CNA_mtx)
    tum_cells <- gsub("\\.","-",tum_cells)
  }else{
    tum_cells <- names(res[,2][res[,2] == "malignant"])
  }
  
  if(SEGM){
    ress <- list(tum_cells, cbind(count_mtx_annot[,c(4,1,3)], count_mtx_relat), norm.cell.names)
  }else if(length(norm.cell.names) < 1){
    ress <- list(tum_cells, cbind(count_mtx_annot[,c(4,1,3)], CNA_mtx), norm.cell.names)
  }else{
    ress <- list(tum_cells, cbind(count_mtx_annot[,c(4,1,3)], CNA_mtx_relat), norm.cell.names)
  }
  
  names(ress) <- c("tum_cells", "CNAmat", "confidentNormal")
  return(ress)
}


