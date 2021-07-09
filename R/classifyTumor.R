
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
baseline.synthetic <- function(norm.mat=norm.mat, min.cells=10, par_cores = 20){ 
  
  d <- parallelDist::parDist(t(norm.mat), threads = par_cores) 
  
  fit <- hclust(d, method="ward.D2")
  
  sCalinsky <- calinsky(fit, d, gMax = 10)
  km <- which.max(sCalinsky)
  #km <- 6
  ct <- cutree(fit, k=km)
  
  while(!all(table(ct)>min.cells)){
    km <- km -1
    ct <- cutree(fit, k=km)
    if(km==2){
      break
    }
  }
  
  expr.relat <- NULL
  syn <- NULL
  for(i in min(ct):max(ct)){
    data.c1 <- norm.mat[, which(ct==i)]
    sd1 <- apply(data.c1,1,sd)
    set.seed(123)
    syn.norm <- sapply(sd1,function(x)(x<- rnorm(1,mean = 0,sd=x)))
    relat1 <- data.c1 - syn.norm
    expr.relat <- rbind(expr.relat, t(relat1))
    syn <- cbind(syn,syn.norm)
  }
  
  reslt <- list(data.frame(t(expr.relat)), data.frame(syn), ct)
  names(reslt) <- c("expr.relat","syn.normal", "cl")
  
  return(reslt)
}



computeCNAmtx <- function(mtx, BR, par_cores){
  
  n <- nrow(mtx)
  
  seg.test <- parallel::mclapply(1:ncol(mtx), function(z){
    x<-numeric(n)
    for (i in 1:(length(BR)-1)){
      x[BR[i]:BR[i+1]] <- mean(mtx[BR[i]:BR[i+1],z])
    }
    return(x)
  }, mc.cores = par_cores)
  
  CNA <- matrix(unlist(seg.test), ncol = ncol(mtx), byrow = FALSE)
  
  return(CNA)
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
    relt <- baseline.synthetic(norm.mat=count_mtx_smooth, par_cores=par_cores)
    count_mtx_relat <- relt$expr.relat
    colnames(count_mtx_relat) <- colnames(relt$expr.relat)
    CL <- relt$cl
    
  } else {
    print("8): measuring baselines (confident normal cells)")
    
    if(length(norm.cell.names) == 1){
      basel <- count_mtx_smooth[, which(colnames(count_mtx_smooth) %in% norm.cell.names)]
    }else{
      basel <- apply(count_mtx_smooth[, which(colnames(count_mtx_smooth) %in% norm.cell.names)],1,median)
    }
    
    d <- parallelDist::parDist(t(count_mtx_smooth),threads =par_cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells
    
    fit <- hclust(d, method="ward.D2")
    #CL <- cutree(fit, km)
    
    sCalinsky <- calinsky(fit, d, gMax = 10)
    km <- which.max(sCalinsky)
    #km <- 6
    print(paste("cluster calisky:", km))
    
    CL <- cutree(fit, km)
    
    while(!all(table(CL)>5)){
      km <- km -1
      CL <- cutree(fit, k=km)
      if(km==2){
        break
      }
    }
    
    ##relative expression using pred.normal cells
    count_mtx_relat <- count_mtx_smooth-basel
    
  }
  
  print(paste("final segmentation: ", nrow(count_mtx_relat), " genes; ", ncol(count_mtx_relat), " cells", sep=""))
  
  CL <- CL[which(names(CL) %in% colnames(count_mtx_relat))]
  CL <- CL[order(match(names(CL), colnames(count_mtx_relat)))]
  
  
  ##### Segmentation with VegaMC #####
  
  if(SEGMENTATION_CLASS & length(norm.cell.names) > 0){
    
    print("9) Segmentation (VegaMC)")
    
    mtx_vega <- cbind(count_mtx_annot[,c(4,1,3)], count_mtx_relat)
    colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
    
    BR <- getBreaksVegaMC(mtx_vega, count_mtx_annot[, 3], sample)
    
    CNA_mtx <- computeCNAmtx(count_mtx_relat, BR, par_cores)
    
  }else{
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
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(CNA_mtx)
    
    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, norm.cell.names))/length(cli)
      cl.ID <- c(cl.ID, pid)
    }
    
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "non malignant"
    com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "malignant"
    names(com.pred) <- names(hc.umap)
    
    ################removed baseline adjustment
    CNA_mtx_relat <- CNA_mtx-apply(CNA_mtx[,which(com.pred=="non malignant")], 1, mean)
    CNA_mtx_relat <- apply(CNA_mtx_relat,2,function(x)(x <- x-mean(x)))
    
    cf.h <- apply(CNA_mtx_relat[,which(com.pred=="non malignant")], 1, sd)
    base <- apply(CNA_mtx_relat[,which(com.pred=="non malignant")], 1, mean)
    
    adjN <- function(j){
      a <- CNA_mtx_relat[, j]
      a[abs(a-base) <= 0.25*cf.h] <- mean(a)
      a
    }
    
    mc.adjN <-  parallel::mclapply(1:ncol(CNA_mtx_relat),adjN, mc.cores = par_cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(CNA_mtx_relat), byrow = FALSE)
    rm(mc.adjN)
    colnames(adj.results) <- colnames(CNA_mtx_relat)
    
    rang <- 0.5*(max(adj.results)-min(adj.results))
    CNA_mtx <- adj.results/rang
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx, method = distance)), method = "ward.D")
    }
    
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(CNA_mtx)
    
    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, norm.cell.names))/length(cli)
      cl.ID <- c(cl.ID, pid)
    }
    
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "non malignant"
    com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "malignant"
    names(com.preN) <- names(hc.umap)
    
    res <- cbind(names(com.preN), com.preN)
    colnames(res) <- c("cell.names", "pred")
    
    print("11) plot heatmap")
    
    plotCNA(count_mtx_annot$seqnames, CNA_mtx, hcc, sample, com.preN, ground_truth)
    
  }
  
  
  if(length(norm.cell.names) < 1){
    tum_cells <- colnames(CNA_mtx)
    tum_cells <- gsub("\\.","-",tum_cells)
  }else{
    tum_cells <- names(res[,2][res[,2] == "malignant"])
  }
  
  ress <- list(tum_cells, cbind(count_mtx_annot[,c(4,1,3)], CNA_mtx), norm.cell.names)
  
  names(ress) <- c("tum_cells", "CNAmat", "confidentNormal")
  return(ress)
}