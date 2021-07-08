
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
  
  sCalinsky <- calinsky(fit, d, gMax = 15)
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
    #i <- i+1
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
      i<- i+1
    }
    return(x)
  }, mc.cores = par_cores)
  
  CNA <- matrix(unlist(seg.test), ncol = ncol(mtx), byrow = FALSE)
  
  return(CNA)
}



#' func
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
classifyTumorCells <- function(count_mtx_smooth, count_mtx_proc, count_mtx_annot, sample, distance="euclidean", par_cores=20, ground_truth = NULL, norm.cell.names = NULL, UP.DR = 0.1, ngene.chr=5, WRITE = FALSE, SEGMENTATION_CLASS = TRUE){
  
  set.seed(1)
  
  if (length(norm.cell.names) < 1){
    print("7): measuring baselines (pure tumor - synthetic normal cells)")
    relt <- baseline.synthetic(norm.mat=count_mtx_smooth, par_cores=par_cores)
    count_mtx_relat <- relt$expr.relat
    colnames(count_mtx_relat) <- colnames(relt$expr.relat)
    CL <- relt$cl
    
  } else {
    print("7): measuring baselines (confident normal cells)")
    
    basel <- apply(count_mtx_smooth[, which(colnames(count_mtx_smooth) %in% norm.cell.names)],1,median)
    
    d <- parallelDist::parDist(t(count_mtx_smooth),threads =par_cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells
    
    fit <- hclust(d, method="ward.D2")
    #CL <- cutree(fit, km)
    
    sCalinsky <- calinsky(fit, d, gMax = 15)
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
  
  
  ###use a smaller set of genes to perform segmentation
  DR2 <- apply(count_mtx_proc,1,function(x)(sum(x>0)))/ncol(count_mtx_proc)
  ##relative expression using pred.normal cells
  count_mtx_relat <- count_mtx_relat[which(DR2>=UP.DR),]
  
  ###filter cells
  anno.mat2 <- count_mtx_annot[which(DR2>=UP.DR), ]
  
  ToRemov3 <- NULL
  for(i in 6:ncol(anno.mat2)){
    cell <- cbind(anno.mat2$seqnames, anno.mat2[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    } else if(length(rle(cell[,1])$length)<22|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i<- i+1
  }
  
  if(length(ToRemov3)==ncol(count_mtx_relat)) stop ("all cells are filtered")
  
  if(length(ToRemov3)>0 & !length((which(colnames(count_mtx_relat) %in% ToRemov3))) == 0L){
    count_mtx_relat <- count_mtx_relat[, -which(colnames(count_mtx_relat) %in% ToRemov3)]
    print(paste("filtered out ", length(ToRemov3), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  }
  
  print(paste("final segmentation: ", nrow(count_mtx_relat), " genes; ", ncol(count_mtx_relat), " cells", sep=""))
  
  CL <- CL[which(names(CL) %in% colnames(count_mtx_relat))]
  CL <- CL[order(match(names(CL), colnames(count_mtx_relat)))]
  
  
  ##### Segmentation with VegaMC #####
  
  if(SEGMENTATION_CLASS & length(norm.cell.names) > 0){
    
    print("8) Segmentation (VegaMC)")
    
    mtx_vega <- cbind(anno.mat2[,c(4,1,3)], count_mtx_relat)
    colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
    
    BR <- getBreaksVegaMC(mtx_vega, anno.mat2[, 3], sample)
    
    logCNA <- computeCNAmtx(count_mtx_relat, BR, par_cores)
    
    results <- list(logCNA, BR)
    names(results) <- c("logCNA","breaks")
    
  }else{
    results <- list(count_mtx_relat, c())
    names(results) <- c("logCNA","breaks")
  }
  
  colnames(results$logCNA) <- colnames(count_mtx_relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:5], results.com)
  
  if(WRITE){
    write.table(RNA.copycat, paste(sample, "CNA_raw_results_gene_by_cell.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
  }
  
  print("8) Adjust baseline")
  
  if(length(norm.cell.names) < 1){
    
    mat.adj <- data.matrix(RNA.copycat[6:ncol(RNA.copycat)])
    
    if(WRITE){
      write.table(cbind(RNA.copycat[,c(5,2,1)], mat.adj), paste(sample, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
    }
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
    }
    
    if(WRITE){
      saveRDS(hcc, file = paste(sample,"clustering_results.rds",sep=""))
    }
    
    #plot heatmap
    print("9) plot heatmap")
    
    plotCNA(RNA.copycat$seqnames, mat.adj, hcc, sample)
    

    
  } else {
    
    
    uber.mat.adj <- data.matrix(RNA.copycat[6:ncol(RNA.copycat)])
    
    ################removed baseline adjustment
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
    }
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)
    
    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, norm.cell.names))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }
    
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "non malignant"
    com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "malignant"
    names(com.pred) <- names(hc.umap)
    
    ################removed baseline adjustment
    results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="non malignant")], 1, mean)
    results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
    results.com.rat.norm <- results.com.rat[,which(com.pred=="non malignant")]; dim(results.com.rat.norm)
    
    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)
    
    adjN <- function(j){
      a <- results.com.rat[, j]
      a[abs(a-base) <= 0.25*cf.h] <- mean(a)
      a
    }
    
    
    mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = par_cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
    rm(mc.adjN)
    colnames(adj.results) <- colnames(results.com.rat)
    
    rang <- 0.5*(max(adj.results)-min(adj.results))
    mat.adj <- adj.results/rang
    
    print("step 8: final prediction ...")
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
    }
    
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)
    
    if(WRITE){
      saveRDS(hcc, file = paste(sample,"clustering_results.rds",sep=""))
    }
    
    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, norm.cell.names))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }
    
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "non malignant"
    com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "malignant"
    names(com.preN) <- names(hc.umap)
    
    print("step 9: saving results...")
    res <- cbind(names(com.preN), com.preN)
    colnames(res) <- c("cell.names", "copykat.pred")
    
    if(WRITE){
      write.table(res, paste(sample, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
      
      ####save copycat CNA
      write.table(cbind(RNA.copycat[,c(5,2,1)], mat.adj), paste(sample, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
    }
    
    ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
    print("step 10: ploting heatmap ...")
    
    plotCNA(RNA.copycat$seqnames, mat.adj, hcc, sample, com.preN, ground_truth)
    
  }
  
  
  if(length(norm.cell.names) < 1){
    tum_cells <- colnames(RNA.copycat[,-c(1:5)])
    tum_cells <- gsub("\\.","-",tum_cells)
  }else{
    tum_cells <- names(res[,2][res[,2] == "malignant"])
  }
  
  ress <- list(tum_cells, cbind(RNA.copycat[,c(4,1,3)], mat.adj), norm.cell.names)
  
  names(ress) <- c("tum_cells", "CNAmat", "confidentNormal")
  return(ress)
}