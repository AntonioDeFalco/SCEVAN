



#' Title
#'
#' @param pred 
#' @param ground_truth 
#'
#' @return
#' @export
#'
#' @examples
computeF1score <- function(pred, ground_truth){
  
  TP <- length(intersect(pred,ground_truth))
  FP <- sum(!(pred %in% ground_truth))
  FN <- sum(!(ground_truth %in% pred))
  
  F1_Score <- TP/(TP + (1/2)*(FP+ FN))
  
  return(F1_Score)
}

calinsky <- function (hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order),
                                                                   10))){
  msg <- ""
  if (is.null(dist)) {
    require(clue)
    dist <- sqrt(as.cl_ultrametric(hhc))
    message(msg <- "The distance matrix not is provided, using the cophenetic matrix")
  }
  else if (attr(dist, "method") != "euclidean") {
    require(clue)
    dist <- sqrt(as.cl_ultrametric(hhc))
    message(msg <- "The distance matrix is not euclidean, using the cophenetic matrix")
  }
  dist <- as.matrix(dist)^2
  A <- -dist/2
  A_bar <- apply(A, 1, mean)
  totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
  n <- length(hhc$order)
  ans <- rep(0, gMax)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1)
        next
      A <- as.matrix(-dist/2)[cclust == k, cclust == k]
      A_bar <- apply(A, 1, mean)
      withinSum <- withinSum + sum(diag(A) - 2 * A_bar +
                                     mean(A))
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  class(ans) <- "calinsky"
  attr(ans, "message") <- msg
  return(ans)
}














library(cluster)

clusGap_parallel <- function (x, FUNcluster, K.max, B = 100, d.power = 1,
                              n.cores = 20, spaceH0 = c("scaledPCA", 
                                                        "original"), verbose = interactive(), ...) 
{
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 
              2, (n <- nrow(x)) >= 1, ncol(x) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
    stop("'B' has to be a positive integer")
  cl. <- match.call()
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
    clus <- if (kk > 1) 
      FUNcluster(X, kk, ...)$cluster
    else rep.int(1L, nrow(X))
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(dist(xs)^d.power/nrow(xs))
    }, 0))
  }
  logW <- E.logW <- SE.sim <- numeric(K.max)
  if (verbose) 
    cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. ", 
        sep = "")
  logW <- unlist(parallel::mclapply(1:K.max, function(k) log(W.k(x, k)), mc.cores = 4))
  if (verbose) 
    cat("done\n")
  spaceH0 <- match.arg(spaceH0)
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  switch(spaceH0, scaledPCA = {
    V.sx <- svd(xs, nu = 0)$v
    xs <- xs %*% V.sx
  }, original = {
  }, stop("invalid 'spaceH0':", spaceH0))
  rng.x1 <- apply(xs, 2L, range)
  logWks <- matrix(0, B, K.max)
  if (verbose) 
    cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
        sep = "")
  for (b in 1:B) {
    z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
                                                 max = M[2]), nn = n)
    z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), 
                original = z1) + m.x
    tmplogWks <- unlist(parallel::mclapply(1:K.max, function(k) log(W.k(z, k)), mc.cores = 4))
    
    logWks[b,1:K.max] <- tmplogWks
    if (verbose) 
      cat(".", if (b%%50 == 0) 
        paste(b, "\n"))
  }
  if (verbose && (B%%50 != 0)) 
    cat("", B, "\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  structure(class = "clusGap", list(Tab = cbind(logW, E.logW, 
                                                gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0, 
                                    n = n, B = B, FUNcluster = FUNcluster))
}


k_clusGap <- function(mtx){
  
  mydist <- function(x){
    
    return (parallelDist::parDist(t(x),threads = 20, method = "euclidean"))
  }
  
  hcc <- hclust(mydist(mtx), method = "ward.D")
  
  mycluster <- function(x, k) list(cluster=cutree(hcc, k=k))
  myclusGap <- clusGap_parallel(mtx,
                                FUN = mycluster, 
                                K.max = 5, 
                                B = 50)
  plot(myclusGap)
  k <- maxSE(f         = myclusGap$Tab[,"gap"],
             SE.f      = myclusGap$Tab[,"SE.sim"],
             method    = "firstSEmax",
             SE.factor = 1)
  return (k)
  
}










#' Title
#'
#' @param tum_cells 
#' @param CNAmat 
#' @param sample 
#'
#' @return
#' @export
#'
#' @examples
subclonesTumorCells <- function(tum_cells, CNAmat, sample){
  
  norm.mat.relat <- CNAmat[,-c(1:3)]
  info_mat <- CNAmat[,c(1:3)]
  
  if(grep(".",(tum_cells)[1]) & grep("-",colnames(norm.mat.relat)[1])){
    tum_cells <- gsub("\\.", "-",tum_cells)
  }else if((grep("-",(tum_cells)[1]) & grep(".",colnames(norm.mat.relat)[1]))){
    tum_cells <- gsub("-", "\\.",tum_cells)
  }
  
  norm.mat.relat <- norm.mat.relat[,tum_cells]
  
  #k <- k_clusGap(norm.mat.relat)
  
  n.cores = 20 
  distance="euclidean"
  dist_mtx <- parallelDist::parDist(t(norm.mat.relat),threads = n.cores, method = distance)
  
  hcc <- hclust(dist_mtx, method = "ward.D")
  #plot(hcc, cex = 0.6, hang = -1)
  hc.clus <- cutree(hcc, h = 15)
  
  n_subclones <- 0
  results.com <- NULL
  breaks_subclones <- NULL
  
  if(length(hc.clus) > 1){
    
    sCalinsky <- calinsky(hcc, dist_mtx, gMax = 10)
    #plot(sCalinsky, type= "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
    #text(1:length(sCalinsky), sCalinsky, paste(1:length(sCalinsky)))
    
    n_subclones <- which.max(sCalinsky)
    
    print(paste("sCalinsky Max", sCalinsky[which.max(sCalinsky)]))
    
    hc.clus <- cutree(hcc,n_subclones)
    print(paste("found", n_subclones, "subclones", sep = " "))
    perc_cells_subclones <- table(hc.clus)/length(hc.clus)
    names(perc_cells_subclones) <- paste0("percentage_cells_subsclone_",names(perc_cells_subclones))
    print(perc_cells_subclones)
    
    if(which.max(sCalinsky)>0.10 & all(perc_cells_subclones > 0.10)){
      
      breaks_subclones <- list()
      
      for (i in 1:n_subclones){
        
        print(paste("Segmentation of subclone : ", i))
        
        tum_cells_sub1 <- names(hc.clus[hc.clus==i])
        
        
        mtx_vega <- cbind(info_mat, norm.mat.relat[,tum_cells_sub1])
        
        colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
        
        breaks_subclones[[i]] <- getBreaksVegaMC(mtx_vega, CNAmat[,3], paste0(sample,"_subclone",i))
      }
      
      BR <- c()
      
      BR <- sort(unique(unlist(breaks_subclones)))
      
      logCNA <- computeCNAmtx(norm.mat.relat, BR, n.cores)
      
      results <- list(logCNA, BR)
      names(results) <- c("logCNA","breaks")
      
      colnames(results$logCNA) <- colnames(norm.mat.relat)
      results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
      
      rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:n_subclones])
      subclones <- rbPal5(2)[as.numeric(factor(hc.clus))]
      cells <- rbind(subclones,subclones)
      
      
      chr <- as.numeric(CNAmat[,2]) %% 2+1
      
      rbPal1 <- colorRampPalette(c('black','grey'))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR,CHR)
      
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
      col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
      
      jpeg(paste(sample,"heatmap_subclones.jpeg",sep=""), height=10*250, width=4000, res=100)
      heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
                ColSideColors=chr1,RowSideColors=cells, Colv=NA, Rowv=TRUE,
                notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
                keysize=1, density.info="none", trace="none",
                cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                symm=F,symkey=F,symbreaks=T,cex=1, main=paste("Subclones ", sample), cex.main=4, margins=c(10,10))
      dev.off()
    }
  }
  
  
  ress <- list(n_subclones, breaks_subclones, results.com)
  
  names(ress) <- c("n_subclones", "breaks_subclones", "logCNA")
  
  
  return(ress)
  
}






createGeneSetNormal <- function(){
  
  #ESTIMATE
  load("/storage/qnap_home/adefalco/singleCell/AllData/ClassTumorCells/SI_geneset.RData") #from https://sourceforge.net/projects/estimateproject/
  geneSet <- c()
  geneSet$stromal <- as.character(unlist(SI_geneset["StromalSignature",-1]))
  geneSet$immune <- as.character(unlist(SI_geneset["ImmuneSignature",-1]))
  
  #FRANCESCA & CANCER CELLS
  load("/storage/qnap_home/caruso/Analisi2021/signatures/scTHI_c8_signatures_968.RData")
  load("/storage/qnap_home/caruso/Analisi2021/CPTAC/FgesSignature/fges_signature.RData")
  
  findSign <- function(Phenotype){
    all_sign <- rownames(signature_Colors[signature_Colors$ALLPhenotypeFinal==Phenotype,])
    ind <- which(names(signature) %in% all_sign)
    subset_sign <- signature[ind]
    num_sin <- length(subset_sign)/2
    subset_sign <- unlist(subset_sign)
    
    i <- 1
    while( (length(unique(subset_sign)) > 200) & (i < num_sin)){
      subset_sign <- subset_sign[duplicated(subset_sign)]
      i <- i + 1
    }
    
    return(unique(subset_sign))
  }
  
  all_sign <- rownames(signature_Colors[signature_Colors$ALLPhenotypeFinal=="Oligodendrocytes",])
  ind <- which(names(signature) %in% all_sign)
  subset_sign <- signature[ind]
  
  geneSet$olig_Myelinating <- union(subset_sign$Anna_PreMyelinatingOligo, subset_sign$Anna_MyelinatingOligo)
  geneSet$olig_Myelinating <- union(geneSet$olig_Myelinating, subset_sign$CNS_Myelinating.Oligodendrocytes) 
  
  geneSet$olig <- union(subset_sign$Anna_OligoLineage, subset_sign$CNS_Oligodendrocytes)
  geneSet$olig <- union(geneSet$olig, subset_sign$CNS_Newly.Formed.Oligodendrocyte)
  
  #geneSet$olig <- findSign("Oligodendrocytes")
  geneSet$Tcell <- findSign("Tcell")
  #geneSet$astro <- findSign("Astrocytes")
  geneSet$macro <-findSign("Macrophages")
  geneSet$micro <-findSign("Microglia")
  geneSet$neuro <-findSign("Neurons")
  
  GO2 <- gmt2GO("~/singleCell/GSEA/c2.cp.reactome.v7.4.symbols.gmt")
  reactome_cellcycle <- GO2$REACTOME_CELL_CYCLE
  
  usethis::use_data(geneSet, reactome_cellcycle, internal = TRUE, overwrite = TRUE)
  
}

