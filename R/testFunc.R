



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
  
  if(isTRUE(grep("\\.",(tum_cells)[1])==1) & isTRUE(grep("-",colnames(norm.mat.relat)[1])==1)){
    tum_cells <- gsub("\\.", "-",tum_cells)
  }else if( isTRUE(grep("-",(tum_cells)[1])==1) & isTRUE(grep("\\.",colnames(norm.mat.relat)[1])==1)){
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
    
    #print(paste("sCalinsky Max", sCalinsky[which.max(sCalinsky)]))
    
    hc.clus <- cutree(hcc,n_subclones)
    
    perc_cells_subclones <- table(hc.clus)/length(hc.clus)
    
    if(which.max(sCalinsky)>0.10 & all(perc_cells_subclones > 0.10)){
      
      print(paste("found", n_subclones, "subclones", sep = " "))
      names(perc_cells_subclones) <- paste0("percentage_cells_subsclone_",names(perc_cells_subclones))
      print(perc_cells_subclones)
      
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
      
      hcc <- hclust(parallelDist::parDist(t(results.com),threads = 20, method = "euclidean"), method = "ward.D")

      plotSubclones(CNAmat[,2], results.com,hcc, n_subclones, sample)
      
      hc.clus <- cutree(hcc,n_subclones)

    }else{
      n_subclones <- 0
      hc.clus <- 0
      print(paste("found", n_subclones, "subclones", sep = " "))
    }
  }
  
  
  ress <- list(n_subclones, breaks_subclones, results.com, hc.clus)
  
  names(ress) <- c("n_subclones", "breaks_subclones", "logCNA", "clustersSub")
  
  
  return(ress)
  
}







analyzeSegm <- function(sample, nSub = 1){
  
  all_segm <- list()
  
  for (i in 1:nSub){
    
    segm <- read.csv(paste0("./output/ ",sample,"_subclone",i," vega_output"), sep = "\t")
    segm <- segm[(segm$L.pv<0.001 | segm$G.pv<0.001) & (segm$Loss.Mean<(-0.30) | segm$Gain.Mean>0.30),c(1,2,3,6,7)]
    
    segm$Alteration <- "D"
    segm$Alteration[segm$G.pv<0.005] <- "A"
    segm <- segm[,c(1,2,3,6)]
    
    segm_new <- c()
    for (ch in unique(segm$Chr)) {
      segm_ch <- segm[segm$Chr==ch,]
      
      br <- 2
      
      while(nrow(segm_ch)>1){
        
        if(br>nrow(segm_ch)){
          break
        }
        
        if( (abs((segm_ch$End[(br-1)] - segm_ch$Start[br])) < 90000000) & (segm_ch$Alteration[(br-1)] == segm_ch$Alteration[br])){
          segm_ch$End[(br-1)] <- segm_ch$End[br]
          segm_ch <- segm_ch[-(br-1),]
        }else{
          br <- br + 1
        }
        
      }
      segm_new <- rbind(segm_new,segm_ch)
    }
    
    all_segm[[paste0(sample,"_subclone", i)]] <- segm_new
    
  }
  
  return(all_segm)
  
}

diffSubclones <- function(sampleAlter, sample, nSub = 2){
  
  all_segm_diff <- list()
  
  segm_sh <- c()
  
  for(sub in 1:nSub){
    
    if(sub==1){
      cl1 <- sampleAlter[[1]]
      cl2 <- sampleAlter[[2]]
    }else{
      cl2 <- sampleAlter[[1]]
      cl1 <- sampleAlter[[2]]
    }
    
    segm_new <- c()
    
    for (ch in sort(unique(union(cl1$Chr,cl2$Chr)))) {
      
      if(sum(cl1$Chr==ch)>0){
        
        cl1_ch <- cl1[cl1$Chr==ch,]
        cl2_ch <- cl2[cl2$Chr==ch,]
        
        for (br in 1:nrow(cl1_ch)) {
          
          FOUND <- FALSE
          AltPres <- which(cl2_ch$Alteration == cl1_ch[br,]$Alteration)
          
            if(length(AltPres)>0){
              
              for(br2 in AltPres){
                if( ((cl1_ch[br,]$Start >= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End <= cl2_ch[br2,]$End)) | ((cl1_ch[br,]$Start <= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End >= cl2_ch[br2,]$End))){
                  FOUND <- TRUE
                  
                  if(sub==1) segm_sh <- rbind(segm_sh,cl1_ch[br,]) 
                  
                  break
                }
              }
            }
            
            if(!FOUND){
              segm_new <- rbind(segm_new,cl1_ch[br,])  
            }
        }
      }
    }
    
    all_segm_diff[[paste0(sample,"_subclone", sub)]] <- segm_new
    
  }
  
  all_segm_diff[[paste0(sample,"clone")]] <- segm_sh
  
  return(all_segm_diff)
}



