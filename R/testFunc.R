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


subclonesTumorCells <- function(tum_cells, CNAmat, relativeSmoothMtx, samp, n.cores, beta_vega, res_proc){
  library(scran)
  
  norm.mat.relat <- CNAmat[,-c(1:3)]
  info_mat <- CNAmat[,c(1:3)]
  
  if(isTRUE(grep("\\.",(tum_cells)[1])==1) & isTRUE(grep("-",colnames(norm.mat.relat)[1])==1)){
    tum_cells <- gsub("\\.", "-",tum_cells)
  }else if( isTRUE(grep("-",(tum_cells)[1])==1) & isTRUE(grep("\\.",colnames(norm.mat.relat)[1])==1)){
    tum_cells <- gsub("-", "\\.",tum_cells)
  }
  
  #save(tum_cells, norm.mat.relat, file = "debug.RData")
  norm.mat.relat <- norm.mat.relat[,tum_cells]
  
  library(igraph)
  dd <- 30
  if(dim(relativeSmoothMtx)[2]<=50){
    dd <- 5
  }
  
  #if(packageVersion("scran")>"1.16.0"){
  #  graph <- scran::buildSNNGraph(relativeSmoothMtx, k = 10, type = "number",  d =dd)#, type = "number")
  #}else{
  graph <- buildSNNGraph(relativeSmoothMtx, k = 10, type = "number", d =dd)
  #}install.packages("igraph")
  lc <- cluster_louvain(graph)
  #nSub <- 3
  #rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:nSub])
  #plot(graph, vertex.color=rbPal5(nSub)[lc$membership])
  
  hc.clus <- membership(lc)
  names(hc.clus) <- colnames(relativeSmoothMtx)
  
  perc_cells_subclones <- table(hc.clus)/length(hc.clus)
  
  removeSubcl <- as.numeric(names(perc_cells_subclones[perc_cells_subclones<0.02]))
  
  if(length(removeSubcl)>0) hc.clus <- hc.clus[!hc.clus %in% removeSubcl]
  
  n_subclones <- length(unique(hc.clus))
  
  results.com <- NULL
  breaks_subclones <- NULL
  
  if(n_subclones > 1 & modularity(lc)>0.13){
    
    perc_cells_subclones <- table(hc.clus)/length(hc.clus)
    
    print(paste("found", n_subclones, "subclones", sep = " "))
    names(perc_cells_subclones) <- paste0("percentage_cells_subsclone_",names(perc_cells_subclones))
    print(perc_cells_subclones)
    
    breaks_subclones <- list()
    
    logCNAl <- list()
    
    CNV <- list()
    
    colName <- c()
    
    for (i in 1:n_subclones){
      
      print(paste("Segmentation of subclone : ", i))
      
      tum_cells_sub1 <- names(hc.clus[hc.clus==i])
      colName <- append(colName,tum_cells_sub1)
      
      mtx_vega <- cbind(info_mat, norm.mat.relat[,tum_cells_sub1])
      
      colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
      
      breaks_subclones[[i]] <- getBreaksVegaMC(mtx_vega, CNAmat[,3], paste0(samp,"_subclone",i), beta_vega = 0.5)
      
      mtx_CNA3 <- computeCNAmtx(norm.mat.relat[,tum_cells_sub1], breaks_subclones[[i]], n.cores, rep(TRUE, length(breaks_subclones[[i]])))
      
      colnames(mtx_CNA3) <- tum_cells_sub1
      rownames(mtx_CNA3) <- rownames(norm.mat.relat)
      
      #save(mtx_CNA3, file = paste0("./output/",sample,"_mtx_CNA3.RData"))
      
      CNV[[i]] <- getCNcall(mtx_CNA3, res_proc$count_mtx_annot, sample = samp, subclone = i, par_cores = n.cores)
      
      mtx_CNA3 <- computeCNAmtx(norm.mat.relat[,tum_cells_sub1], breaks_subclones[[i]], n.cores, CNV[[i]]$Call != 2)
      
      logCNAl[[i]] <- mtx_CNA3
    }
    
    BR <- c()
    
    BR <- sort(unique(unlist(breaks_subclones)))
    
    paste0(samp,"_subclone",i)
    
    logCNA <- do.call(cbind, logCNAl)
    
    results <- list(logCNA, BR)
    names(results) <- c("logCNA","breaks")
    
    colnames(results$logCNA) <- colName #colnames(norm.mat.relat)
    results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
    
    hcc <- hclust(parallelDist::parDist(t(results.com),threads = 20, method = "euclidean"), method = "ward.D")
    
    plotSubclones(CNAmat[,2], results.com,hcc, n_subclones, samp)
    save(results.com, file = paste0("./output/",samp,"_CNAmtxSubclones.RData"))
    
  }else{
    n_subclones <- 0
    hc.clus <- 0
    print(paste("found", n_subclones, "subclones", sep = " "))
  }
  
  ress <- list(n_subclones, breaks_subclones, results.com, hc.clus)
  
  names(ress) <- c("n_subclones", "breaks_subclones", "logCNA", "clustersSub")
  
  return(ress)
}



ReSegmSubclones <- function(tum_cells, CNAmat, samp, hc.clus, n.cores, beta_vega){
  
    norm.mat.relat <- CNAmat[,-c(1:3)]
    info_mat <- CNAmat[,c(1:3)]
  
    if(isTRUE(grep("\\.",(tum_cells)[1])==1) & isTRUE(grep("-",colnames(norm.mat.relat)[1])==1)){
      tum_cells <- gsub("\\.", "-",tum_cells)
    }else if( isTRUE(grep("-",(tum_cells)[1])==1) & isTRUE(grep("\\.",colnames(norm.mat.relat)[1])==1)){
      tum_cells <- gsub("-", "\\.",tum_cells)
    }
    
    norm.mat.relat <- norm.mat.relat[,tum_cells]
    
    #n.cores = 20 
    #distance="euclidean"
    
  perc_cells_subclones <- table(hc.clus)/length(hc.clus)
  
  n_subclones <- length(unique(hc.clus))
  
  print(paste("found", n_subclones, "subclones", sep = " "))
  names(perc_cells_subclones) <- paste0("percentage_cells_subsclone_",names(perc_cells_subclones))
  print(perc_cells_subclones)
  
  breaks_subclones <- list()
  
  logCNAl <- list()
  
  for (i in 1:n_subclones){
    
    print(paste("Segmentation of subclone : ", i))
    
    tum_cells_sub1 <- names(hc.clus[hc.clus==i])
    
    mtx_vega <- cbind(info_mat, norm.mat.relat[,tum_cells_sub1])
    
    colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
    
    breaks_subclones[[i]] <- getBreaksVegaMC(mtx_vega, CNAmat[,3], paste0(samp,"_subclone",i), beta_vega = beta_vega)
    
    subSegm <- read.csv(paste0("./output/ ",paste0(samp,"_subclone",i)," vega_output"), sep = "\t")
    
    segmAlt <- abs(subSegm$Mean)>0.10 | (subSegm$G.pv<0.5 | subSegm$L.pv<0.5)
    #segmAlt <- rep(TRUE, length(subSegm$Mean))
    
    logCNAl[[i]] <- computeCNAmtx(norm.mat.relat[,tum_cells_sub1], breaks_subclones[[i]], n.cores, segmAlt)
  }
  
  BR <- c()
  
  BR <- sort(unique(unlist(breaks_subclones)))
  
  paste0(samp,"_subclone",i)
  
  logCNA <- do.call(cbind, logCNAl)
  
  results <- list(logCNA, BR)
  names(results) <- c("logCNA","breaks")
  
  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
  
  hcc <- hclust(parallelDist::parDist(t(results.com),threads = 20, method = "euclidean"), method = "ward.D")
  
  plotSubclones(CNAmat[,2], results.com, hcc, n_subclones, samp)
  save(results.com, file = paste0("./output/",samp,"_CNAmtxSubclones.RData"))
  
  ress <- list(n_subclones, breaks_subclones, results.com, hc.clus)
  
  names(ress) <- c("n_subclones", "breaks_subclones", "logCNA", "clustersSub")

  return(ress)
}

analyzeSegm <- function(samp, nSub = 1){
  
  all_segm <- list()
  
  if(nSub > 0){
    for (i in 1:nSub){
      
      segm <- read.csv(paste0("./output/ ",samp," _ ",i," _CN.seg"), sep = "\t")
      all_segm[[paste0(samp,"_subclone", i)]] <- getPossibleSpecAltFromSeg(segm)
      
    }
  }else{
    segm <- read.csv(paste0("./output/ ",samp," _  _CN.seg"), sep = "\t")
    all_segm <- getPossibleSpecAltFromSeg(segm)
  }
  
  return(all_segm)
  
}

analyzeSegm2 <- function(samp, nSub = 1){
  
  all_segm <- list()
  
  for (i in 1:nSub){
    
    segm <- read.csv(paste0("./output/ ",samp,"_subclone",i," vega_output"), sep = "\t")
    all_segm[[paste0(samp,"_subclone", i)]] <- getPossibleSpecAltFromSeg(segm)
    
  }
  
  return(all_segm)
  
}

getPossibleSpecAltFromSeg <- function(segm, name){
  
  colnames(segm)[4] <- "Mean"
  segm$Mean <- segm$Mean-2
  
  segm <- segm[segm$Mean!=0,]
  
  #segm <- segm[abs(segm$Mean)>=0.10 | (segm$G.pv<=0.5 | segm$L.pv<=0.5),]
  segm_new <- c()
  
  if(dim(segm)[1]>0){
    
    # segm$Alteration <- "D"
    #segm$Alteration[segm$G.pv<5*10^{-5}] <- "A"
    # segm$Alteration[segm$Mean>0.01] <- "A"
    #  
    #  segm <- segm[,c(1,2,3, ncol(segm))]
    
    colnames(segm)[2] <- "Start"
    colnames(segm)[4] <- "Alteration"
    
    for (ch in unique(segm$Chr)) {
      segm_ch <- segm[segm$Chr==ch,]
      
      br <- 2
      
      while(nrow(segm_ch)>1){
        
        if(br>nrow(segm_ch)){
          break
        }
        
        if( (abs((segm_ch$End[(br-1)] - segm_ch$Start[br])) < 10000000) & (abs((segm_ch$Alteration[(br-1)] - segm_ch$Alteration[br]))<2)){
          if(segm_ch$Alteration[(br-1)]<0){
            altt <- min(segm_ch$Alteration[(br-1)],segm_ch$Alteration[br])
          }else{
            altt <- max(segm_ch$Alteration[(br-1)],segm_ch$Alteration[br])
          }
          
          segm_ch$End[(br-1)] <- segm_ch$End[br]
          segm_ch <- segm_ch[-(br),]
          segm_ch$Alteration[(br-1)] <- altt
        }else{
          br <- br + 1
        }
        
      }
      segm_new <- rbind(segm_new,segm_ch)
    }
    
  }
  
  
  
  return(segm_new)
  
} 


getPossibleSpecAltFromSeg2 <- function(segm, name){
  
  segm <- segm[abs(segm$Mean)>=0.10 | (segm$G.pv<=0.5 | segm$L.pv<=0.5),]
  segm_new <- c()
  
  if(dim(segm)[1]>0){
    
    segm$Alteration <- "D"
    #segm$Alteration[segm$G.pv<5*10^{-5}] <- "A"
    segm$Alteration[segm$Mean>0.01] <- "A"
    segm <- segm[,c(1,2,3, ncol(segm))]
    
    for (ch in unique(segm$Chr)) {
      segm_ch <- segm[segm$Chr==ch,]
      
      br <- 2
      
      while(nrow(segm_ch)>1){
        
        if(br>nrow(segm_ch)){
          break
        }
        
        if( (abs((segm_ch$End[(br-1)] - segm_ch$Start[br])) < 10000000) & (segm_ch$Alteration[(br-1)] == segm_ch$Alteration[br])){
          segm_ch$End[(br-1)] <- segm_ch$End[br]
          segm_ch <- segm_ch[-(br),]
        }else{
          br <- br + 1
        }
        
      }
      segm_new <- rbind(segm_new,segm_ch)
    }
    
  }

  
  
  return(segm_new)
  
} 


diffSubclones <- function(sampleAlter, samp, nSub = 2){
  
  #save(sampleAlter, samp, nSub , file ="diffSubclones.RData")
  
  all_segm_diff <- list()
  
  vectSubcl <- 1:nSub
  
  for(sub in vectSubcl){
    segm_sh <- c()
    segm_new <- c()
    segm_sh_sub <- c()
    
    cl1 <- sampleAlter[[sub]]
    
    for (ch in 1:22) {
      
      if(sum(cl1$Chr==ch)>0){
        
        cl1_ch <- cl1[cl1$Chr==ch,]
        remInd <- c()
        for (br in 1:nrow(cl1_ch)) {
          
          
          sh_sub <- c()
          FOUND <- 0
          FOUND_small <- 0
          
          for(sub2 in vectSubcl[-sub]){
            
            cl2 <- sampleAlter[[sub2]]
            
            cl2_ch <- cl2[cl2$Chr==ch,]
            
            AltPres <- which(abs(cl2_ch$Alteration - cl1_ch[br,]$Alteration)<2)
            
            if(length(AltPres)>0){
              
              for(br2 in AltPres){
                FOUND_br2 <- FALSE
                
                if( abs(cl1_ch[br,]$Start - cl2_ch[br2,]$Start)<1000000 | abs(cl1_ch[br,]$End - cl2_ch[br2,]$End)<1000000){
                  
                  if( ((cl1_ch[br,]$Start >= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End <= cl2_ch[br2,]$End)) | ((cl1_ch[br,]$Start <= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End >= cl2_ch[br2,]$End))){
                    
                    cl_ch_new <- cl1_ch[br,]
                    
                    if(cl1_ch[br,]$Alteration<0){
                      cl_ch_new$Alteration <- min(cl1_ch[br,]$Alteration, cl2_ch[br2,]$Alteration)
                    }else{
                      cl_ch_new$Alteration <- max(cl1_ch[br,]$Alteration, cl2_ch[br2,]$Alteration)
                    }
                    
                    cl_ch_new$Start <- max(cl1_ch[br,]$Start, cl2_ch[br2,]$Start)
                    cl_ch_new$End <- min(cl1_ch[br,]$End, cl2_ch[br2,]$End)
                    
                    if( ((cl1_ch[br,]$End - cl1_ch[br,]$Start) / (cl2_ch[br2,]$End - cl2_ch[br2,]$Start) < 0.40 ) | ((cl2_ch[br2,]$End - cl2_ch[br2,]$Start) / (cl1_ch[br,]$End - cl1_ch[br,]$Start)  < 0.40 )){
                      
                      FOUND_small <- 1
                      
                    }else{
                      sh_sub <- append(sh_sub, sub2)
                      
                      FOUND <- FOUND + 1
                      FOUND_br2 <- TRUE
                      
                      remInd <- append(remInd, br2)
                      
                      break 
                    }
                  }
                  
                  
                }
              }
            }
            
            sampleAlter[[sub2]] <- sampleAlter[[sub2]][]
            
          }
          
          if(FOUND==(nSub-1)){
            segm_sh <- rbind(segm_sh,cl_ch_new) 
          }else if(FOUND>0){
            segm_sh_sub <- rbind(segm_sh_sub,cbind(cl_ch_new,sh_sub))
          }else if(FOUND==0 & FOUND_small<1){
            segm_new <- rbind(segm_new,cl1_ch[br,])  
          }
        }
        
      }
      
      all_segm_diff[[paste0(samp,"_subclone", sub)]] <- segm_new
      all_segm_diff[[paste0(samp,"_clone", sub)]] <- segm_sh
      all_segm_diff[[paste0(samp,"_share", sub)]] <- segm_sh_sub
    }
    
  }
  
  all_segm_diff
  
  return(all_segm_diff)
}  






diffSubclones2 <- function(sampleAlter, samp, nSub = 2){
  
  #save(sampleAlter, samp, nSub , file ="diffSubclones.RData")
  
  all_segm_diff <- list()
  
  vectSubcl <- 1:nSub
  
  for(sub in vectSubcl){
    segm_sh <- c()
    segm_new <- c()
    segm_sh_sub <- c()
    
    cl1 <- sampleAlter[[sub]]
    
    for (ch in 1:22) {
      
      if(sum(cl1$Chr==ch)>0){
        
        cl1_ch <- cl1[cl1$Chr==ch,]
        remInd <- c()
            for (br in 1:nrow(cl1_ch)) {
              
              
              sh_sub <- c()
              FOUND <- 0
              FOUND_small <- 0
              
              for(sub2 in vectSubcl[-sub]){
                
                cl2 <- sampleAlter[[sub2]]
                
                cl2_ch <- cl2[cl2$Chr==ch,]
           
                AltPres <- which(cl2_ch$Alteration == cl1_ch[br,]$Alteration)
          
               if(length(AltPres)>0){
              
                    for(br2 in AltPres){
                      FOUND_br2 <- FALSE
                      
                                if( abs(cl1_ch[br,]$Start - cl2_ch[br2,]$Start)<10000000 | abs(cl1_ch[br,]$End - cl2_ch[br2,]$End)<10000000){
                                    
                                  if( ((cl1_ch[br,]$Start >= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End <= cl2_ch[br2,]$End)) | ((cl1_ch[br,]$Start <= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End >= cl2_ch[br2,]$End))){
                                    
                                    cl_ch_new <- cl1_ch[br,]
                                    cl_ch_new$Start <- max(cl1_ch[br,]$Start, cl2_ch[br2,]$Start)
                                    cl_ch_new$End <- min(cl1_ch[br,]$End, cl2_ch[br2,]$End)
                                     
                                    if( ((cl1_ch[br,]$End - cl1_ch[br,]$Start) / (cl2_ch[br2,]$End - cl2_ch[br2,]$Start) < 0.40 ) | ((cl2_ch[br2,]$End - cl2_ch[br2,]$Start) / (cl1_ch[br,]$End - cl1_ch[br,]$Start)  < 0.40 )){
                                      
                                      FOUND_small <- 1
                     
                                    }else{
                                      sh_sub <- append(sh_sub, sub2)
                                      
                                      FOUND <- FOUND + 1
                                      FOUND_br2 <- TRUE
                                      
                                      remInd <- append(remInd, br2)

                                      break 
                                    }
                                    }
                                    
                                  
                                  }
                    }
              }
            
                sampleAlter[[sub2]] <- sampleAlter[[sub2]][]
                
              }
              
              if(FOUND==(nSub-1)){
                segm_sh <- rbind(segm_sh,cl_ch_new) 
              }else if(FOUND>0){
                segm_sh_sub <- rbind(segm_sh_sub,cbind(cl_ch_new,sh_sub))
              }else if(FOUND==0 & FOUND_small<1){
                segm_new <- rbind(segm_new,cl1_ch[br,])  
              }
    }
    
    }
    
    all_segm_diff[[paste0(samp,"_subclone", sub)]] <- segm_new
    all_segm_diff[[paste0(samp,"_clone", sub)]] <- segm_sh
    all_segm_diff[[paste0(samp,"_share", sub)]] <- segm_sh_sub
    }
    
  }
  
  all_segm_diff
  
  return(all_segm_diff)
}

testSpecificAlteration <- function(listAltSubclones, nSub = 2, samp){
  
  #save(count_mtx, mtx_annot, listAltSubclones, clust_subclones, nSub, samp, file = "testSpecificAlteration.RData")
  
  subclonesAlt <- list()

  segm_sh <- c()  
  
  nSubl <- grep("subclone",names(listAltSubclones))
  
  for(sub in nSubl){
    
    findInd <- regexpr(pattern ='subclone',names(listAltSubclones)[sub])
    subInd <- substr(names(listAltSubclones)[sub],findInd[1]+8,nchar(names(listAltSubclones)[sub]))
    
    subclonesAlt[[paste0(samp,"_subclone", subInd)]] <- listAltSubclones[[sub]] 
    
  }


  clone <- c()
  for(i in grep("_clone",names(listAltSubclones))){
    clone <- rbind(clone, listAltSubclones[[i]])
  }
  
  for (i in 1:nrow(clone)){
    
    duplShared <- (clone[,]$Chr == clone[i,]$Chr) & (clone[,]$Start == clone[i,]$Start) & (clone[,]$End == clone[i,]$End)
    duplShared[is.na(duplShared)] <- FALSE
    if(sum(duplShared)>1){
      clone <- clone[-which(duplShared)[-1],]
    }
  }
  
  for (i in 1:nrow(clone)){
    
    duplShared <- (clone[,]$Chr == clone[i,]$Chr) & (abs(clone[,]$Start - clone[i,]$Start)<10000000) & (abs(clone[,]$End - clone[i,]$End)<10000000)
    duplShared[is.na(duplShared)] <- FALSE
    
    if(sum(duplShared)>1){
      clone[which(duplShared)[1],]$Start <- min(clone[duplShared,]$Start)
      clone[which(duplShared)[1],]$End <-  max(clone[duplShared,]$End)
      clone <- clone[-which(duplShared)[-1],]
    }
  }
  
  clone <- clone[order(clone$Chr,clone$Start),]
  
  
  #subclonesAlt[[paste0(samp,"_clone")]] <- rbind(segm_sh)#, listAltSubclones[[grep("_clone",names(listAltSubclones))]])

  subclonesAlt[[paste0(samp,"_clone")]] <- clone
  
  #subclonesAlt[[paste0(samp,"_clone")]]$Mean <- 0
  
  if(nSub>2) subclonesAlt[[paste0(samp,"_shareSubclone")]] <- testSpecificSubclonesAlteration(listAltSubclones, nSub, samp)
  
  subclonesAlt <- lapply(subclonesAlt, function(x) x[x$Alteration!=0,])
  
  #subclonesAlt[[paste0(samp,"_clone")]] <- subclonesAlt[[paste0(samp,"_clone")]][order(as.numeric(as.character(subclonesAlt[[paste0(samp,"_clone")]]$Chr))),]

  nSubl <- grep("subclone",names(subclonesAlt))
  
  for(sub in nSubl){
    subclonesAlt[[sub]] <- subclonesAlt[[sub]][order(abs(subclonesAlt[[sub]]$Alteration), decreasing = TRUE),]
  }
  subclonesAlt[[paste0(samp,"_clone")]] <- subclonesAlt[[paste0(samp,"_clone")]][order(abs(subclonesAlt[[paste0(samp,"_clone")]]$Alteration), decreasing = TRUE),]
  
  return(subclonesAlt)
}


testSpecificSubclonesAlteration <- function(listAltSubclones, nSub = 2, samp){
  
  subclonesAltClone <- list()
  subclonesAlt <- list()
  
  segm_sh <- c()  
  
  nSubl <- grep("share",names(listAltSubclones))

  for(sub in nSubl){
    
    for (i in 1:nrow(listAltSubclones[[sub]])){
      
      duplShared <- (listAltSubclones[[sub]][,]$Chr == listAltSubclones[[sub]][i,]$Chr) & (listAltSubclones[[sub]][,]$Start == listAltSubclones[[sub]][i,]$Start) & (listAltSubclones[[sub]][,]$End == listAltSubclones[[sub]][i,]$End)
      otherSim <- listAltSubclones[[sub]][duplShared,]
      duplShared[is.na(duplShared)] <- FALSE
      if(sum(duplShared)>1){
        if(sum(duplicated(otherSim$sh_sub))>0){
          listAltSubclones[[sub]][which(duplShared)[1],]$sh_sub <- gsub(", ","-",toString(otherSim$sh_sub[-which(duplicated(otherSim$sh_sub))]))
        }else{
          listAltSubclones[[sub]][which(duplShared)[1],]$sh_sub <- gsub(", ","-",toString(otherSim$sh_sub))
        }
        
        listAltSubclones[[sub]] <- listAltSubclones[[sub]][-which(duplShared)[-1],]
      }
    }
    
  }
  
  for(sub in nSubl){
    
    listAltSubclones[[sub]]$Mean <- 0
    
    segm_new <- c() 
    
    for (i in 1:nrow(listAltSubclones[[sub]])){
      
      # subsetChr <- mtx_annot[mtx_annot$seqnames == listAltSubclones[[sub]][i,]$Chr,] 
      # subsetCna <- count_mtx[mtx_annot$seqnames == listAltSubclones[[sub]][i,]$Chr,] 
      # posSta <- which(subsetChr$end == listAltSubclones[[sub]][i,]$Start)
      # posEnd <- which(subsetChr$end == listAltSubclones[[sub]][i,]$End)
      # 
       findInd <- regexpr(pattern ='share',names(listAltSubclones)[sub])
       subInd <- substr(names(listAltSubclones)[sub],findInd[1]+5,nchar(names(listAltSubclones)[sub]))
      # 
      # #subInd <- substr(names(listAltSubclones)[sub],nchar(names(listAltSubclones)[sub]),nchar(names(listAltSubclones)[sub]))
      # subIndSh <- as.numeric(listAltSubclones[[sub]][i,]$sh_sub)
       subInd <- as.numeric(subInd)
      # #subIndSh <- as.numeric(otherSim$sh_sub)
      # 
      # subClone1 <- subsetCna[posSta:posEnd, names(clust_subclones[clust_subclones %in% c(subInd,subIndSh)])]
      # subClone2 <- subsetCna[posSta:posEnd, names(clust_subclones[!(clust_subclones %in% c(subInd,subIndSh))])]
      # 
      # subClone1 <- apply(subClone1, 1, mean)
      # subClone2 <- apply(subClone2, 1, mean)
      # 
      # if( listAltSubclones[[sub]][i,]$Alteration == "A" ){
      #   test <- t.test(subClone1,subClone2, alternative = "greater")
      # }else{
      #   test <- t.test(subClone1,subClone2, alternative = "less")
      # }
      # 
      # if(test$p.value<10^{-10} & abs(mean(subClone1)-mean(subClone2))>=0.05  & abs(mean(subClone1))>=0.10){
      #   
      #   listAltSubclones[[sub]][i,]$Mean <- mean(subClone1)
      #   
         listAltSubclones[[sub]][i,]$sh_sub <- paste(subInd, listAltSubclones[[sub]][i,]$sh_sub, sep = "-")
      #   
         segm_new <- rbind(segm_new, listAltSubclones[[sub]][i,])
      # }else{
      #   listAltSubclones[[sub]][i,]$Mean <- mean(append(subClone1,subClone2))
      #  segm_sh <- rbind(segm_sh, listAltSubclones[[sub]][i,])
      # }
      
    }
    
    if(length(segm_sh)>0)  subclonesAltClone[[paste0(samp,"_clone", subInd)]] <- segm_sh
    if(length(segm_new)>0)  subclonesAlt[[paste0(samp,"_share", subInd)]] <- segm_new

  }
  
  if(length(subclonesAlt)>0){
    
  subclonesAlt2 <- do.call(rbind,subclonesAlt)
  subclonesAlt2 <- subclonesAlt2[order(subclonesAlt2$Chr,subclonesAlt2$Start),]
  subclonesAlt2 <-subclonesAlt2[!duplicated(subclonesAlt2[,c(1,2,3,4)]),]
  
  for (i in 1:nrow(subclonesAlt2)){
    
    duplShared <- (subclonesAlt2[,]$Chr == subclonesAlt2[i,]$Chr) & (abs(subclonesAlt2[,]$Start - subclonesAlt2[i,]$Start)<10000000) & (abs(subclonesAlt2[,]$End - subclonesAlt2[i,]$End)<10000000) & (subclonesAlt2[,]$sh_sub == subclonesAlt2[i,]$sh_sub) & (abs(subclonesAlt2[,]$Alteration - subclonesAlt2[i,]$Alteration)<2) #(subclonesAlt2[,]$Alteration == subclonesAlt2[i,]$Alteration) 
    duplShared[is.na(duplShared)] <- FALSE
    
    if(sum(duplShared)>1){
      subclonesAlt2[which(duplShared)[1],]$Start <- min(subclonesAlt2[duplShared,]$Start)
      subclonesAlt2[which(duplShared)[1],]$End <-  max(subclonesAlt2[duplShared,]$End)
      subclonesAlt2 <- subclonesAlt2[-which(duplShared)[-1],]
    }
  }
  }else{
    subclonesAlt2 <- NULL
  }
  
  return(subclonesAlt2)
}

genesDE <- function(count_mtx, count_mtx_annot, clustersSub, samp, specAlt, par_cores = 20){

  #save(count_mtx, count_mtx_annot, clustersSub, samp, specAlt, par_cores, file = paste0(samp,"genesDE.RData"))
  
  library(ggrepel)
  library(ggplot2)
  
  for(nsub in 1:length(specAlt)){
    for (i in 1:nrow(specAlt[[nsub]])){

      chrr <- specAlt[[nsub]][i,]$Chr
      startpos <- specAlt[[nsub]][i,]$Start
      endpos <- specAlt[[nsub]][i,]$End
      
      strr <- which(count_mtx_annot$seqnames==chrr & count_mtx_annot$start==startpos)
      #strr <- which(count_mtx_annot$seqnames==chrr & count_mtx_annot$end==startpos)
      endd <- which(count_mtx_annot$seqnames==chrr & count_mtx_annot$end==endpos)
      
      #top_genes <- rownames(count_mtx)
      top_genes <- count_mtx_annot$gene_name[strr:endd]
      
      
      subInd <- substr(names(specAlt)[nsub],nchar(names(specAlt)[nsub]),nchar(names(specAlt)[nsub]))
      
      cells_sub1 <- names(clustersSub[clustersSub==subInd])
      cells_sub2 <- names(clustersSub[clustersSub!=subInd])
    
      parDE <- function(g){
        geneID <- top_genes[g]
        H1 = as.numeric(count_mtx[geneID, cells_sub1])
        H2 = as.numeric(count_mtx[geneID, cells_sub2])
        res1 = t.test(H1,H2)
        p_value = res1$p.value
        fc = mean(H1)-mean(H2)
        
        return(data.frame(p_value,fc,geneID))
      }
      
      #test.mc <-parallel::mclapply(1:length(top_genes), parDE, mc.cores = par_cores)
      
      if(Sys.info()["sysname"]=="Windows"){
        cl <- parallel::makeCluster(getOption("cl.cores", par_cores))
        test.mc <- parallel::parLapply(cl, 1:length(top_genes), parDE)
        parallel::stopCluster(cl)
      }else{
        test.mc <-parallel::mclapply(1:length(top_genes), parDE, mc.cores = par_cores)
      }
      
      fact_spec2 <- do.call(rbind, test.mc)
      
      fact_spec2["p_value"][fact_spec2["p_value"]==0] <- (10^-1)*min(fact_spec2["p_value"][fact_spec2["p_value"]!=0])
      fact_spec2["p_value"] <- -log10(fact_spec2["p_value"])
      
      topGenes <- fact_spec2[
        with(fact_spec2, order(abs(fc), p_value, decreasing = c(TRUE,TRUE))),
      ][1:min(50,nrow(fact_spec2)),]
      
      FOUND_SIGN_DE <- FALSE

      txtRepel <- c()
      
      if(nrow(subset(topGenes, fc > 0.5 & p_value>1.30103))>0) {
         txtRepel <- append(txtRepel,geom_text_repel(data= subset(topGenes, fc > 0.5 & p_value>1.30103), aes(fc, p_value, label = geneID, colour = "blue", size = 30)))
         FOUND_SIGN_DE <- TRUE
      } 
      if(nrow(subset(topGenes, fc < -0.5 & p_value>1.30103))>0) {
          txtRepel <- append(txtRepel,geom_text_repel(data= subset(topGenes, fc < -0.5 & p_value>1.30103), aes(fc, p_value, label = geneID, colour = "red", size = 30)))
          FOUND_SIGN_DE <- TRUE
      } 
      
      if(FOUND_SIGN_DE){
      
      png(paste("./output/",samp,"-DE", "chr",chrr,"-",startpos,"-",endpos, "_subclones.png",sep=""), height=850, width=1250, res=150)
      
      p1 <- ggplot(fact_spec2, aes(fc, p_value, label = geneID)) + geom_point() + txtRepel +
        xlab("log2 Fold Change") + ylab("-log10 pvalue") + ggtitle(paste(samp,"- DE", "chr",chrr,":",startpos,":",endpos)) + theme_bw(base_size = 16) + 
        theme(legend.position="none") 
      
      plot(p1)
      
      dev.off()
  
      }
    }
  }
  
}

pathwayAnalysis <- function(count_mtx, count_mtx_annot, clustersSub, samp, par_cores = 20){

  library(forcats)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  
  nSub <- length(unique(clustersSub))
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:nSub])
  subclones <- rbPal5(nSub)

  
  for(sub in unique(clustersSub)){
  cells_sub1 <- names(clustersSub[clustersSub==sub])
  cells_sub2 <- names(clustersSub[clustersSub!=sub])
  
  H1 = apply(count_mtx[,cells_sub1],1, mean)
  H2 = apply(count_mtx[,cells_sub2],1, mean)
  rankData = H1-H2
  
  names(rankData) <- rownames(count_mtx)

  #pathwaysH <- gmt2GO("/storage/qnap_home/adefalco/singleCell/GSEA/c2.cp.reactome.v7.4.symbols.gmt")

  fgseaRes <- fgseaMultilevel(pathwaysH,rankData , minSize=15, maxSize = 500, nproc = 1, nPermSimple = 10000, eps = 0)
  
  fgseaRes$pathway <- gsub("REACTOME_","",fgseaRes$pathway)
  
  fgseaRes <-  fgseaRes %>% dplyr::filter(padj < 0.05)
  
  topUp <- fgseaRes %>% 
    dplyr::filter(ES > 0) %>% 
    top_n(30, wt=-padj)
  #topDown <- fgseaRes %>% 
  #  dplyr::filter(ES < 0) %>% 
  #  top_n(30, wt=-padj)
  #topPathways <- bind_rows(topUp, topDown) %>% 
  topPathways <- topUp %>% 
    arrange(-NES)
  
  #save(fgseaRes,topPathways, file = paste(samp,"pathwayAnalysis_subclones",sub,".RDATA"))
  png(paste("./output/",samp,"pathwayAnalysis_subclones",sub,".png",sep=""), width = 1600, height = 1080, units = "px", res=100)

  colnames(fgseaRes)[3] <- "pvalue"
  
  p1 <- ggplot(fgseaRes[fgseaRes$pathway %in% topPathways$pathway,], aes(x = NES, y = fct_reorder(pathway, NES))) + 
    geom_bar(stat='identity', aes(fill = pvalue)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="gray") +
    scale_fill_gradient(low=subclones[sub], high = "gray")  +
    ylab(NULL) +
    ggtitle(paste(samp,"- Subclones", sub," Pathway Analysis"))
  plot(p1)
  
  dev.off()
  }
  
}



    
annoteBandOncoHeat <- function(mtx_annot,diffSub, nSub){
  
  diffSub <- diffSub[unlist(lapply(diffSub, function(x) nrow(x)>0))]
  
  for(elem in 1:length(diffSub)){
    if(nrow(diffSub[[elem]])!=0){
      for(r in 1:nrow(diffSub[[elem]])){
        subset <- mtx_annot[mtx_annot$seqnames == diffSub[[elem]][r,]$Chr,]
        posSta <- which(subset$start == diffSub[[elem]][r,]$Start)
        posEnd <- which(subset$end == diffSub[[elem]][r,]$End)
        geneToAnn <- subset[posSta:posEnd, ]$gene_name
        found_genes <- intersect(geneToAnn,biomartGeneInfo$geneSymbol)
        min_band <- biomartGeneInfo[which(biomartGeneInfo$geneSymbol %in% found_genes[1]),]$band 
        max_band <- biomartGeneInfo[which(biomartGeneInfo$geneSymbol %in% found_genes[length(found_genes)]),]$band 
        band_str <- paste(min_band,max_band, sep = "-")
        diffSub[[elem]][r,1] <- paste0(diffSub[[elem]][r,1], " (", band_str, ") ")
      }
    }
  }
  
  oncoHeat <- data.frame(row.names = paste0("Subclone",1:nSub))
  IndSub <- grep("subclone",names(diffSub))
  
  if(length(IndSub)>0){
    
    for(i in IndSub){
      
      for(j in 1:nrow(diffSub[[i]])){
        
        indNsub <- as.numeric(substr(names(diffSub)[i],nchar(names(diffSub)[i]),nchar(names(diffSub)[i])))
         vect <- rep(0, nSub)
        # if(diffSub[[i]][j,]$Alteration=="A"){
        #   vect[indNsub] <- 1
        # }else{
        #   vect[indNsub] <- -1
        # }
        
        vect[indNsub] <- diffSub[[i]][j,]$Alteration
        
        oncoHeat[diffSub[[i]][j,]$Chr]  <- vect
      }
      
    }
  }
  
  
  IndSub <- grep("_clone",names(diffSub))
  
  if(length(IndSub)>0){
    
    for(i in IndSub){
      
      for(j in 1:nrow(diffSub[[i]])){
        
        #if(abs(diffSub[[i]][j,]$Mean)>0.10){
          
          indNsub <- as.numeric(substr(names(diffSub)[i],nchar(names(diffSub)[i]),nchar(names(diffSub)[i])))
          
          # if(diffSub[[i]][j,]$Alteration=="A"){
          #   vect <- rep(1, nSub)
          # }else{
          #   vect <- rep(-1, nSub)
          # }
          # 
          vect <- rep(diffSub[[i]][j,]$Alteration, nSub)
          
          oncoHeat[diffSub[[i]][j,]$Chr]  <- vect
        #}
      }
      
    }
  }
  
  IndSub <- grep("shareSubclone",names(diffSub))
  
  if(length(IndSub)>0){
    
    
    for(i in IndSub){
      
      for(j in 1:nrow(diffSub[[i]])){
        
        indNsub <- as.numeric(substr(names(diffSub)[i],nchar(names(diffSub)[i]),nchar(names(diffSub)[i])))
        
        indNsub <- as.integer(unlist(strsplit(diffSub[[i]][j,]$sh_sub,"-")))
        
        vect <- rep(0, nSub)
        # if(diffSub[[i]][j,]$Alteration=="A"){
        #   vect[indNsub] <- 1
        # }else{
        #   vect[indNsub] <- -1
        # }
        vect[indNsub] <- diffSub[[i]][j,]$Alteration
        
        oncoHeat[diffSub[[i]][j,]$Chr]  <- vect
      }
      
    }
  }
  
  return(oncoHeat)
  
}




