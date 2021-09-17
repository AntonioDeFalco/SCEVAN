#' computeF1score compute F1 score
#'
#' @param pred inferred classification 
#' @param ground_truth  ground truth of classification
#'
#' @return F1 Score
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


#' subclonesTumorCells  Check the presence of subclonal structures in tumor cells, by determining the optimal number of clusters present in the CNA matrix of 
#' tumor cells using the Calinski-Harabasz Index as a criterion \cite{Calinski}. For each possible subclone the joint segmentation algorithm is 
#' applied and the separate segmentation results are analyzed to see if there are any significant alterations specific to one subclone compared to
#' the others.
#'
#' @param tum_cells malignant cells
#' @param CNAmat CNA matrix
#' @param samp sample name
#'
#' @return
#' n_subclones number of subclones
#' breaks_subclones breakpoints of subclones
#' logCNA CNA matrix
#' clustersSub clustering of subclones
#' @export
#'
#' @examples
subclonesTumorCells <- function(tum_cells, CNAmat, samp){
  
  norm.mat.relat <- CNAmat[,-c(1:3)]
  info_mat <- CNAmat[,c(1:3)]
  
  if(isTRUE(grep("\\.",(tum_cells)[1])==1) & isTRUE(grep("-",colnames(norm.mat.relat)[1])==1)){
    tum_cells <- gsub("\\.", "-",tum_cells)
  }else if( isTRUE(grep("-",(tum_cells)[1])==1) & isTRUE(grep("\\.",colnames(norm.mat.relat)[1])==1)){
    tum_cells <- gsub("-", "\\.",tum_cells)
  }
  
  norm.mat.relat <- norm.mat.relat[,tum_cells]
  
  n.cores = 20 
  distance="euclidean"
  dist_mtx <- parallelDist::parDist(t(norm.mat.relat),threads = n.cores, method = distance)
  
  hcc <- hclust(dist_mtx, method = "ward.D")
  hc.clus <- cutree(hcc, h = 15)
  
  n_subclones <- 0
  results.com <- NULL
  breaks_subclones <- NULL
  
  if(length(hc.clus) > 1){
    
    sCalinsky <- calinsky(hcc, dist_mtx, gMax = 10)
    n_subclones <- which.max(sCalinsky)
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
        
        breaks_subclones[[i]] <- getBreaksVegaMC(mtx_vega, CNAmat[,3], paste0(samp,"_subclone",i))
      }
      
      BR <- c()
      
      BR <- sort(unique(unlist(breaks_subclones)))
      
      logCNA <- computeCNAmtx(norm.mat.relat, BR, n.cores)
      
      results <- list(logCNA, BR)
      names(results) <- c("logCNA","breaks")
      
      colnames(results$logCNA) <- colnames(norm.mat.relat)
      results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
      
      hcc <- hclust(parallelDist::parDist(t(results.com),threads = 20, method = "euclidean"), method = "ward.D")

      plotSubclones(CNAmat[,2], results.com,hcc, n_subclones, samp)
      
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




analyzeSegm <- function(samp, nSub = 1){
  
  all_segm <- list()
  
  for (i in 1:nSub){
    
    segm <- read.csv(paste0("./output/ ",samp,"_subclone",i," vega_output"), sep = "\t")
    all_segm[[paste0(samp,"_subclone", i)]] <- getPossibleSpecAltFromSeg(segm)
    
  }
  
  return(all_segm)
  
}

getPossibleSpecAltFromSeg <- function(segm, name){
  
  segm <- segm[abs(segm$Mean)>0.10,]
  segm_new <- c()
  
  if(dim(segm)[1]>0){
    
    segm$Alteration <- "D"
    #segm$Alteration[segm$G.pv<5*10^{-5}] <- "A"
    segm$Alteration[segm$Mean>0.10] <- "A"
    segm <- segm[,c(1,2,3, ncol(segm))]
    
    for (ch in unique(segm$Chr)) {
      segm_ch <- segm[segm$Chr==ch,]
      
      br <- 2
      
      while(nrow(segm_ch)>1){
        
        if(br>nrow(segm_ch)){
          break
        }
        
        if( (abs((segm_ch$End[(br-1)] - segm_ch$Start[br])) < 20000000) & (segm_ch$Alteration[(br-1)] == segm_ch$Alteration[br])){
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
                if( ((cl1_ch[br,]$Start >= cl2_ch[br2,]$Start) & (cl1_ch[br,]$Start <= cl2_ch[br2,]$End)) | ((cl1_ch[br,]$End <= cl2_ch[br2,]$End) & (cl1_ch[br,]$End >= cl2_ch[br2,]$Start)) | ((cl2_ch[br2,]$Start >= cl1_ch[br,]$Start) & (cl2_ch[br2,]$Start <= cl1_ch[br,]$End)) | ((cl2_ch[br2,]$End <= cl1_ch[br,]$End) & (cl2_ch[br2,]$End >= cl1_ch[br,]$Start))){
                  #if( ((cl1_ch[br,]$Start >= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End <= cl2_ch[br2,]$End)) | ((cl1_ch[br,]$Start <= cl2_ch[br2,]$Start) | (cl1_ch[br,]$End >= cl2_ch[br2,]$End))){
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
    
    all_segm_diff[[paste0(samp,"_subclone", sub)]] <- segm_new
    
  }
  
  all_segm_diff[[paste0(samp,"_clone")]] <- segm_sh
  
  return(all_segm_diff)
}


testSpecificAlteration <- function(count_mtx, mtx_annot, listAltSubclones, clust_subclones, nSub = 2, samp){
  
  subclonesAlt <- list()

  segm_sh <- c()  
  
  nSubl <- grep("subclone",names(listAltSubclones))

  for(sub in nSubl){
    
    listAltSubclones[[sub]]$Mean <- 0

    segm_new <- c()  
    
    for (i in 1:nrow(listAltSubclones[[sub]])){

      subsetChr <- mtx_annot[mtx_annot$seqnames == listAltSubclones[[sub]][i,]$Chr,] 
      subsetCna <- count_mtx[mtx_annot$seqnames == listAltSubclones[[sub]][i,]$Chr,] 
      posSta <- which(subsetChr$end == listAltSubclones[[sub]][i,]$Start)
      posEnd <- which(subsetChr$end == listAltSubclones[[sub]][i,]$End)
      
      subClone1 <- subsetCna[posSta:posEnd, names(clust_subclones[clust_subclones==1])]
      subClone2 <- subsetCna[posSta:posEnd, names(clust_subclones[clust_subclones==2])]

      subClone1 <- apply(subClone1,1, mean)
      subClone2 <- apply(subClone2,1, mean)
      test <- t.test(subClone1,subClone2)
      
      print(listAltSubclones[[sub]][i,])
      print(test)
      #print(abs(mean(subClone1)-mean(subClone2)))
      #& abs(mean(subClone1)-mean(subClone2))>0.10
      if(test$p.value<10^{-15}){
        
        subCna <- as.matrix(subsetCna[posSta:posEnd,names(clust_subclones[clust_subclones==(nSub+1)-sub])])
        
        meanS <- mean(subCna)
        
        listAltSubclones[[sub]][i,]$Mean <- meanS
        
        segm_new <- rbind(segm_new, listAltSubclones[[sub]][i,])

      }else{
        segm_sh <- rbind(segm_sh, listAltSubclones[[sub]][i,])
     }
      
    }
    
    if(length(segm_new)>0)  subclonesAlt[[paste0(samp,"_subclone", sub)]] <- segm_new
  }

  listAltSubclones[[paste0(samp,"_clone")]]$Mean <- 0
  subclonesAlt[[paste0(samp,"_clone")]] <- rbind(segm_sh, listAltSubclones[[grep("_clone",names(listAltSubclones))]])

  subclonesAlt
  
  remov_alt <- c()
  for (i in 1:nrow(subclonesAlt[[paste0(samp,"_clone")]])){
    
    subsetChr <- mtx_annot[mtx_annot$seqnames == subclonesAlt[[paste0(samp,"_clone")]][i,]$Chr,] 
    subsetCna <- count_mtx[mtx_annot$seqnames == subclonesAlt[[paste0(samp,"_clone")]][i,]$Chr,] 
    posSta <- which(subsetChr$end == subclonesAlt[[paste0(samp,"_clone")]][i,]$Start)
    posEnd <-which(subsetChr$end == subclonesAlt[[paste0(samp,"_clone")]][i,]$End)
    
    subCna <- as.matrix(subsetCna[posSta:posEnd,names(clust_subclones)])
    
    meanS <- mean(subCna)
    
    subclonesAlt[[paste0(samp,"_clone")]][i,]$Mean <- meanS
    
      if(abs(meanS)<0.05){
        remov_alt <- append(remov_alt,i)
      }
  }
  if(length(remov_alt>0)) subclonesAlt[[paste0(samp,"_clone")]] <- subclonesAlt[[paste0(samp,"_clone")]][-remov_alt,]
  
  #subclonesAlt[[paste0(samp,"_clone")]] <- subclonesAlt[[paste0(samp,"_clone")]][order(as.numeric(as.character(subclonesAlt[[paste0(samp,"_clone")]]$Chr))),]

    for(sub in nSubl){
    subclonesAlt[[sub]] <- subclonesAlt[[sub]][order(abs(subclonesAlt[[sub]]$Mean), decreasing = TRUE),]
  }
  subclonesAlt[[paste0(samp,"_clone")]] <- subclonesAlt[[paste0(samp,"_clone")]][order(abs(subclonesAlt[[paste0(samp,"_clone")]]$Mean), decreasing = TRUE),]
  
  return(subclonesAlt)
}

annoteBand <- function(mtx_annot,diffSub){
  
  diffSub <- diffSub[unlist(lapply(diffSub, function(x) nrow(x)>0))]
  
  for(elem in 1:length(diffSub)){
    if(nrow(diffSub[[elem]])!=0){
      for(r in 1:nrow(diffSub[[elem]])){
        subset <- mtx_annot[mtx_annot$seqnames == diffSub[[elem]][r,]$Chr,]
        posSta <- which(subset$end == diffSub[[elem]][r,]$Start)
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
  
  if(length(grep("subclone",names(diffSub)))>0){
    
    vectAlt_all <- lapply(diffSub, function(x) apply(x, 1, function(x){ gsub(" ","",x); paste0(x[1],x[4])}))
    vectAlt_all <- lapply(vectAlt_all, function(x) { do.call(paste, c(as.list(unique(x)), sep = "\n"))})
    
    if(length(grep("subclone",names(diffSub)))>1){
      vectAlt1 <- vectAlt_all[[1]]
      vectAlt2 <- vectAlt_all[[2]]
      vectAltsh <- vectAlt_all[[3]]
    }else{
      vectAlt1 <- vectAlt_all[[1]]
      vectAlt2 <- c("")
      vectAltsh <- vectAlt_all[[2]]
    }
    
    vectAlt <- list(vectAlt1,vectAlt2,vectAltsh)
    
  }else{
    vectAlt <- list()
  }
  
  return(vectAlt)
  
}



genesDE <- function(count_mtx, count_mtx_annot, clustersSub, samp, specAlt, par_cores = 20){

  library(ggrepel)
  library(ggplot2)
  
  for(nsub in 1:length(specAlt)){
    for (i in 1:nrow(specAlt[[nsub]])){

      chrr <- specAlt[[nsub]][i,]$Chr
      startpos <- specAlt[[nsub]][i,]$Start
      endpos <- specAlt[[nsub]][i,]$End
      strr <- which(count_mtx_annot$seqnames==chrr & count_mtx_annot$end==startpos)
      endd <- which(count_mtx_annot$seqnames==chrr & count_mtx_annot$end==endpos)
      
      #top_genes <- rownames(count_mtx)
      top_genes <- count_mtx_annot$gene_name[strr:endd]
      
      cells_sub1 <- names(clustersSub[clustersSub==2])
      cells_sub2 <- names(clustersSub[clustersSub==1])
    
      parDE <- function(g){
        geneID <- top_genes[g]
        H1 = as.numeric(count_mtx[geneID, cells_sub1])
        H2 = as.numeric(count_mtx[geneID, cells_sub2])
        res1 = t.test(H1,H2)
        p_value = res1$p.value
        fc = mean(H1)-mean(H2)
        
        return(data.frame(p_value,fc,geneID))
      }
      test.mc <-parallel::mclapply(1:length(top_genes), parDE, mc.cores = par_cores)
      fact_spec2 <- do.call(rbind, test.mc)
      
      fact_spec2["p_value"][fact_spec2["p_value"]==0] <- (10^-1)*min(fact_spec2["p_value"][fact_spec2["p_value"]!=0])
      fact_spec2["p_value"] <- -log10(fact_spec2["p_value"])
      
      topGenes <- fact_spec2[
        with(fact_spec2, order(abs(fc), p_value, decreasing = c(TRUE,TRUE))),
      ][1:min(50,nrow(fact_spec2)),]
      
      print(topGenes)
      
      png(paste("./output/",samp,"- DE", "chr",chrr,":",startpos,":",endpos, "_subclones.png",sep=""), height=850, width=1250, res=150)

      txtRepel <- c()
      
      if(nrow(subset(topGenes, fc > 0))>0) {
        txtRepel <- append(txtRepel,geom_text_repel(data= subset(topGenes, fc > 0),aes(fc, p_value, label = geneID, colour = "blue", size = 30)))
      } 
      if(nrow(subset(topGenes, fc < 0))>0) {
          txtRepel <- append(txtRepel,geom_text_repel(data= subset(topGenes, fc < 0),aes(fc, p_value, label = geneID, colour = "red", size = 30)))
      } 
      
      p1 <- ggplot(fact_spec2, aes(fc, p_value, label = geneID)) + geom_point() + txtRepel +
        xlab("log2 Fold Change") + ylab("-log10 pvalue") + ggtitle(paste(samp,"- DE", "chr",chrr,":",startpos,":",endpos)) + theme_bw(base_size = 16) + 
        theme(legend.position="none") 
      
      plot(p1)
      
      dev.off()
  
    }
  }
  
}

pathwayAnalysis <- function(count_mtx, count_mtx_annot, clustersSub, samp, par_cores = 20){

  library(forcats)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  
  cells_sub1 <- names(clustersSub[clustersSub==2])
  cells_sub2 <- names(clustersSub[clustersSub==1])
  
  H1 = apply(count_mtx[,cells_sub1],1, mean)
  H2 = apply(count_mtx[,cells_sub2],1, mean)
  rankData = H1-H2
  
  names(rankData) <- rownames(count_mtx)

  #pathwaysH <- gmt2GO("/storage/qnap_home/adefalco/singleCell/GSEA/c2.cp.reactome.v7.4.symbols.gmt")

  fgseaRes <- fgseaMultilevel(pathwaysH,rankData , minSize=15, maxSize = 500, nproc = 1, nPermSimple = 10000, eps = 0)
  
  fgseaRes$pathway <- gsub("REACTOME_","",fgseaRes$pathway)
  
  fgseaRes <-  fgseaRes %>% dplyr::filter(padj < 0.1)
  
  topUp <- fgseaRes %>% 
    dplyr::filter(ES > 0) %>% 
    top_n(30, wt=-padj)
  topDown <- fgseaRes %>% 
    dplyr::filter(ES < 0) %>% 
    top_n(30, wt=-padj)
  topPathways <- bind_rows(topUp, topDown) %>% 
    arrange(-NES)
  

  png(paste("./output/",samp,"pathwayAnalysis_subclones.png",sep=""), width = 1600, height = 1080, units = "px", res=100)

  colnames(fgseaRes)[3] <- "pvalue"
  
  p1 <- ggplot(fgseaRes[fgseaRes$pathway %in% topPathways$pathway,], aes(x = NES, y = fct_reorder(pathway, NES))) + 
    geom_bar(stat='identity', aes(fill = pvalue)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="gray") +
    scale_fill_gradient(low="darkgreen", high = "gray")  +
    ylab(NULL) +
    ggtitle(paste(samp,"- Subclones Pathway Analysis"))
  plot(p1)
  
  dev.off()
  
}
