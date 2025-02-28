

getCorrelationCN <- function(CNVref, CNVcomp){
  #TODO - add the standard output dir as var for the function
  
  df_1 <- as.data.frame(approx(CNVcomp$Pos,CNVcomp$Mean, seq(min(CNVref$Pos,CNVcomp$Pos),max(CNVref$Pos,CNVcomp$Pos), by = 1), ties = "ordered"))
  df_2 <- as.data.frame(approx(CNVref$Pos,CNVref$Mean, seq(min(CNVref$Pos,CNVcomp$Pos), max(CNVref$Pos,CNVcomp$Pos), by = 1), ties = "ordered"))
  
  testRes <- cor.test(df_1$y,df_2$y)
  testRes
}

# Chr Start End Mean
modifySEG <- function(segDF){
  test <- apply(segDF, 1, function(x) {
    start <- x[c(1,2,4)]
    names(start)[2] <- "Pos"
    end <- x[c(1,3,4)]
    names(end)[2] <- "Pos"
    as.data.frame(rbind(start,end))
  })
  segDF <- do.call(rbind,test)
  segDF
}

modifyPOS <- function(CNV, organism = "human"){
  
  if(organism == "human"){
    totChr <- 22
    add_chr <- sizeGRCh38
  }else{
    totChr <- 19
    add_chr <- sizeGRCm39
  }
  
  
  #add_chr <- read.table("/home/adefalco/singleCell/AllData/ClassTumorCells/VegaMC/sizeGRCh38.csv", header = TRUE)
  add_chr <- (add_chr$Size/1000)[1:(totChr-1)]
  
  CNV$Pos <- CNV$Pos/1000
  
  extr_chr <- unlist(lapply(1:totChr, function(x) max(which(CNV$Chr==x))))
  
  for (i in 1:(totChr-1)){
    CNV[(extr_chr[i]+1):extr_chr[(i+1)],]$Pos <- (CNV[(extr_chr[i]+1):extr_chr[(i+1)],]$Pos + sum(add_chr[1:i]))
  }
  
  CNV
}

addInterval <- function(segm, organism = "human", ref = 0){
  
  #add_chr <- read.table("/home/adefalco/singleCell/AllData/ClassTumorCells/VegaMC/sizeGRCh38.csv", header = TRUE)
  
  if(organism == "human"){
    totChr <- 22
    add_chr <- sizeGRCh38
  }else{
    totChr <- 19
    add_chr <- sizeGRCm39
  }
  
  
  toAdd <- c()
  for(x in 1:totChr){
    subset <- segm[segm$Chr==x,]
    if(nrow(subset)>0){
      if(min(subset$Pos)!=0) toAdd <- rbind(toAdd,data.frame(Chr = x, Pos = 0, End = min(subset$Pos)-1, Mean = ref))
      if(max(subset$End) < add_chr[add_chr$Chr==x,]$Size) toAdd <- rbind(toAdd,data.frame(Chr = x, Pos = max(subset$End)+1, End = add_chr[add_chr$Chr==x,]$Size, Mean = 0))
      if(nrow(subset)>1){
        for(i in 2:nrow(subset)){
          if(subset[i,]$Pos-1 - subset[i-1,]$End+1 > 2){
            toAdd <- rbind(toAdd,data.frame(Chr = x, Pos = subset[i-1,]$End+1, End = subset[i,]$Pos-1, Mean = ref))
          }
        }
      }
    }else{
      toAdd <- rbind(toAdd,data.frame(Chr = x, Pos = 0, End = add_chr[add_chr$Chr==x,]$Size, Mean = ref))
    }
  }
  
  segm <- rbind(segm,toAdd)
  segm <- segm[order(segm$Chr,segm$Pos),]
  segm
}

getModifyPosSeg <- function(x, organism = "human", ref = 0){
    mod <- addInterval(x, organism, ref)
    mod <- modifySEG(mod)
    mod <- modifyPOS(mod, organism)
    mod
}

#CNV c( "Chr","Start","End","Mean")
plotSegmentation <- function(CNV, organism = "human", modifyPosSeg = TRUE, CN = TRUE){
  
  if(ncol(CNV)>4){
    if(CN){
      CNV <- CNV[,-5]
      colnames(CNV)[4] <- "Mean"
      ref <- 2 
    }else{
      CNV <- CNV[,-4]
      colnames(CNV)[4] <- "Mean"
      ref <- 0
    }
  }else{
  
    indMean <- grepl("Alteration",colnames(CNV))
    if(any(indMean)){
      colnames(CNV)[indMean] <- "Mean"
      ref <- 2
    }else{
      ref <- 0
    }
  }
  
    #CNV$Call <- CNV$Call-2
    
  if(organism == "human"){
    totChr <- 22
  }else{
    totChr <- 19
  }
  
  if(modifyPosSeg) CNV <- getModifyPosSeg(CNV, organism = organism, ref = ref)
  
  par(cex=1, cex.main = 1.5, cex.lab = 1.5,xaxs="i")
  with(data = CNV,
       expr = {
         plot(x = Pos,
              y = Mean,
              pch = NA_integer_,ylab = "Copy number",  xlab="CHR", xaxt='n', type="l", xgap.axis = ref)
         polygon(x = c(min(Pos), Pos, max(Pos)),
                 y = c(ref, Mean, ref),
                 col = "red")
         clip(x1 = min(Pos),
              x2 = max(Pos),
              y1 = min(Mean),
              y2 = ref)
         polygon(x = c(min(Pos), Pos, max(Pos)),
                 y = c(ref, Mean, ref),
                 col = "blue")
         
       })
  
  abline(h = ref, col = "gray60", lwd = 1)
  
  extr_chr <- CNV[unlist(lapply(1:totChr, function(x) max(which(CNV$Chr==x)))),]$Pos
  
  extr_chr <- append(1, extr_chr)
  axis(1, extr_chr[2:(totChr+1)] - diff(extr_chr)/2, labels = 1:totChr, las = 1, line = 0.2, tick = 0,
       cex.axis = 1.0, gap.axis = 0)
  abline(v=extr_chr, col="black", lwd = 2)
}

getScevanCNVfinal <- function(sample , path = "", subclone = NULL){
  if(is.null(subclone)){
    #TODO
    # change from ./output to output dir
    CNV_infer <- read.table(paste0(path,"./output/",sample, "_Clonal_CN.seg"), sep="\t", header=TRUE, as.is=TRUE)
  }else{
    CNV_infer <- read.table(paste0(path,"./output/",sample,"_subclone",subclone,"_CN.seg"), sep="\t", header=TRUE, as.is=TRUE)
  }
  CNV_infer
}


getScevanCNV <- function(sample , path = "" , filter = FALSE, beta = ""){
  #TODO the read.table functions also refer to the ./output folder, meanign that 
  # those need to be changed around as well
  CNV_infer <- read.table(paste0(path,"./output/ ",sample, beta," vega_output"), sep="\t", header=TRUE, as.is=TRUE)
  
  if(filter) CNV_infer[!((CNV_infer$Mean<(-0.10) | CNV_infer$L.pv<0.01) | (CNV_infer$Mean>0.10 | CNV_infer$G.pv<0.01)),]$Mean <- 0 
  
  CNV_infer <- CNV_infer[,c(1,2,3,5)]
  colnames(CNV_infer)[1:4] <- c("Chr" , "Pos" , "End", "Mean")
  CNV_infer$Chr <- as.integer(gsub("chr","",CNV_infer$Chr))
  CNV_infer
}

heatmapConsensusPlot <- function(CNV,sample,file, CN = TRUE){
 
  if(ncol(CNV)>4){
    if(CN){
      CNV <- CNV[,-5]
      colnames(CNV)[4] <- "Mean"
      ref <- 2
      CNV$Mean <- CNV$Mean - ref
    }else{
      CNV <- CNV[,-4]
      colnames(CNV)[4] <- "Mean"
      CNV[abs(CNV$Mean)<0.05,]$Mean <- 0
      ref <- 0
    }
  }else{
    
    indMean <- grepl("Alteration",colnames(CNV))
    if(any(indMean)){
      colnames(CNV)[indMean] <- "Mean"
      ref <- 2
    }else{
      CNV[abs(CNV$Mean)<0.05,]$Mean <- 0
      ref <- 0
    }
  }
  
    x <- getModifyPosSeg(CNV, ref = ref)
  
    dfL <- as.data.frame(approx(x$Pos,x$Mean, seq(min(x$Pos),max(x$Pos), length.out = 1000000), ties = "ordered"))
    colnames(dfL) <- c("Pos","Mean")
    dfL$Mean[is.na(dfL$Mean)] <- 0
    
    dfL$Chr <- rep(1, nrow(dfL))
    for (ch in 1:21){
      maxCh <- max(x[x$Chr==ch,]$Pos)
      dfL[dfL$Pos>maxCh,]$Chr <- ch+1
    }
    
    chr <- as.numeric(dfL$Chr) %% 2+1
    rbPal1 <- colorRampPalette(c('black','grey'))
    CHR <- rbPal1(2)[as.numeric(chr)]
    chr1 <- cbind(CHR,CHR)
    
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    cells <- rbind(rbPal5(1),rbPal5(1))
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(n = 999)
    
    
    heatmap.3(t(cbind(dfL$Mean,dfL$Mean)),Rowv = FALSE, Colv = FALSE, dendrogram = "none", chr_lab = dfL$Chr, keysize=1, density.info="none", trace="none",
              cexRow=3.0,cexCol=2.0,cex.main=3.0,cex.lab=3.0,
              ColSideColors=chr1,
              symm=F,symkey=F,symbreaks=T,cex=3.0, main="", cex.main=4, margins=c(10,20), key=FALSE,
              notecol="black",col=my_palette, RowSideColors=t(cells))
    
    
}

plotCNclonal_old <- function(sample,ClonalCN, organism = "human"){
  
  if(ClonalCN) {
    fileNames <- c("ClonalCNProfile","onlytumor")
  }else{
    fileNames <- c("onlytumor")
  }
  
  for(name in fileNames){
    
    if(name == "ClonalCNProfile"){
      file = "_coarse-grained"
    }else{
      file = "_fine-grain_"
    }
    
    #TODO
    # change from ./output to output_dir
    segm <- getScevanCNV(paste0(sample,name))
    png(paste("./output/",sample,file,"ClonalCNProfile.png",sep=""), height=1050, width=2250, res=250) 
    plotSegmentation(CNV = segm, organism = organism)
    dev.off()
    png(paste("./output/",sample,file,"consensus.png",sep=""), height=650, width=3150, res=180)
    heatmapConsensusPlot(segm,sample,file)
    dev.off()
    }
}
  


plotCNclonal <- function(sample, ClonalCN, organism = "human"){
  
  if(ClonalCN) {
    fileNames <- c("ClonalCNProfile","onlytumor")
  }else{
    fileNames <- c("onlytumor")
  }
  
  for(name in fileNames){
    
    if(name == "ClonalCNProfile"){
      file = "_coarse-grained"
    }else{
      file = "_fine-grain_"
    }
    
    if(file == "_coarse-grained"){
      segm <- getScevanCNVfinal(sample)
    }else{
      segm <- getScevanCNV(paste0(sample,name))
    }
    
    #TODO
    # change from ./output to output_dir
    png(paste("./output/",sample,file,"ClonalCNProfile.png",sep=""), height=1050, width=2250, res=250) 
    plotSegmentation(CNV = segm, organism = organism)
    dev.off()
    png(paste("./output/",sample,file,"consensus.png",sep=""), height=650, width=3150, res=180)
    heatmapConsensusPlot(segm,sample,file)
    dev.off()
  }
}

