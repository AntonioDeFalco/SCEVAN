
#' annotate genes with reference to hg38.
#'
#' @param mat data matrix; genes in rows; cell names in columns.
#' @param ID.type gene id type: Symbol or Ensemble.
#' @param full.anno annotation file for all known genes, automatically loaded in copycat.
#'
#' @return annotations of each genes in rows with chrom and positions.
#'
#' @examples
#' test.anno.mat <- annotateGenes.hg20(mat=matx, ID.type="ENSEMBLE_id", full.anno = full.anno)
#' @export
annotateGenes.hg20 <- function(mat, ID.type="S"){
  
  load("~/singleCell/AllData/ClassTumorCells/copyKAT/sysdata.rda")
  
  print("start annotation ...")
  
  if(substring(ID.type,1,1) %in% c("E", "e")){
    shar <- intersect(rownames(mat), full.anno$ensembl_gene_id)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$ensembl_gene_id) %in% shar),]
    anno <- anno[!duplicated(anno$hgnc_symbol),]
    anno <- anno[order(match(anno$ensembl_gene_id, rownames(mat))),]
    data <- cbind(anno, mat)
    
  }else if(substring(ID.type,1,1) %in% c("S", "s")) {
    
    shar <- intersect(rownames(mat), full.anno$hgnc_symbol)
    mat <- mat[which(rownames(mat) %in% shar),]
    anno <- full.anno[which(as.vector(full.anno$hgnc_symbol) %in% shar),]
    anno <- anno[!duplicated(anno$hgnc_symbol),]
    anno <- anno[order(match(anno$hgnc_symbol, rownames(mat))),]
    data <- cbind(anno, mat)
  }
}



preprocessingMtx <- function(rawmat, id.type ="SYMBOL", ngene.chr=5, LOW.DR=0.05, UP.DR=0.1, n.cores=20, SMOOTH = TRUE){
  
  load("~/singleCell/AllData/ClassTumorCells/copyKAT/sysdata.rda")
  
  start_time <- Sys.time()
  set.seed(1)
  
  print(paste(" raw data - genes: ", nrow(rawmat), " cells: ", ncol(rawmat), sep=""))
  
  print("1) Filter: cells > 200 genes")
  
  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
  
  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells past filtering ", ncol(rawmat), " cells", sep=""))
  }
  
  print(paste0("2) Filter: genes > ", LOW.DR*100, "% of cells"))
  
  der <- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)
  
  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past filtering", sep=""))
  }
  
  if(nrow(rawmat) < 7000){
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  
  print("3) Annotations gene coordinates")
  
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE),]
  
  print(paste(nrow(anno.mat)," genes annotated", sep=""))
  
  print("4) Filter: genes involved in the cell cycle")
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
  if(length(toRev)>0){
    anno.mat <- anno.mat[-toRev, ]
  }
  
  print(paste(nrow(anno.mat)," genes past filtering ", sep=""))
  
  print(paste0("5)  Filter: cells > ", ngene.chr, "genes per chromosome "))
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<23|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }
  
  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")
  
  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }
  
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat) <-  colnames(rawmat3)
  
  print(paste("A total of ", ncol(norm.mat), " cells, ", nrow(norm.mat), " genes after preprocessing", sep=""))
  
  ##### smooth data ##### 
  if(SMOOTH){
    print("6) Smoothing data with dlm")
    dlm.sm <- function(c){
      model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
      x <- dlm::dlmSmooth(norm.mat[, c], model)$s
      x<- x[2:length(x)]
      x <- x-mean(x)
    }
    
    test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
    norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)
    rm(test.mc)
    colnames(norm.mat.smooth) <- colnames(norm.mat)
  }else{
    norm.mat.smooth <- norm.mat
    rm(norm.mat)
  }
  
  rownames(norm.mat.smooth) <- anno.mat$hgnc_symbol
  
  res <- list(norm.mat.smooth, rawmat3, anno.mat)
  names(res) <- c("norm.mat.smooth","rawmat3", "anno.mat")
  
  end_time<- Sys.time()
  print(end_time -start_time)
  
  return(res)
  
}


