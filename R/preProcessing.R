#' annotateGenes Annotate genes with genomic coordinates with reference to hg38 using Ensembl based annotation package
#'
#' @param mtx Count matrix with genes on row names (Ensemble or Symbol)
#'
#' @return Annotated matrix
#'
#' @examples
#' count_mtx_annot <- annotateGenes(count_mtx)
#' @export
annotateGenes <- function(mtx, organism = "human"){
  library(dplyr)

  if(organism == "human"){
    edb <- EnsDB_Hsapiens_v86 #From EnsDb.Hsapiens.v86
  }else{
    edb <- EnsDb_Mmusculus_v79 #From EnsDb.Hsapiens.v86
  }
  
  chr <- 1:22  
  edb <- (edb[edb$seqnames %in% chr,c(1,2,3,6,7)])
  edb$seqnames <- as.numeric(as.character(edb$seqnames))
  
    if (any(rownames(mtx) %in% edb$gene_name)){
      use_geneID <- "gene_name"
    }else{
      use_geneID <- "gene_id"
    }
  
  genes_inters <- intersect(rownames(mtx), edb[[use_geneID]])
  mtx <- mtx[which(rownames(mtx) %in% genes_inters),]
  edb <- edb[which(as.vector(edb[[use_geneID]]) %in% genes_inters),]
  edb <- edb[!duplicated(edb$gene_name),]
  edb <- edb[order(match(edb[[use_geneID]], rownames(mtx))),]
  if(class(mtx)[1]=="dgCMatrix"){
    mtx_annot <- cbind(edb, as.matrix(mtx))
  }else{
    mtx_annot <- cbind(edb, mtx)
  }
  
  return(mtx_annot)
}


#' preprocessingMtx  Pre-processing steps: Cells with less than 200 genes and the genes expressed in less than 1% of cells are removed. Genes are annotated and sorted 
#' according to genomic coordinates. Highly confident normal cells are sought in the matrix. Genes involved in the cell cycle pathway are removed.  Log-Freeman–Tukey transformation to stabilize variance 
#' and a polynomial dynamic linear modeling (DLM) to smooth out the outliers. 
#'
#' @param count_mtx raw count matrix
#' @param ngenes_chr minimum number of genes per chromosome (optional)
#' @param perc_genes percentage of cells in which each gene is to be expressed (optional)
#' @param par_cores number of cores (optional)
#' @param SMOOTH Boolean value to perform smoothing (optional)
#' @param findConfident Boolean value to search for normal cells (default TRUE)
#' @param AdditionalGeneSets List of additional signatures to be used to search for normal cells (optional)
#' @param SCEVANsignatures Boolean value TRUE to use internal SCEVAN signatures for normal cells or FALSE to use only signatures specified in AdditionalGeneSets (default TRUE)
#'
#' @return 
#' count_mtx_smooth processed and smoothed matrix
#' count_mtx_annot annotated matrix
#'
#' @examples
#' count_mtx_annot <- annotateGenes(count_mtx)
#' @export
preprocessingMtx <- function(count_mtx, sample, ngenes_chr=5, perc_genes=0.1, par_cores=20, findConfident = TRUE, AdditionalGeneSets = NULL, SCEVANsignatures = TRUE, organism = "human"){
  
  set.seed(123)
  
  print(paste(" raw data - genes: ", nrow(count_mtx), " cells: ", ncol(count_mtx), sep=""))
  
  print("1) Filter: cells > 200 genes")
  
  genes.raw <- apply(count_mtx, 2, function(x)(sum(x>0)))
  
  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  
  if(sum(genes.raw<100)>1){
    count_mtx <- count_mtx[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells past filtering ", ncol(count_mtx), " cells", sep=""))
  }
  
  der <- apply(count_mtx,1,function(x)(sum(x>0)))/ncol(count_mtx)
  
  if( sum(der > perc_genes) > 7000){
    print(paste0("2) Filter: genes > ", perc_genes*100, "% of cells"))
    count_mtx <- count_mtx[which(der > perc_genes), ]; 
  }else{
    perc_genes <- perc_genes - 0.05
    print("low data quality")
    print(paste0("2) Filter: genes > ", perc_genes*100, "% of cells"))
    count_mtx <- count_mtx[which(der > perc_genes), ]; 
  }
  print(paste(nrow(count_mtx)," genes past filtering", sep=""))

  #norm_cell <- getConfidentNormalCells(count_mtx, par_cores = par_cores)
  
  print("3) Annotations gene coordinates")
  
  count_mtx_annot <- annotateGenes(count_mtx, organism) 
  
  count_mtx <- count_mtx_annot[,-c(1:5)]
  rownames(count_mtx) <- count_mtx_annot$gene_name
  
  if(findConfident){
    norm_cell <- getConfidentNormalCells(count_mtx, sample, par_cores = par_cores, AdditionalGeneSets = AdditionalGeneSets, SCEVANsignatures = SCEVANsignatures, organism = organism)
  }else{
    norm_cell <- NULL
  }
  
  rm(count_mtx)
  
  count_mtx_annot <- count_mtx_annot[
    with(count_mtx_annot, order(as.numeric(as.character(seqnames)), as.numeric(as.character(end))), decreasing = FALSE),
  ]
  
  print(paste(nrow(count_mtx_annot)," genes annotated", sep=""))
  
  print("4) Filter: genes involved in the cell cycle")
  
  if(organism == "human"){
    totChr <- 22
    cellcycle <- reactome_cellcycle #From EnsDb.Hsapiens.v86
  }else{
    totChr <- 19
    cellcycle <- reactome_cellcycle_Mmusculus #From EnsDb.Hsapiens.v86
  }
  
  HLAs <- count_mtx_annot$gene_name[grep("^HLA-", count_mtx_annot$gene_name)]
  toRev <- which(count_mtx_annot$gene_name %in% c(as.vector(cellcycle), HLAs))
  if(length(toRev)>0){
    count_mtx_annot <- count_mtx_annot[-toRev, ]
  }
  
  print(paste(nrow(count_mtx_annot)," genes past filtering ", sep=""))
  
  print(paste0("5)  Filter: cells > ", ngenes_chr, "genes per chromosome "))
  cellsFilt <- NULL
  
  
  for(i in 6:ncol(count_mtx_annot)){
    cellChr <- cbind(count_mtx_annot$seqnames, count_mtx_annot[,i])
    cellChr <- cellChr[cellChr[,2]!=0,]
    if(length(rle(cellChr[,1])$length)<totChr|min(rle(cellChr[,1])$length)< ngenes_chr){
      cellsFilt <- c(cellsFilt, colnames(count_mtx_annot)[i])
    }
  }
  
  if(length(cellsFilt)==(ncol(count_mtx_annot)-5)) stop("all cells are filtered")

  if(length(cellsFilt)>0){
    count_mtx_annot <-count_mtx_annot[, -which(colnames(count_mtx_annot) %in% cellsFilt)]
  }
  
  if((ncol(count_mtx_annot)-5)<15) stop("Bad sample low cells < 15")
  
  print("6) Log Freeman Turkey transformation")
  
  count_mtx_proc <- data.matrix(count_mtx_annot[, 6:ncol(count_mtx_annot)])
  count_mtx_annot <- count_mtx_annot[, 1:5]
  count_mtx_norm <- log(sqrt(count_mtx_proc)+sqrt(count_mtx_proc+1))
  count_mtx_norm <- apply(count_mtx_norm,2,function(x)(x <- x-mean(x)))
  colnames(count_mtx_norm) <-  colnames(count_mtx_proc)
  rm(count_mtx_proc)
  
  print(paste("A total of ", ncol(count_mtx_norm), " cells, ", nrow(count_mtx_norm), " genes after preprocessing", sep=""))
  
  rownames(count_mtx_norm) <- count_mtx_annot$gene_name
  
  norm_cell <- norm_cell[names(norm_cell) %in% colnames(count_mtx_norm)]
  
  res <- list(count_mtx_norm, count_mtx_annot, norm_cell)
  names(res) <- c("count_mtx_norm", "count_mtx_annot", "norm_cell")
  
  return(res)
  
}


