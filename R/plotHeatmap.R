# heatmap3 function modified from from "https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R"
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      hcr,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      chr_lab,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",
                      labels_gene=NULL,
                      ...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    ddc <- as.dendrogram(hcr)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    ddc <- as.dendrogram(hcr)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  
  totChr <- max(chr_lab)
  extr_chr <- unlist(lapply(1:totChr, function(x) max(which(chr_lab==x))))
  abline(v=extr_chr, col="black", lwd = 3)
  
  extr_chr <- append(1, extr_chr)
  axis(1, extr_chr[2:(totChr+1)] - diff(extr_chr)/2, labels = 1:totChr, las = 1, line = 0.2, tick = 0,
       cex.axis = cexCol, gap.axis = 0)
  
  if(!is.null(labels_gene)){
    for(i in 1:nrow(labels_gene)){
      axis(3, labels_gene[i,1], labels = labels_gene[i,2], las = 1, line = 0.2, tick = TRUE,
           cex.axis = 1.5, gap.axis = 0)
    }
  }
  
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}


plotCNA <- function(chr_lab, mtx_CNA, hcc, samp, pred = NULL, ground_truth = NULL){
  
  
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
  
  chr <- as.numeric(chr_lab) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  
  if (ncol(mtx_CNA)< 3000){
    h <- 10
  } else {
    h <- 15
  }
  
  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
  
  png(paste("./output/",samp,"heatmap.png",sep=""), height=h*250, width=4000, res=100)
  
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  prediction <- rep(rbPal5(1), ncol(mtx_CNA))
  cells <- rbind(prediction,prediction)
  
  if(length(pred)>0){

  prediction <- rbPal5(2)[as.numeric(factor(pred))]
  cells <- rbind(prediction,prediction)
  
  if(length(ground_truth) > 0){
    classCorr <- ground_truth[names(pred)]
    classCorr[classCorr!="malignant"] <- "non malignant"
    
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    Ground_Truth <- rbPal5(2)[as.numeric(factor(classCorr))]
    cells <- rbind(prediction,Ground_Truth)
  }
  
  }
  
  heatmap.3(t(mtx_CNA),dendrogram="r", hcr = hcc,
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE, chr_lab = chr_lab,
            keysize=1, density.info="none", trace="none",
            cexRow=3.0,cexCol=3.0,cex.main=3.0,cex.lab=3.0,
            symm=F,symkey=F,symbreaks=T,cex=3.0, main=paste("Heatmap ", samp), cex.main=4, margins=c(10,10))
  
  #legend("topright", paste("pred.",names(table(pred)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
  dev.off()
  
}



plotSubclones <- function(chr_lab, mtx_CNA, hcc, n_subclones, samp, par_cores=20){

chr <- as.numeric(chr_lab) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

hc.clus <- cutree(hcc,n_subclones)
rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:n_subclones])
subclones <- rbPal5(n_subclones)[as.numeric(factor(hc.clus))]
cells <- rbind(subclones,subclones)

if (ncol(mtx_CNA)< 3000){
  h <- 10
} else {
  h <- 15
}

png(paste("./output/",samp,"heatmap_subclones.png",sep=""), height=h*250, width=4000, res=100)

heatmap.3(t(mtx_CNA),dendrogram="r", hcr = hcc,
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE, chr_lab = chr_lab,
          keysize=1, density.info="none", trace="none",
          cexRow=3.0,cexCol=3.0,cex.main=3.0,cex.lab=3.0,
          symm=F,symkey=F,symbreaks=T,cex=3.0, main=paste("Heatmap ", samp), cex.main=4, margins=c(10,10))

dev.off()

}

plotTSNE <- function(raw_count_mtx, CNAmat , filt_genes, tum_cells, clustersSub, samp){
  
  #save(raw_count_mtx, CNAmat , filt_genes, tum_cells, clustersSub, file = paste0(samp,"plotTSNE.RData"))
  
  library(Rtsne)
  library(ggplot2)
  set.seed(1)
  newmtx <- raw_count_mtx[filt_genes,]
  
  pred <- rep("normal", length(colnames(newmtx)))
  names(pred) <- colnames(newmtx)
  pred[tum_cells] <- "tumor"
  
  if(class(newmtx)[1]=="dgCMatrix"){
    newmtx <- t(as.matrix(newmtx))
  }else{
    newmtx <- t(newmtx)
  }
  
  #tsne <- Rtsne(newmtx)
  
  tsne <- tryCatch(
    expr = {
      Rtsne(newmtx, check_duplicates = FALSE)
    },
    error = function(e){ 
      Rtsne(newmtx, check_duplicates = FALSE, perplexity = 15)
    }
  )
  
  df <- data.frame(x = tsne$Y[,1],
                   y = tsne$Y[,2],
                   CellType = pred)
  
  png(paste("./output/",samp,"tsne_scRNA.png",sep=""), height=1650, width=1650, res=200)
  
  pp <- ggplot(df, aes(x, y, colour = CellType)) +
    geom_point() + theme_bw() + scale_color_manual(breaks = c("tumor", "normal"),
                                                   values=c("red", "green"))
  plot(pp)
  dev.off()
  if(length(unique(clustersSub))>0){

    #tsne <- Rtsne(t(as.matrix(CNAmat[,tum_cells])))
    
    tsne <- tryCatch(
      expr = {
        #Rtsne(t(as.matrix(CNAmat[,tum_cells])))
        Rtsne(t(as.matrix(CNAmat[,names(clustersSub)])))
      },
      error = function(e){ 
        #Rtsne(t(as.matrix(CNAmat[,tum_cells])), perplexity = 15)
        Rtsne(t(as.matrix(CNAmat[,names(clustersSub)])), perplexity = 15)
      }
    )
    
    pred <- paste0("subclone_",clustersSub)
    names(pred) <- colnames(tum_cells)
    
    df <- data.frame(x = tsne$Y[,1],
                     y = tsne$Y[,2],
                     Subclones = pred)
    png(paste("./output/",samp,"tsne_CNA.png",sep=""), height=1650, width=1650, res=200)
    
    pp <- ggplot(df, aes(x, y, colour = Subclones)) +
      geom_point() + theme_bw() + scale_color_brewer(palette="Paired")
    
    plot(pp)
    dev.off()
  }
}

plotCNsubclones <- function(samp) {

  segmentation <- read.table(paste0("./output/ ",samp,"_subclone1 vega_output"), sep="\t", header=TRUE, as.is=TRUE)
  segmentation2 <- read.table(paste0("./output/ ",samp,"_subclone2 vega_output"), sep="\t", header=TRUE, as.is=TRUE)
  
  segmPos <- function(segmentation){
    segmentation$Start <- segmentation$Start/1000
    segmentation$End <- segmentation$End/1000
    extr_chr <- unlist(lapply(1:22, function(x) max(which(segmentation$Chr==x))))
    add_chr <- segmentation$End[extr_chr[1:(length(extr_chr)-1)]]
    sum_chr <- unlist(lapply(1:22, function(x) sum(segmentation$Chr==x)))
    for (i in 1:21){
      segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start + sum(add_chr[1:i]))
      segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End + sum(add_chr[1:i]))
    }
    
    x <- segmentation[1,]
    segm <- cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5]))
    for(i in 2:nrow(segmentation)){
      x <- segmentation[i,]
      abspos <- rbind(segm, cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5])))
      segm <- abspos
    }  
    colnames(segm) <- c("CHR","Pos","Mean")
    return(segm)
  }
  
  segm <- segmPos(segmentation)
  segm2 <- segmPos(segmentation2)
  
  df <- as.data.frame(approx(segm$Pos,segm$Mean, seq(min(segm$Pos), max(segm$Pos), length.out = 1000), ties = "ordered"))
  df2 <- as.data.frame(approx(segm2$Pos,segm2$Mean, seq(min(segm2$Pos), max(segm2$Pos), length.out = 1000), ties = "ordered"))
  colnames(df) <- c("Pos","Mean")
  colnames(df2) <- c("Pos","Mean")
  
  plot(df, ylab = "Copy number",  xlab="CHR", xaxt='n', type="l", xlim = c(min(segm[,2])+105000,max(segm[,2])-105000))#, ylim = c(min(segm[,3]-0.3),max(segm[,3])+0.3))
  
  plotSeg <- function(segm, col_lin){
    segmentationGain <- segm
    segmentationGain[segmentationGain$Mean<0,]$Mean <- 0
    lines(segmentationGain$Pos,segmentationGain$Mean, type="l", col=col_lin[1], lwd=3, pch=19)
    
    segmentationLoss <- segm
    segmentationLoss[segmentationLoss$Mean>0,]$Mean <- 0
    lines(segmentationLoss$Pos,segmentationLoss$Mean, type="l", col=col_lin[2], lwd=3, pch=19)
  }
  
  
  df_diff <- df
  thresh <- 0.08
  df_diff$Mean[abs(df$Mean - df2$Mean)<thresh] <- rowMeans(cbind(df$Mean[abs(df$Mean - df2$Mean)<thresh],df2$Mean[abs(df$Mean - df2$Mean)<thresh]))
  df_diff$Mean[abs(df$Mean - df2$Mean)>thresh] <- 0
  mmed <- function(x,n=5){runmed(x,n)} #Median
  df_diff$Mean <- mmed(df_diff$Mean, 11)
  plot(df_diff$Pos, df_diff$Mean, ylab = "Copy number",  xlab="CHR", xaxt='n', type="l", xlim = c(min(segm[,2])+105000,max(segm[,2])-105000))#, ylim = c(min(segm[,3]-0.3),max(segm[,3])+0.3))
  plotSeg(df_diff, c("red","blue"))
  
  df_Spec <- df
  df_Spec$Mean[abs(df$Mean - df2$Mean)<thresh] <- 0
  df_Spec$Mean <- mmed(df_Spec$Mean, 11)
  plotSeg(df_Spec, c("darkgreen","darkgreen"))
  
  
  df_Spec2 <- df2
  df_Spec2$Mean[abs(df$Mean - df2$Mean)<thresh] <- 0
  df_Spec2$Mean <- mmed(df_Spec2$Mean, 11)
  plotSeg(df_Spec2, c("purple","purple"))
  
  abline(h = 0, col = "gray60", lwd = 5)
  
  extr_chr <- segm[unlist(lapply(1:22, function(x) max(which(segm$CHR==x)))),]$Pos
  
  abline(v=extr_chr, col="black", lwd = 1)
  extr_chr <- append(1, extr_chr)
  axis(1, extr_chr[2:23] - diff(extr_chr)/2, labels = 1:22, las = 1, line = 0.2, tick = 0,
       cex.axis = 1.0, gap.axis = 0)
  
  legend("topright", inset=c(0,0), legend=c("Shared Gain", "Shared Loss", "Subclone 1", "Subclone 2"),
         col=c("red", "blue", "darkgreen", "purple"), lty=1, cex=0.7)

}

plotCN <- function(segmentation) {
  
  segmPos <- function(segmentation){
    segmentation$Start <- segmentation$Start/1000
    segmentation$End <- segmentation$End/1000
    extr_chr <- unlist(lapply(1:22, function(x) max(which(segmentation$Chr==x))))
    add_chr <- segmentation$End[extr_chr[1:(length(extr_chr)-1)]]
    sum_chr <- unlist(lapply(1:22, function(x) sum(segmentation$Chr==x)))
    for (i in 1:21){
      segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start + sum(add_chr[1:i]))
      segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End + sum(add_chr[1:i]))
    }
    
    x <- segmentation[1,]
    segm <- cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5]))
    for(i in 2:nrow(segmentation)){
      x <- segmentation[i,]
      abspos <- rbind(segm, cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5])))
      segm <- abspos
    }  
    colnames(segm) <- c("CHR","Pos","Mean")
    return(segm)
  }
  
  segm <- segmPos(segmentation)

  df <- as.data.frame(approx(segm$Pos,segm$Mean, seq(min(segm$Pos), max(segm$Pos), length.out = 1000), ties = "ordered"))
  colnames(df) <- c("Pos","Mean")

  plot(df, ylab = "Copy number",  xlab="CHR", xaxt='n', type="l", xlim = c(min(segm[,2])+105000,max(segm[,2])-105000))#, ylim = c(min(segm[,3]-0.3),max(segm[,3])+0.3))
  
  plotSeg <- function(segm, col_lin){
    segmentationGain <- segm
    segmentationGain[segmentationGain$Mean<0,]$Mean <- 0
    lines(segmentationGain$Pos,segmentationGain$Mean, type="l", col=col_lin[1], lwd=3, pch=19)
    
    segmentationLoss <- segm
    segmentationLoss[segmentationLoss$Mean>0,]$Mean <- 0
    lines(segmentationLoss$Pos,segmentationLoss$Mean, type="l", col=col_lin[2], lwd=3, pch=19)
  }
  

  
  plotSeg(df, c("red","blue"))
 
  abline(h = 0, col = "gray60", lwd = 5)
  
  extr_chr <- segm[unlist(lapply(1:22, function(x) max(which(segm$CHR==x)))),]$Pos
  
  abline(v=extr_chr, col="black", lwd = 1)
  extr_chr <- append(1, extr_chr)
  axis(1, extr_chr[2:23] - diff(extr_chr)/2, labels = 1:22, las = 1, line = 0.2, tick = 0,
       cex.axis = 1.0, gap.axis = 0)
  
}



segmPos <- function(segmentation){
  segmentation$Start <- segmentation$Start/1000
  segmentation$End <- segmentation$End/1000
  extr_chr <- unlist(lapply(1:22, function(x) max(which(segmentation$Chr==x))))
  
  add_chr <- sizeGRCh38
  add_chr <- (add_chr$Size/1000)[1:21]
  for (i in 1:21){
    segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start + sum(add_chr[1:i]))
    segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End + sum(add_chr[1:i]))
  }
  
  x <- segmentation[1,]
  segm <- cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5]))
  for(i in 2:nrow(segmentation)){
    x <- segmentation[i,]
    abspos <- rbind(segm, cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5])))
    segm <- abspos
  }  
  colnames(segm) <- c("CHR","Pos","Mean")
  return(segm)
}


segmPosSpec <- function(segmentation){
  
  segmentation <- segmentation[order(segmentation$Chr,segmentation$Start),]
  
  add_chr <- sizeGRCh38
  add_chr <- (add_chr$Size/1000)[1:21]
  
  extr_chr <- unlist(lapply(1:22, function(x) max(which(segmentation$Chr==x))))
  if(sum(extr_chr<0)>0){
    abs_chr <- which(extr_chr<0)
    for(i in 1:length(abs_chr)){
      segmentation <-  rbind(segmentation,c(as.character(abs_chr[i]),"10000" ,"100000","D",0))
    }
  }
  segmentation$Chr <- as.numeric(segmentation$Chr)
  segmentation$Start <- as.numeric(segmentation$Start)
  segmentation$End <- as.numeric(segmentation$End)
  
  segmentation <- segmentation[order(segmentation$Chr,segmentation$Start),]
  
  segmentation$Start <- segmentation$Start/1000
  segmentation$End <- segmentation$End/1000
  
  extr_chr <- unlist(lapply(1:22, function(x) max(which(segmentation$Chr==x))))
  
  for (i in 1:21){
    segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$Start + sum(add_chr[1:i]))
    segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End <- (segmentation[(extr_chr[i]+1):extr_chr[(i+1)],]$End + sum(add_chr[1:i]))
  }
  
  x <- segmentation[1,]
  segm <- cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5]))
  for(i in 2:nrow(segmentation)){
    x <- segmentation[i,]
    abspos <- rbind(segm, cbind(rbind(as.numeric(x[1]),as.numeric(x[1])),rbind(as.numeric(x[2]),as.numeric(x[3])), rbind(x[5],x[5])))
    segm <- abspos
  }  
  colnames(segm) <- c("CHR","Pos","Mean")
  return(segm)
}


plotCNAline <- function(segmList, segmListSpec, samp, nSub, colors_samp = NULL){
  
  segmListSpec <- lapply(segmListSpec, function(x) segmPosSpec(x))
  
  segmSUB_pos <- lapply(1:nSub, function(x) {
    subb <- paste0("subclone",x)
    if(sum(grepl(subb,names(segmListSpec)))>0) segmSUB1_pos <- segmListSpec[[grep(subb,names(segmListSpec))]]
  })
  
  for(i in 1:nSub){
    segmList[[i]][!(abs(segmList[[i]]$Mean)>0.10 | (segmList[[i]]$L.pv<=0.5 | segmList[[i]]$G.pv<=0.5)),]$Mean <- 0
  }

  segmList <- lapply(segmList, function(x) segmPos(x))
  
  minPos <- 1
  add_chr <- sizeGRCh38
  add_chr <- (add_chr$Size/1000)
  maxPos <- sum(add_chr)
  
  dfL <- lapply(segmList, function(x) {
    df <- as.data.frame(approx(x$Pos,x$Mean, seq(minPos,maxPos, length.out = 1000000), ties = "ordered"))
    colnames(df) <- c("Pos","Mean")
    df$Mean[is.na(df$Mean)] <- 0
    #df$Mean[abs(df$Mean)<0.10] <- 0
    return(df)
  })
  
  segm2_pos <- segmList[[grep("subclone1",names(dfL))]]
  df <- lapply(1:nSub, function(x) {
    subb <- paste0("subclone",x)
    dfL[[grep(subb,names(dfL))]]
  })
  
  
  df_VEGAchr <- as.data.frame(approx(segm2_pos$Pos,segm2_pos$CHR, seq(minPos,maxPos, length.out = 1000000), ties = "ordered"))
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 12, name = "RdBu")))(n = 999)
  
  chr <- as.numeric(df_VEGAchr$y) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  
  #res <- cor.test(df_1$Mean,df_2$Mean)
  
  if(length(colors_samp)==0) colors_samp <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:nSub])
  cells <- rbind(colors_samp(nSub),colors_samp(nSub))
  
  png(paste("./output/",samp,"consensus.png",sep=""), height=750, width=2850, res=180)
  heatmap.3(t(do.call(cbind,lapply(df, function(x) x$Mean))),Rowv = FALSE, Colv = FALSE, dendrogram = "none", chr_lab = 1:22, keysize=1, density.info="none", trace="none",
            cexRow=3.0,cexCol=2.0,cex.main=3.0,cex.lab=3.0,
            ColSideColors=chr1,
            symm=F,symkey=F,symbreaks=T,cex=3.0, main=paste("Sample:", samp), cex.main=4, margins=c(10,10), key=FALSE,
            notecol="black",col=my_palette, RowSideColors=cells)
  dev.off()
  
  df_share <- df[[1]]
  df_share$Mean <- apply(do.call(rbind,lapply(df, function(x) x$Mean)),2,mean)
  
  names(df) <- paste0("subclone",1:nSub)
  
  toREMdf_all <- lapply(1:nSub, function(x){
    
    toREMdf <- rep(FALSE,length(df[[x]]$Mean))
    
    if(sum(grepl(paste0("subclone",x),names(segmListSpec)))>0) {
      posOK <- segmSUB_pos[[x]]$Pos[abs(as.numeric(segmSUB_pos[[x]]$Mean))>0]
      
      for(i in 1:(length(posOK)/2)){
        #print(posOK[(i*2)-1])
        #print(posOK[(i*2)])
        toREMdf <- toREMdf | ((df[[x]]$Pos <= posOK[(i*2)] & df[[x]]$Pos >= posOK[(i*2)-1]))
      }
      sum(toREMdf)
    }
    
    toREMdf
    
  })
  
  for(i in 1:nSub){
    df[[i]]$Mean[!toREMdf_all[[i]]] <- 0
  }
  
  toREMdf_all <- Reduce("|", toREMdf_all)
  
  df_share$Mean[toREMdf_all] <- 0
  
  png(paste("./output/",samp,"plotCNline.png",sep=""), height=1050, width=2250, res=250)
  
  #df_share <- dfL[[grep("_clone",names(dfL))]]
  
  
  plot(df_share, ylab = "Copy number",  xlab="CHR", xaxt='n', type="l", xlim = c(min(segm2_pos[,2])+105000,max(segm2_pos[,2])-105000), main = samp, ylim = c(min(df_share$Mean,df[[1]]$Mean, df[[2]]$Mean)-0.05,max(df_share$Mean,df[[1]]$Mean, df[[2]]$Mean)+0.05))
  
  plotSeg <- function(segm, col_lin){
    segmentationGain <- segm
    if(sum(segmentationGain$Mean>0)>0){
      if(sum(segmentationGain$Mean<0)>0){
        segmentationGain[segmentationGain$Mean<0,]$Mean <- 0
      }
      lines(segmentationGain$Pos,segmentationGain$Mean, type="l", col=col_lin[1], lwd=3, pch=19)
      
    }
    
    segmentationLoss <- segm
    if(sum(segmentationLoss$Mean<0)>0){
      if(sum(segmentationLoss$Mean>0)>0){
        segmentationLoss[segmentationLoss$Mean>0,]$Mean <- 0
      }
      lines(segmentationLoss$Pos,segmentationLoss$Mean, type="l", col=col_lin[2], lwd=3, pch=19)
    }
    
  }
  
  plotSeg(df_share, c("red","blue"))
  
  if(length(colors_samp)==0) colors_samp <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:nSub])
  subclones <- colors_samp(nSub)
  
  lapply(1:nSub, function(x) plotSeg(df[[x]], c(subclones[x],subclones[x])))
  #lapply(df, function(x) plotSeg(x, c("purple","purple")))
  
  abline(h = 0, col = "gray60", lwd = 5)
  
  extr_chr <- segm2_pos[unlist(lapply(1:22, function(x) max(which(segm2_pos$CHR==x)))),]$Pos
  
  abline(v=extr_chr, col="black", lwd = 1)
  extr_chr <- append(1, extr_chr)
  axis(1, extr_chr[2:23] - diff(extr_chr)/2, labels = 1:22, las = 1, line = 0.2, tick = 0,
       cex.axis = 1.0, gap.axis = 0)
  
  #legend("topright", inset=c(0,0), legend=c("Shared Gain", "Shared Loss", "Subclone1", "Subclone2"),
  #       col=c("red", "blue", "darkgreen", "purple"), lty=1, cex=0.7)
  dev.off()
  
}



plotCNAlineOnlyTumor <- function(samp){
  
  segmList <- read.table(paste0("./output/ ",samp,"onlytumor vega_output"), sep="\t", header=TRUE, as.is=TRUE)
  
  segmList[!(abs(segmList$Mean)>0.10 | (segmList$L.pv<=0.5 | segmList$G.pv<=0.5)),]$Mean <- 0
  
  segmList <- segmPos(segmList)
  
  minPos <- 1
  add_chr <- sizeGRCh38
  add_chr <- (add_chr$Size/1000)
  maxPos <- sum(add_chr)
  
  x <- segmList
  dfL <- as.data.frame(approx(x$Pos,x$Mean, seq(minPos,maxPos, length.out = 1000000), ties = "ordered"))
  colnames(dfL) <- c("Pos","Mean")
  dfL$Mean[is.na(dfL$Mean)] <- 0
  
  segm2_pos <- segmList
  
  df_VEGAchr <- as.data.frame(approx(segm2_pos$Pos,segm2_pos$CHR, seq(minPos,maxPos, length.out = 1000000), ties = "ordered"))
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 12, name = "RdBu")))(n = 999)
  
  chr <- as.numeric(df_VEGAchr$y) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  cells <- rbind(rbPal5(1),rbPal5(1))
  
  png(paste("./output/",samp,"consensus.png",sep=""), height=750, width=2850, res=180)
  heatmap.3(t(cbind(dfL$Mean,dfL$Mean)),Rowv = FALSE, Colv = FALSE, dendrogram = "none", chr_lab = df_VEGAchr$y, keysize=1, density.info="none", trace="none",
            cexRow=3.0,cexCol=2.0,cex.main=3.0,cex.lab=3.0,
            ColSideColors=chr1,
            symm=F,symkey=F,symbreaks=T,cex=3.0, main=paste("Sample:", samp), cex.main=4, margins=c(10,10), key=FALSE,
            notecol="black",col=my_palette, RowSideColors=t(cells))
  dev.off()
  
  df_share <- dfL

  png(paste("./output/",samp,"plotCNline.png",sep=""), height=1050, width=2250, res=250)
  
  plot(df_share, ylab = "Copy number",  xlab="CHR", xaxt='n', type="l", xlim = c(min(segm2_pos[,2])+105000,max(segm2_pos[,2])-105000), main = samp)
  
  plotSeg <- function(segm, col_lin){
    segmentationGain <- segm
    if(sum(segmentationGain$Mean>0)>0){
      if(sum(segmentationGain$Mean<0)>0){
        segmentationGain[segmentationGain$Mean<0,]$Mean <- 0
      }
      lines(segmentationGain$Pos,segmentationGain$Mean, type="l", col=col_lin[1], lwd=3, pch=19)
      
    }
    
    segmentationLoss <- segm
    if(sum(segmentationLoss$Mean<0)>0){
      if(sum(segmentationLoss$Mean>0)>0){
        segmentationLoss[segmentationLoss$Mean>0,]$Mean <- 0
      }
      lines(segmentationLoss$Pos,segmentationLoss$Mean, type="l", col=col_lin[2], lwd=3, pch=19)
    }
    
  }
  
  plotSeg(df_share, c("red","blue"))
  
  abline(h = 0, col = "gray60", lwd = 5)
  
  extr_chr <- segm2_pos[unlist(lapply(1:22, function(x) max(which(segm2_pos$CHR==x)))),]$Pos
  
  abline(v=extr_chr, col="black", lwd = 1)
  extr_chr <- append(1, extr_chr)
  axis(1, extr_chr[2:23] - diff(extr_chr)/2, labels = 1:22, las = 1, line = 0.2, tick = 0,
       cex.axis = 1.0, gap.axis = 0)
  
  dev.off()
  
}

plotOncoHeat <- function(oncoHeat, nSub, samp, annotdf, mycolors, organism = "human"){
  
  legendBreaks <- sort(unique(as.numeric(as.matrix(oncoHeat))), decreasing = TRUE)
  legendLabels <-c("AMP","GAIN","","LOSS","DEL")
  i <- 1
  okLabels <- c()
  for(state in 2:-2){
    if(state %in% legendBreaks){
      okLabels <- c(okLabels,i)
    }
    i <- i+1
  }
  legendLabels <- legendLabels[okLabels]
  legendColors <- rev(c("blue","#add8e6","white","#FF7F7F","red"))
  legendColors <- rev(legendColors[okLabels])
  
  if(organism == "human"){
    listPos <- unlist(lapply(colnames(oncoHeat), function(x) strsplit(x, split = " ")[[1]][2]))
    
    dff <- data.frame(chr = as.numeric(substr(colnames(oncoHeat), 1,2)), let = substr(listPos,2,2), pos= as.numeric(substr(listPos,3,4)))
    
    dff[dff$let=="p",]$pos <- -dff[dff$let=="p",]$pos
    
    oncoHeat <- oncoHeat[,order(dff$chr, dff$let, dff$pos,decreasing = c(FALSE,FALSE,FALSE))]
  }else{
    listPos <- strsplit(colnames(oncoHeat),"-")
    ch <- as.numeric(unlist(lapply(listPos, function(x) x[1])))
    start <- as.numeric(unlist(lapply(listPos, function(x) x[2])))
    oncoHeat <- oncoHeat[,order(ch,start)]
  }
  
  png(paste("./output/",samp,"OncoHeat2.png",sep=""), height=2350, width=1050, res=150)
  pp <- pheatmap::pheatmap(t(oncoHeat), color = legendColors, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotdf, annotation_colors = mycolors, legend_breaks = legendBreaks, legend_labels = legendLabels,cellwidth = 30, annotation_legend = TRUE, fontsize = 14, labels_col = rep("",nrow(oncoHeat)))  
  h = grid::convertHeight(sum(pp$gtable$heights), "in", TRUE)
  w = grid::convertWidth(sum(pp$gtable$widths), "in", TRUE)
  ggplot2::ggsave(paste("./output/",samp,"OncoHeat.png",sep=""), device = "png", pp$gtable, width=w, height=h, dpi=300)
  dev.off()
}


plotOncoHeatSubclones <- function(oncoHeat, nSub, samp, perc_subclones, organism = "human"){
  
  if(!is.null(perc_subclones)){
    annotdf <- data.frame(row.names = rownames(oncoHeat), 
                          Subclone = rep(paste0("Subclone", seq(nSub), " (",round(perc_subclones*100,digits=2), "%)")) )  
  }else{
    annotdf <- data.frame(row.names = rownames(oncoHeat), 
                          Sample = rownames(oncoHeat) )  
  }
  
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:nSub])
  subclones <- rbPal5(nSub)
  names(subclones) <- unique(annotdf$Subclone)
  mycolors <- list(Subclone = subclones)
  
  plotOncoHeat(oncoHeat, nSub, samp, annotdf, mycolors, organism)
}


plotOncoHeat2 <- function(oncoHeat, nSub, samp, annotdf, mycolors){
  oncoHeat <- oncoHeat[,order(as.numeric(substr(colnames(oncoHeat), 1,2)),decreasing = FALSE)]
  
  png(paste("./output/",samp,"OncoHeat2.png",sep=""), height=2250, width=1450, res=150)
  pp <- pheatmap::pheatmap(t(oncoHeat), color = c("blue","white","red"), cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotdf, annotation_colors = mycolors, legend_breaks = c(1,0,-1), legend_labels = c("AMP","","DEL"),cellwidth = 30, annotation_legend = TRUE, fontsize = 14, labels_col = rep("",nrow(oncoHeat)))  
  h = grid::convertHeight(sum(pp$gtable$heights), "in", TRUE)
  w = grid::convertWidth(sum(pp$gtable$widths), "in", TRUE)
  ggplot2::ggsave(paste("./output/",samp,"OncoHeat.png",sep=""), device = "png", pp$gtable, width=w, height=h, dpi=300)
  dev.off()
}


modifyNameCNV <- function(CNV, CN = TRUE){
  
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
  
  return(list(CNV,ref))
}


plotCloneTree <- function(sample,res_subclones){

  library(tidytree)
  library(ape)
  library(ggtree)
  library(ggplot2)
  
  segmList <- lapply(1:res_subclones$n_subclones, function(x) read.table(paste0("./output/",sample,"_subclone",x,"_CN.seg"), sep="\t", header=TRUE, as.is=TRUE))
  #names(segmList) <- paste0("subclone",1:res_subclones$n_subclones)
  
  names(segmList) <- paste0(" ",1:res_subclones$n_subclones," ")
  
  segmList <- lapply(segmList, function(x) modifyNameCNV(x, CN = FALSE)[[1]])
  segmList <- lapply(segmList, function(x) getModifyPosSeg(x))
  
  #segmList <- lapply(segmList, function(x) segmPos(x))
  
  minPos <- 1
  add_chr <- sizeGRCh38
  add_chr <- (add_chr$Size/1000)
  maxPos <- sum(add_chr)
  
  dfL <- lapply(segmList, function(x) {
    df <- as.data.frame(approx(x$Pos,x$Mean, seq(minPos,maxPos, length.out = 1000000), ties = "ordered"))
    colnames(df) <- c("Pos","Mean")
    df$Mean[is.na(df$Mean)] <- 0
    #df$Mean[abs(df$Mean)<0.10] <- 0
    return(df)
  })
  
  dfLmean <- lapply(dfL, function(x) x$Mean)
  mtx <- do.call(rbind, dfLmean)
  
  distt <- parallelDist::parDist(mtx,threads = 10, method = "euclidean")

  PhyTree <- ape::nj(distt)
  
  colors_samp <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:res_subclones$n_subclones])
  colors <- colors_samp(res_subclones$n_subclones)
  
  branches <- lapply(1:res_subclones$n_subclones, function(x) x)
  names(branches) <- paste0("Subclone",1:res_subclones$n_subclones)
  tree <- groupOTU(PhyTree,branches)
  
  colors <- colors_samp(res_subclones$n_subclones)
  
  png(paste("./output/",sample,"CloneTree.png",sep=""), height=1650, width=1650, res=200)
  
  pp <- ggtree(tree, layout="daylight", size = 2) + 
    ggtitle(paste0(sample,"-Clone Tree")) + 
    geom_tiplab(aes(fill=group), geom = "label", size = 4) +
    scale_fill_manual(values=colors) + 
    theme_tree2(legend.position = "none") + 
    geom_nodepoint(color="#606060", alpha=1/3, size=10)  
  
  #plot(pp)
  
  plot(pp + xlim_expand(c(-5, 5), 'Dot'))
  
  dev.off()
}
