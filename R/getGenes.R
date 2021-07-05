getGenes <-
function(out_file_name, segmentation, pval_threshold, baf, ensembl_dataset,
        mart_database, html, correction){
## This function generates the html file reporting the list
## (with the respective details) of each gene overlapping a significant region

    eps <- 0.00000001;
    message("Connecting to BioMart \"", mart_database, "\" Database")
    if(mart_database=="ensembl"){
        ensembl <- useMart(mart_database, dataset=ensembl_dataset)
    }else{
        ensembl <- useMart(mart_database, archive=TRUE,
                dataset=ensembl_dataset)
    }
    ensembl_url <- "http://www.ensembl.org/Gene/Summary?g="

    del_genes <- matrix(nrow=0, ncol=12)
    amp_genes <- matrix(nrow=0, ncol=12)
    if(baf==TRUE){
        loh_genes <- matrix(nrow=0, ncol=12)
    }

    message("Downloading Gene Information from \"", ensembl_dataset,
        "\" Database")
    chromosomes <- sort(unique(segmentation[,1]))
    for(chr in chromosomes){
        message(".", appendLF = FALSE)
        curr_seg <- segmentation[which(as.character(segmentation[,1])
                        ==as.character(chr)),]
        if(class(curr_seg)!="matrix"){
            curr_seg <- t(as.matrix(curr_seg))
        }
        for(j in 1:nrow(curr_seg)){
            curr_reg <- curr_seg[j,]
            if(min(c(as.numeric(curr_reg[6]),as.numeric(curr_reg[7]),
            as.numeric(curr_reg[8])))<=pval_threshold){

                genes_details <- getBM(attributes=c("ensembl_gene_id",
                        "external_gene_id", "chromosome_name",
                        "start_position", "end_position",
                        "band",    "strand", "description"),
                        filters=c("chromosome_name", "start", "end"),
                        values=list(curr_reg[1], curr_reg[2], curr_reg[3]),
                        mart=ensembl)
               
                if(nrow(genes_details)>0){
                    if(as.numeric(curr_reg[6])<=pval_threshold){
                        new <- cbind(genes_details,
			    round(as.numeric(curr_reg[6]),5), curr_reg[9], 
			    as.numeric(curr_reg[13]),
			    round(as.numeric(curr_reg[16]),5)
			    )
                        del_genes <- rbind(del_genes, new)
                    }
                    if(as.numeric(curr_reg[7])<=pval_threshold){
                        new <- cbind(genes_details,
			    round(as.numeric(curr_reg[7]),5), curr_reg[10], 
			    as.numeric(curr_reg[14]),
			    round(as.numeric(curr_reg[17]),5)
			    )
                        amp_genes <- rbind(amp_genes, new)
                    }
                    if(baf==TRUE && as.numeric(curr_reg[8])<=pval_threshold){
                        new <- cbind(genes_details,
			    round(as.numeric(curr_reg[8]),5), curr_reg[11], 
			    as.numeric(curr_reg[15]),
			    round(as.numeric(curr_reg[18]),5)
			    )
                        loh_genes <- rbind(loh_genes, new)
                    }
                }
            }
        }
    }

    colnames(del_genes) <- c("Ensembl Gene ID", "External Gene ID",
			     "Chromosome", "Gene Start", "Gene End", "Cytoband",
			     "Strand",  "Description", "p-val", "Frequency",
			     "Mean", "Focal Score")
    delfile <- paste(substr(out_file_name, 1, nchar(out_file_name)-5),
		    "DelGenes", sep="")
   
    ## Matrix containing the set of genes overlapping lost regions
    del_genes[,1] <- paste("<a href=\"", ensembl_url, del_genes[,1], "\">",
                del_genes[,1] ,"</a>", sep="")
    del_genes[,4:5] <- format(del_genes[,4:5], big.mark=".")
    

    colnames(amp_genes) <- c("Ensembl Gene ID", "External Gene ID",
			     "Chromosome", "Gene Start", "Gene End", "Cytoband",
			     "Strand",  "Description", "p-val", "Frequency",
			     "Mean", "Focal Score")
    ampfile <- paste(substr(out_file_name, 1, nchar(out_file_name)-5),
		    "AmpGenes", sep="")
   
    ## Matrix containing the set of genes overlapping amplified regions
    amp_genes[,1] <- paste("<a href=\"", ensembl_url, amp_genes[,1], "\">",
                amp_genes[,1] ,"</a>", sep="")
    amp_genes[,4:5] <- format(amp_genes[,4:5], big.mark=".")
  

    ## Matrix containing the set of genes overlapping LOH regions
    if(baf==TRUE){
	colnames(loh_genes) <- c("Ensembl Gene ID", "External Gene ID",
			     "Chromosome", "Gene Start", "Gene End", "Cytoband",
			     "Strand",  "Description", "p-val", "Frequency",
			     "Mean", "Focal Score")
	lohfile <- paste(substr(out_file_name, 1, nchar(out_file_name)-5),
		    "LOHGenes", sep="")
	
        loh_genes[,1] <- paste("<a href=\"", ensembl_url, loh_genes[,1], "\">",
                loh_genes[,1] ,"</a>", sep="")
        loh_genes[,4:5] <- format(loh_genes[,4:5], big.mark=".")
       
    }

    ## Create the Header
    inFile <- system.file("template/HeaderGenes.html", package="VegaMC")
    inFile <- file(inFile, "r")
    outFile <- file(out_file_name, "w")
    index=""
    if(html==TRUE){
        html_genes <- paste(substr(out_file_name, 1, nchar(out_file_name)-10),
                ".html", sep="")
    }
    values <- c(TITLE=out_file_name, DATE=date(), REGFILE=html_genes)
    values <- c(values, NDEL=nrow(del_genes))
    values <- c(values, NAMP=nrow(amp_genes))
    if(baf==TRUE){
        values <- c(values, NLOH=nrow(loh_genes))
    }else{
        values <- c(values, NLOH="NA")
    }
    copySubstitute(inFile, outFile, values)
    close(inFile)

   
    if(nrow(del_genes)>0){
	del_genes <- as.matrix(del_genes)
	uniq_gene <- unique(del_genes[,1])
	## There are duplicate genes
	if(length(uniq_gene) != nrow(del_genes)){
	    for(i in 1:length(uniq_gene)){
		pos <- which(del_genes[,1]==uniq_gene[i])
		if(length(pos)>1){
		    pv <- round(min(as.numeric(del_genes[pos,9]))+eps, 5)
		    del_genes[pos[1],9] <- pv
				
		    fr <- as.character(del_genes[pos,10])
		    fr <- substr(fr, 1, nchar(fr)-1)
		    del_genes[pos[1],10] <- 
			paste(round(max(as.numeric(fr)), 1), "%", sep="")

		    m <- round(mean(as.numeric(del_genes[pos,11]))+eps, 5)
		    del_genes[pos[1],11] <- m
				
		    fs <- round(mean(as.numeric(del_genes[pos,12]))+eps, 5)
		    del_genes[pos[1],12] <- fs
		
		    del_genes <- del_genes[-pos[2:length(pos)],]
		}
		
	    }
	}

	del_genes <-
	    del_genes[order(as.numeric(substr(del_genes[,10], 1,
	    nchar(del_genes[,10])-1)), decreasing=TRUE),]
	
	 write.table(del_genes, delfile, sep="\t", col.names=TRUE,
		row.names=FALSE, eol="\n", quote=FALSE) 

	del_genes[,c(9,11,12)] <-
	format(round(as.numeric(del_genes[,c(9,11,12)]), 5), decimal.mark=",")
	
	
	## Create the Header for the Table of Deleted Genes
        inFile <- system.file("template/HeaderTableAberrantGenes.html",
                    package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, TABNAME="Deleted", ABERR="Loss")
	if(correction==TRUE){
	    values <- c(values, SIGN="q")
	}else{
	    values <- c(values, SIGN="p")
	}
        copySubstitute(inFile, outFile, values)
        close(inFile)
        ## Insert a new Line for each Gene
        for(i in 1:nrow(del_genes)){
	    inFile <- system.file("template/AberrantGeneTableRow.html",
		    package="VegaMC")
	    inFile <- file(inFile, "r")
	    values <- c(ENSID=as.character(del_genes[i,1]))
	    values <- c(values, GENEID=as.character(del_genes[i,2]))
	    values <- c(values, CHR=as.character(del_genes[i,3]))
	    values <- c(values, BPSTART=as.character(del_genes[i,4]))   
	    values <- c(values, BPEND=as.character(del_genes[i,5]))
	    if(is.na(del_genes[i,6])){
		values <- c(values, CYTO="")
	    }else{
		values <- c(values, CYTO=as.character(del_genes[i,6]))
	    }
	    if(is.na(del_genes[i,7])){
		values <- c(values, STRAND="")
	    }else{
		values <- c(values, STRAND=as.character(del_genes[i,7]))
	    }
	    if(is.na(del_genes[i,8])){
		values <- c(values, DESCR="")
	    }else{
		values <- c(values, DESCR=as.character(del_genes[i,8]))
	    }
	    values <- c(values, ABERRPVAL=as.character(del_genes[i,9]))
	    values <- c(values, ABERRPRC=as.character(del_genes[i,10]))
	    values <- c(values, MEAN=as.character(del_genes[i,11]))
	    values <- c(values, FOC=as.character(del_genes[i,12]))
	    copySubstitute(inFile, outFile, values)
	    close(inFile)   
        }
        ## Create the Footer for the Table of Deleted Genes
        inFile <- system.file("template/FooterTable.html", package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, VAL="NO")
        copySubstitute(inFile, outFile, values)
        close(inFile)
    }
  
    if(nrow(amp_genes)>0){
        amp_genes <- as.matrix(amp_genes)
	uniq_gene <- unique(amp_genes[,1])
	## There are duplicate genes
	if(length(uniq_gene) != nrow(amp_genes)){
	    for(i in 1:length(uniq_gene)){
		pos <- which(amp_genes[,1]==uniq_gene[i])
		if(length(pos)>1){
		    pv <- round(min(as.numeric(amp_genes[pos,9]))+eps, 5)
		    amp_genes[pos[1],9] <- pv
				
		    fr <- as.character(amp_genes[pos,10])
		    fr <- substr(fr, 1, nchar(fr)-1)
		    amp_genes[pos[1],10] <- 
			paste(round(max(as.numeric(fr)), 1), "%", sep="")

		    m <- round(mean(as.numeric(amp_genes[pos,11]))+eps, 5)
		    amp_genes[pos[1],11] <- m
				
		    fs <- round(mean(as.numeric(amp_genes[pos,12]))+eps, 5)
		    amp_genes[pos[1],12] <- fs
		
		    amp_genes <- amp_genes[-pos[2:length(pos)],]
		}
		
	    }
	}

	amp_genes <-
	    amp_genes[order(as.numeric(substr(amp_genes[,10], 1,
	    nchar(amp_genes[,10])-1)), decreasing=TRUE),]

	 write.table(amp_genes, ampfile, sep="\t", col.names=TRUE, 
		row.names=FALSE, eol="\n", quote=FALSE) 

	amp_genes[,c(9,11,12)] <-
	format(round(as.numeric(amp_genes[,c(9,11,12)]), 5), decimal.mark=",")
        ## Create the Header for the Table of Amplified Genes
        inFile <- system.file("template/HeaderTableAberrantGenes.html",
                package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, TABNAME="Amplified", ABERR="Gain")
	if(correction==TRUE){
	    values <- c(values, SIGN="q")
	}else{
	    values <- c(values, SIGN="p")
	}
        copySubstitute(inFile, outFile, values)
        close(inFile)
	
        ## Insert a new Line for each Gene
        for(i in 1:nrow(amp_genes)){
        inFile <- system.file("template/AberrantGeneTableRow.html",
                package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(ENSID=as.character(amp_genes[i,1]))
        values <- c(values, GENEID=as.character(amp_genes[i,2]))
        values <- c(values, CHR=as.character(amp_genes[i,3]))
        values <- c(values, BPSTART=as.character(amp_genes[i,4]))   
        values <- c(values, BPEND=as.character(amp_genes[i,5]))
        if(is.na(amp_genes[i,6])){
            values <- c(values, CYTO="")
        }else{
            values <- c(values, CYTO=as.character(amp_genes[i,6]))
        }
        if(is.na(amp_genes[i,7])){
            values <- c(values, STRAND="")
        }else{
            values <- c(values, STRAND=as.character(amp_genes[i,7]))
        }
        if(is.na(amp_genes[i,8])){
            values <- c(values, DESCR="")
        }else{
            values <- c(values, DESCR=as.character(amp_genes[i,8]))
        }
        values <- c(values, ABERRPVAL=as.character(amp_genes[i,9]))
        values <- c(values, ABERRPRC=as.character(amp_genes[i,10]))
	values <- c(values, MEAN=as.character(amp_genes[i,11]))
	values <- c(values, FOC=as.character(amp_genes[i,12]))
        copySubstitute(inFile, outFile, values)
        close(inFile)   
        }
        ## Create the Footer for the Table of Amplified Genes
        inFile <- system.file("template/FooterTable.html", package="VegaMC")
        inFile <- file(inFile, "r")
        values <- c(values, VAL="NO")
        copySubstitute(inFile, outFile, values)
        close(inFile)
    }
    
    if(baf==TRUE){
        if(nrow(loh_genes)>0){
            loh_genes <- as.matrix(loh_genes)
	    loh_gene <- unique(loh_genes[,1])
	    ## There are duplicate genes
	    if(length(uniq_gene) != nrow(loh_genes)){
		for(i in 1:length(uniq_gene)){
		pos <- which(loh_genes[,1]==uniq_gene[i])
		if(length(pos)>1){
		    pv <- round(min(as.numeric(loh_genes[pos,9]))+eps, 5)
		    loh_genes[pos[1],9] <- pv
				
		    fr <- as.character(loh_genes[pos,10])
		    fr <- substr(fr, 1, nchar(fr)-1)
		    loh_genes[pos[1],10] <- 
			paste(round(max(as.numeric(fr)), 1), "%", sep="")

		    m <- round(mean(as.numeric(loh_genes[pos,11]))+eps, 5)
		    loh_genes[pos[1],11] <- m
				
		    fs <- round(mean(as.numeric(loh_genes[pos,12]))+eps, 5)
		    loh_genes[pos[1],12] <- fs
		
		    loh_genes <- loh_genes[-pos[2:length(pos)],]
		}
		
	    }
	    }

	    loh_genes <-
		loh_genes[order(as.numeric(substr(loh_genes[,10], 1,
			    nchar(loh_genes[,10])-1)), decreasing=TRUE),]

	    write.table(loh_genes, lohfile, sep="\t", col.names=TRUE, 
                    row.names=FALSE, eol="\n", quote=FALSE) 

	    loh_genes[,c(9,11,12)] <-
	format(round(as.numeric(loh_genes[,c(9,11,12)]), 5), decimal.mark=",")
            ## Create the Header for the Table of LOH Genes
            inFile <- system.file("template/HeaderTableAberrantGenes.html",
                    package="VegaMC")
            inFile <- file(inFile, "r")
            values <- c(values, TABNAME="LOH", ABERR="LOH")
	    if(correction==TRUE){
		values <- c(values, SIGN="q")
	    }else{
		values <- c(values, SIGN="p")
	    }
            copySubstitute(inFile, outFile, values)
            close(inFile)

            ## Insert a new Line for each Gene
            for(i in 1:nrow(loh_genes)){
                inFile <- system.file("template/AberrantGeneTableRow.html",
                    package="VegaMC")
                inFile <- file(inFile, "r")
                values <- c(ENSID=as.character(loh_genes[i,1]))
                values <- c(values, GENEID=as.character(loh_genes[i,2]))
                values <- c(values, CHR=as.character(loh_genes[i,3]))
                values <- c(values, BPSTART=as.character(loh_genes[i,4]))   
                values <- c(values, BPEND=as.character(loh_genes[i,5]))
                if(is.na(loh_genes[i,6])){
                    values <- c(values, CYTO="")
                }else{
                    values <- c(values, CYTO=as.character(loh_genes[i,6]))
                }
                if(is.na(loh_genes[i,7])){
                    values <- c(values, STRAND="")
                }else{
                    values <- c(values, STRAND=as.character(loh_genes[i,7]))
                }
                if(is.na(loh_genes[i,8])){
                    values <- c(values, DESCR="")
                }else{
                    values <- c(values, DESCR=as.character(loh_genes[i,8]))
                }
                values <- c(values, ABERRPVAL=as.character(loh_genes[i,9]))
                values <- c(values, ABERRPRC=as.character(loh_genes[i,10]))
		values <- c(values, MEAN=as.character(loh_genes[i,11]))
		values <- c(values, FOC=as.character(loh_genes[i,12]))
                copySubstitute(inFile, outFile, values)
                close(inFile)   
            }
            ## Create the Footer for the Table of LOH Genes
            inFile <- system.file("template/FooterTable.html",
                    package="VegaMC")
            inFile <- file(inFile, "r")
            values <- c(values, VAL="NO")
            copySubstitute(inFile, outFile, values)
            close(inFile)
        }
    }

    ## Create the Footer
    inFile <- system.file("template/Footer.html", package="VegaMC")
    inFile <- file(inFile, "r")
    values <- c(BIOINFO="http://bioinformatics.biogem.it")
    copySubstitute(inFile, outFile, values)
    close(inFile)
    close(outFile)

}

