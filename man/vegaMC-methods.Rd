\name{vegaMC-methods}
\docType{methods}
\alias{vegaMC}
\alias{vegaMC-methods}
\alias{vegaMC,character-method}
\alias{vegaMC}
\title{
    Class aggregator of VegaMC.
}

\description{
    VegaMC enables the detection of driver chromosomal imbalances 
    (deletions, amplifications and loss of heterozygosities (LOHs)) from 
    array comparative genomic hybridization (aCGH) data. VegaMC performs 
    a joint segmentation of aCGH data. Segmented regions are then used 
    into a statistical framework to distinguish between driver and 
    passenger mutations. In this way, significant imbalances can be 
    detected by the associated p-value. VegaMC has been implemented to be 
    easily integrated with the output produced by PennCNV and with the Genoset 
    eSet Objects. VegaMC produces in output two web pages allowing a rapid 
    navigation between both detected regions and altered genes. In the web page 
    summarizing the altered genes, the user finds the link to the respective 
    Ensembl gene web page.
}
\usage{
    vegaMC(dataset, output_file_name="output", beta=0.5, 
    min_region_bp_size=1000, correction=FALSE, 
    loss_threshold=-0.2, gain_threshold=0.2, 
    baf=TRUE, loh_threshold=0.75, loh_frequency=0.8, bs=1000, 
    pval_threshold=0.05, html=TRUE, getGenes=TRUE, 
    mart_database="ensembl", ensembl_dataset="hsapiens_gene_ensembl")
}

\section{Methods}{
	\describe{
		\item{\code{signature(dataset = "character")}}{
			This method allows to run VegaMC on a data file in PennCNV format. 
		}

		\item{\code{signature(dataset = "GenoSet")}}{
			This method allows to run VegaMC on a GenoSet object of 
			\code{genoset} package. 
		}
	}
}

\arguments{

    \item{dataset}{Dataset file following the PennCNV format: The first three columns 
    describe the name, the chromosome and the position respectively. 
    The other columns of the matrix report the LRR and the BAF 
    (if available) of each sample. Note that observations must be ordered by the
    respective genomic position.}

    \item{output_file_name}{(Default code{output}) File name used 
    to save the results.}

    \item{beta}{(Default 0.5) This parameter is used to compute the 
    stop condition. It is used to calculate the maximum jump allowed in 
    scale parameter updating. If \code{beta}=0 then the resulting 
    segmentation will be composed of a region for each probe (all regions 
    will contain just a probe). In contrast, if \code{beta} is very large, 
    then the segmentation will contain just a region for each chromosome.}

    \item{min_region_bp_size}{(Default 1000) VegaMC deletes from the list 
    the regions shorter then this size (in bp).}

    \item{correction}{(Default FALSE) If this parameter is TRUE multiple 
    testing corrections is performed.}

    \item{loss_threshold}{(Default -0.2) Values used to mark a region as a 
    deletion (loss). If the wighted mean of a region is lower than this 
    threshold, then the region is marked as a deletion (loss).}

    \item{gain_threshold}{(Default 0.2) Values used to mark a region as an 
    amplification (gain). If the wighted mean of a region is greater than 
    this threshold, then the region is marked as an amplification (gain).}

    \item{baf}{(Default \code{TRUE}) This parameter specifies if the 
    dataset contains BAF measurements (default \code{TRUE}. If BAF is 
    available, then VegaMC is able to compute LOH imbalances.}

    \item{loh_threshold}{(Default 0.75) Threshold used to distinguish 
    between homozygous and heterozygous genotypes. If the BAF is greater 
    than \code{loh_threshold} or lower then (1-\code{loh_threshold}) 
    then the respective probe is considered to be homozygous.}

    \item{loh_frequency}{(Default 0.8) Minimum fraction of homozygous 
    probes needed for marking a region as LOHs. Regions with a fraction 
    of homozygous probes greater than this threshold are marked as LOH.}

    \item{bs}{(Default 1000) Number of permutation bootstraps performed 
    to compute the null distribution.}

    \item{pval_threshold}{(Default 0.05) Significance level used to 
    reject the null hypothesis. If the p-value of an aberration 
    (loss, gain, LOH) is not greater than this threshold, then the region 
    is considered to be significant and, consequently, it is considered a 
    driver mutation.}

    \item{html}{(Default \code{TRUE}) If this value is \code{TRUE}, then 
    in output will be produced an html file called 
    \code{output_file_name}.html in which a summary of all detected regions 
    is reported.}

    \item{getGenes}{(Default \code{TRUE}) If this value is \code{TRUE}, 
    then in output will be produced an html file called 
    \code{output_file_name}Genes.html in which the list of all genes 
    overlapping the significant regions is reported.}

    \item{mart_database}{(Default \code{ensembl}) BioMart database name 
    you want to connect to. Possible database names can be retrieved 
    with the function listMarts of biomaRt package.}

    \item{ensembl_dataset}{(Default \code{hsapiens_gene_ensembl}) BioMart 
    dataset used to get information from Ensembl BioMart database. 
    In order to obtain the list of all available dataset use the command 
    \code{listDatasets} of biomaRt package.}
}


\value{
    After the execution of this function, a matrix containing all 
    information on the detected regions is returned. This object is a 
    matrix having a row for each detected regions described by the 
    following columns:

    \item{Chromosome}{The chromosome in which the region is located.}
    \item{bp Start}{The position in which the region starts (in bp).}
    \item{bp End}{The position in which the region ends (in bp).}
    \item{Region Size}{The size of the regions (in bp).}
    \item{Mean}{The weighted mean of the region computed on all samples.}
    \item{Loss p-value}{The p-value associated to the probability to have 
    a driver deletion.}
    \item{Gain p-value}{The p-value associated to the probability to have 
    a driver amplification.}
    \item{LOH p-value}{The p-value associated to the probability to have 
    a driver LOH.}
    \item{\% Loss}{The percentage of samples showing a deletion for this 
    region.}
    \item{\% Gain}{The percentage of samples showing an amplification for 
    this region.}
    \item{\% LOH}{The percentage of samples showing a LOH for this region.}
    \item{Probe Size}{The number of probes composing the region.}
    \item{Loss Mean}{Mean of LRR computed only on the samples that show a loss.}
    \item{Gain Mean}{Mean of LRR computed only on the samples that show a gain.}
    \item{LOH Mean}{Mean of LRR computed only on the samples that show a LOH.}
    \item{Focal-score Loss}{Focal Score associated to deletion.}
    \item{Focal-score Gain}{Focal Score associated to amplification.}
    \item{Focal-score LOH}{Focal Score associated to LOH.}
    This matrix is automatically saved in the current work directory as a 
    tab delimited file. For default the name used to asave the file is 
    'output'.
}
\keyword{methods}
\references{
    Morganella S., and Ceccarelli M. VegaMC: a R/bioconductor package for fast downstream analysis of large array comparative genomic hybridization datasets. Bioinformatics, 28(19):2512-4 (2012).
}

\author{
    Sandro Morganella
}
\examples{
    ## Run VegaMC 
    file.copy(system.file("example/breast_Affy500K.txt", package="VegaMC"), ".")
    results <- vegaMC("breast_Affy500K.txt", "results", html=FALSE, 
                        getGenes=FALSE)
                                      
}
