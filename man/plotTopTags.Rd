\name{plotTopTags}
\alias{plotTopTags}

\title{
Coverage plots 
}
\description{
Using genomic coordinates and BAM files this function is useful for make coverage plots
}
\usage{
plotTopTags(auxdf, genome, targetsPlot, output.dir)
}

\arguments{
\item{auxdf}{ A data frame: row.naMes=bin names, gene coordinates, bin  coordinates and  event name columns
}
\item{genome}{
TxDb genome
}
\item{targetsPlot}{
A dataframe containting: bam files name, condition (y axe tag), color for each condition}
\item{output.dir}{ Name of directory where plots are supossed to be exported }
}

\value{
Coverage plots in png format of selected events
}
\author{
Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome <- loadDb(chr14)
features <- binGenome(genome) 
targets <- data.frame(bam=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                       condition=c(rep("CT",4),rep("KD",4)))
bam <- loadBAM(targets)
counts <- readCounts(features, bam, l=100L, maxISize = 50000)
pair <- c("CT","KD")
group <- c(rep("CT", 4),rep("KD", 4))
du_HNRNPC <- DUreport(counts, targets, pair, group)  
bins <- binsDU(du_HNRNPC)
topTagsBins <- which(bins$bin.fdr <= 0.1 & 
                     abs(bins$logFC) >=0.58)
targetsPlot <- data.frame(bam=targets$bam, 
                        sample=targets$condition, 
                        color=c(rep("blue", 4),rep("red", 4)), 
                        stringsAsFactors=FALSE)

auxdf<-bins[topTagsBins,]
#for simplicity, just one: LRR1:E005

plotTopTags(auxdf["LRR1:E005",], 
            genome, 
            targetsPlot, 
            output.dir="testPlots")

}


