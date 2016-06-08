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

\examples{\dontrun{ library(ASpli)
library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome<-loadDb(chr14)
features<-binGenome(genome) 
targets <- data.frame(row.names=names(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
                        bam=matrix(unlist(RNAseqData.HNRNPC.bam.chr14_BAMFILES), ncol = 1, byrow = TRUE), 
                        condition=c(rep("CT",4),rep("KD",4)))
bam<-loadBAM(targets, cores=8)
counts<-readCounts(features, bam, cores=8, l=100L)
pair<-c("CT","KD")
group<-c(rep("KD", 4),rep("CT", 4))
du_HNRNPC<-DUreport(counts, targets, pair, group)  
topTagsBins<-which(du_HNRNPC@bins$bin.fdr <= 0.1 & 
                       abs(du_HNRNPC@bins$logFC) >=0.58)
targetsPlot<-data.frame(bam=c("CT.bam","KD.bam"), 
                          sample=c("ct", "kd"), 
                          color=c("blue", "black"), 
                          stringsAsFactors=FALSE)
auxdf<-du_HNRNPC@bins[topTagsBins[1:10],c(2,6,7,8,9,10)]; head(auxdf)
library(Gviz)
plotTopTags(auxdf, 
            genome, 
            targetsPlot, 
            output.dir="testPlots")
    }}


