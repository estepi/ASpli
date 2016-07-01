\name{DUreport_DEXSeq}
\alias{DUreport_DEXSeq}
\title{
Differential gene expression and differential bin/junction usage estimation
}
\description{
Estimate differential expression at gene level and differential usage at bin and junction level.
}
\usage{
DUreport_DEXSeq(counts, targets, pair, group, minGenReads, minBinReads, minRds,threshold)
}

\arguments{
\item{counts}{
An object of class ASpliCounts
}
\item{targets}{
A dataframe containing sample, bam and condition columns 
}
\item{pair}{
vector of length two, either numeric or character, providing the pair of groups to be compared
}
\item{group}{Factorial vector with tags for each sample}
\item{minGenReads}{Default 10 reads}
\item{minBinReads}{Default 5 reads}
\item{minRds}{Default 0.05}
\item{threshold}{Minimun number of junction. Default 5}
}


\value{ An ASpliDU object with results at \code{genes}, \code{bins} and \code{junctions} level }
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
counts <- readCounts(features, bam, l=100L, maxISize=50000)
group <- factor(c(rep("CT",4), 
                rep("KD",4)))
pair <- c("CT","KD")  
du <- DUreport_DEXSeq(counts, targets, pair, group)
writeDU(du, output.dir="only_du")
}

\seealso{
\code{
\link{DEXSeq}, \link{edgeR}
Accesors: \link{genesDE}, \link{binsDU},\link{junctionsDU}
Export: \link{writeDU}
}
}


