\name{DUreport}
\alias{DUreport}

\title{
Differential gene expression and differential bin/junction usage estimation
}
\description{
Estimate differential expression at gene level and differential usage at bin and junction level.
}
\usage{
DUreport(counts, targets, pair, group, minGenReads, minBinReads, minRds,ignoreExternal,threshold)
}

\arguments{
\item{counts}{
An object of class ASpli-counts
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
\item{ignoreExternal}{Ignore Exon Bins at the beggining or end of the transcript. Default TRUE}
\item{threshold}{Minimun number of junction. Default 5}
}


\author{
Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz
}
\examples{
\dontrun{ 
library(ASpli)
library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome<-loadDb(chr14)
features<-binGenome(genome)
targets <- data.frame(row.names=names(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
                      bam=matrix(unlist(RNAseqData.HNRNPC.bam.chr14_BAMFILES), 
                                 ncol = 1, byrow = TRUE), 
                      condition=c(rep("CT",4),rep("KD",4)))
bam<-loadBAM(targets, cores=8)
counts<-readCounts(features, bam, cores=8, l=100L, maxISize=50000)
group<-factor(c(rep("CT",4), 
                rep("KD",4)))
pair<-c("CT","KD")  
du<-DUreport(counts, targets, pair, group)
writeDU(du, output.dir="only_du")
} 

}
\seealso{
\code{
Accesors: \link{genesDE}, \link{binsDU},\link{junctionsDU}
Export: \link{writeDU}
}
}

