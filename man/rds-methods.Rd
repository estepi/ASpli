\name{rds-methods}
\docType{methods}
\alias{rds-methods}
\alias{rds,ASpli-counts-method}
\title{ Methods for Function \code{rds}}
\description{ Divide number of reads by length of the features}
\usage{
counts<-rds(counts, targets)
}
\arguments{
    \item{counts}{An ASpli-counts object }
    \item{targets}{ targets dataframe }
}
\section{Methods}{
\describe{
\item{\code{signature(counts = "ASpli-counts")}}{
}
}}
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }
\seealso{
\code{
Accesors: \link{rdsg},\link{rdsb} Export:\link{writeRds}}
}

\examples{
library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome<-loadDb(chr14)
features<-binGenome(genome)
targets <- data.frame(row.names=names(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
                      bam=matrix(unlist(RNAseqData.HNRNPC.bam.chr14_BAMFILES), ncol = 1, byrow = TRUE), 
                      condition=c(rep("CT",4),rep("KD",4)))

bam<-loadBAM(targets)
counts<-readCounts(features, bam, l=100L, maxISize=50000)#OK
countsRds<-rds(counts, targets)
writeRds(counts, output.dir="rds")
}



