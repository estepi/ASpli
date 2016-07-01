\name{Counts accesors}
\alias{countsb}
\alias{countse1i}
\alias{countsg}
\alias{countsie2}
\alias{countsj}
\alias{rdsg}
\alias{rdsb}
\title{
Accessors for ASpliCounts object
}
\description{
Accessors for ASpliCounts object
}
\usage{
countsb(x)
countse1i(x)
countsg(x)
countsie2(x)
countsj(x)
rdsg(x)
rdsb(x)
}
\arguments{
  \item{x}{   
  An ASpliCounts object 
  }
}
\value{
Returns dataframes with counts by sample and genomic metadata
}
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }
\examples{
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome <- loadDb(chr14)
features <- binGenome(genome)
library(RNAseqData.HNRNPC.bam.chr14)
targets <- data.frame(bam=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                       condition=c(rep("CT",4),rep("KD",4)))
bam <- loadBAM(targets)
counts <- readCounts(features, bam,  l=100L, maxISize=50000);
countsb(counts)
countse1i(counts)
countsg(counts)
countsie2(counts)
countsj(counts)
rdsg(counts)
rdsb(counts)
}


