\name{AS accesors}
\alias{altPSI}
\alias{esPSI}
\alias{irPIR}
\alias{joint}
\alias{junctionsPIR}
\alias{junctionsPSI}
\title{
Accessors for ASpli-AS object
}
\description{
Accessors for ASpli-AS object
}
\usage{
altPSI(x)
esPSI(x)
irPIR(x)
joint(x)
junctionsPIR(x)
junctionsPSI(x)
}
\arguments{
  \item{x}{   
  An ASpli-AS object 
  }
}
\value{
Returns dataframes with genomic metadata and  PSI and PIR metrics
}
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }
\examples{
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome<-loadDb(chr14)
features<-binGenome(genome)
library(RNAseqData.HNRNPC.bam.chr14)
targets <- data.frame(row.names=names(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
                      bam=matrix(unlist(RNAseqData.HNRNPC.bam.chr14_BAMFILES), 
                                 ncol = 1, byrow = TRUE), 
                      condition=c(rep("CT",4),rep("KD",4)))
bam<-loadBAM(targets)
counts<-readCounts(features, bam,  l=100L, maxISize=50000);
group<-factor(c(rep("CT",4),rep("KD",4)))
pair<-c("CT","KD")  
as<-AsDiscover(counts, targets, features, bam, threshold=5, l=100, pair=pair)
altPSI(as)
esPSI(as)
irPIR(as)
joint(as)
junctionsPIR(as)
junctionsPSI(as)
}

