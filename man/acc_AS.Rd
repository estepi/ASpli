\name{AS accesors}
\alias{altPSI}
\alias{esPSI}
\alias{irPIR}
\alias{joint}
\alias{junctionsPIR}
\alias{junctionsPSI}
\title{
Accessors for ASpliAS object
}
\description{
Accessors for ASpliAS object
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
  An ASpliAS object 
  }
}
\value{
Returns dataframes with genomic metadata and  PSI and PIR metrics
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
group <- factor(c(rep("CT",4),rep("KD",4)))
pair <- c("CT","KD")  
as <- AsDiscover(counts, targets, features, bam, threshold=5, l=100, pair=pair)
altPSI(as)
esPSI(as)
irPIR(as)
joint(as)
junctionsPIR(as)
junctionsPSI(as)
}

