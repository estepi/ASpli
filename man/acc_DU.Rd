\name{DU accesors}
\alias{genesDE}
\alias{binsDU}
\alias{junctionsDU}
\title{
Accessors for ASpliDU object
}
\description{
Accessors for ASpliDU object
}
\usage{
genesDE(x)
binsDU(x)
junctionsDU(x)

}
\arguments{
  \item{x}{   
  An ASpliDU object 
  }
}
\value{
Returns dataframes with genomic metadata and logFC and pvalue
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
du <- DUreport(counts, targets, pair, group)
genesDE(du)
binsDU(du)
junctionsDU(du)
}


