\name{DU accesors}
\alias{genesDE}
\alias{binsDU}
\alias{junctionsDU}
\title{
Accessors for ASpli-DU object
}
\description{
Accessors for ASpli-DU object
}
\usage{
genesDE(x)
binsDU(x)
junctionsDU(x)

}
\arguments{
  \item{x}{   
  An ASpli-DU object 
  }
}
\value{
Returns dataframes with genomic metadata and logFC and pvalue
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
du<-DUreport(counts, targets, pair, group)
genesDE(du)
binsDU(du)
junctionsDU(du)
}


