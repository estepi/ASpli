\name{loadBAM}
\alias{loadBAM}
\title{
Load BAM files
}

\description{
Load BAM files into R session using targets object especification 
}
\usage{loadBAM(targets, cores)}
\arguments{
\item{targets}{A dataframe containing sample, bam and condition columns }
 \item{cores}{Number of proccesors to use}
}
\value{
A list of GAlignments. Each element of the list correspond to a BAM file (or sample)
}
\author{
Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz 
}

\examples{
library(RNAseqData.HNRNPC.bam.chr14)
targets <- data.frame(bam=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                       condition=c(rep("CT",4),rep("KD",4)))
targets
bam <- loadBAM(targets)

}
