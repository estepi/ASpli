\name{readCounts}
\alias{readCounts}
\title{
Summarize read overlaps
}
\description{
Summarize read overlaps against all feature levels }
\usage{
readCounts(features, bam, cores, l, maxISize, minAnchor)
}
\arguments{
\item{features}{
An object of class ASpliFeatures. It is a list of GRanges at gene, bin and junction level
}
\item{bam}{
List of bam files
}
\item{l}{
Read length of sequenced library. It is used for compute E1I and IE2 read summarization
}
\item{maxISize}{maximum intron expected size. Junctions longer than this size will be dicarded  }
\item{cores}{Number of cores to use. Default 1}
\item{minAnchor}{ Percentage of read thath sould be aligned in exon-intron boundary  }
}

\value{
An object of class ASpliCounts. Each slot is a dataframe containing features metadata  and read counts. Summarization is reported at gene, bin, junction and intron flanking regions (E1I, IE2)
}
\author{
Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz 


}
\seealso{
\code{
Accesors: \link{countsg}, 
\link{countsb}, \link{countsj},\link{countse1i}, \link{countsie2},\link{rdsg},\link{rdsb} Export: \link{writeCounts}}
}
\examples{
library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome <- loadDb(chr14)
features <- binGenome(genome)
targets <- data.frame(bam=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                       condition=c(rep("CT",4),rep("KD",4)))
bam <- loadBAM(targets)
counts <- readCounts(features, bam, l=100L, maxISize=50000)#OK
writeCounts(counts,output.dir="only_counts")
}

