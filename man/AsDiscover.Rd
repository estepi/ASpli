\name{AsDiscover}
\alias{AsDiscover}

\title{
Report PSI and PIR using experimental junctions
}
\description{
Given a bin, it is possible to calculate PSI/PIR metric using junctions to estimate changes in the use of it along different conditions.
}


\usage{
AsDiscover(counts, 
          targets, 
          features, 
          bam, 
          l,
          pair)
}
\arguments{
\item{counts}{
An object of class ASpli-counts.
}
\item{targets}{
A dataframe containing sample, bam and condition columns 
}
\item{features}{
An object of class ASpli-features.
}
\item{bam}{
A list with BAM files
}
\item{l}{
Read length of sequenced read. Default 100L
}
\item{pair}{
Vector of length two, either numeric or character, providing the pair of groups to be compared
}
}
\value{
An object of class ASpli-AS 
\item{irPIR}{reports: event, e1i counts (J1), ie1 counts (J2), j_within (J3), PIR by condition. J1, J2, J3 sum of junctions (J1, J2, J3) by condition. 
}

\item{altPSI}{ reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum of junctions (J1, J2, J3) by condition. 

}

\item{esPSI}{ reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum of junctions (J1, J2, J3) by condition. 

}

\item{junctionsPIR}{
PIR metric for each experimental junction using e1i and ie2 counts. Exclusion junction is the junction itself. This output helps to discover new introns as well as new retention events
}

\item{junctionsPSI}{
Given a junction, it is possible to analyze if it shares start, end or both with another junction. If so, is because there is more than one way for/of splicing. Using strand information it is possible to classify those pair of junctions into Alt5'ss, Alt3'ss or ES. Ratio between them along samples is reported.
}

}
\author{
Estefania Mancini, Marcelo Yanovsky and Ariel Chernomoretz 
}

\examples{
\dontrun{ 
library(ASpli)
library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome<-loadDb(chr14)
features<-binGenome(genome)
targets <- data.frame(row.names=names(RNAseqData.HNRNPC.bam.chr14_BAMFILES),
                      bam=matrix(unlist(RNAseqData.HNRNPC.bam.chr14_BAMFILES), ncol = 1, byrow = TRUE), 
                    condition=c(rep("CT",4),rep("KD",4)))
bam<-loadBAM(targets)
counts<-readCounts(features, bam, cores=8, l=100L, maxISize=50000)
group=factor(c(rep("CT",4),rep("KD",4)))
pair<-c("CT","KD")  
as<-AsDiscover(counts, targets, features, bam, l=100L, pair=pair)
writeAS(as=as, output.dir="only_as")
} 
}

\seealso{
\code{
Accesors: \link{irPIR}, \link{altPSI}, \link{esPSI}, \link{junctionsPIR},\link{junctionsPSI} 
Export: \link{writeAS}
}
}