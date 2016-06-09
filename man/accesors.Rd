\name{Accesors}
\alias{featuresg}
\alias{featuresb}
\alias{featuresj}
\alias{countsg}
\alias{countsb}
\alias{countsj}
\alias{countse1i}
\alias{countsie2}
\alias{rdsb}
\alias{rdsg}
\alias{irPIR}
\alias{altPSI}
\alias{esPSI}
\alias{joint}
\alias{junctionsPIR}
\alias{junctionsPSI}
\alias{genesDE}
\alias{binsDU}
\alias{junctionsDU}
\title{
Accessors for ASpli objects for every feature
}
\description{ Extract ASpli objects slots }
\usage{
  GeneRanges<-featuresg(features)
  BinRanges<-featuresb(features)
  JunctionRanges<-featuresj(features)
  GeneCounts<-countsg(counts)
  BinCounts<-countsb(counts)
  JunctionCounts<-countsj(counts)
  E1iCounts<-countse1i(counts)
  IE2Counts<-countsie2(counts)
  geneRds<-rdsg(counts)
  binRds<-rdsb(counts)
  JunctionsIrPIR<-irPIR(as)
  JunctionsAltPSI<-altPSI(as)
  JunctionsEsPSI<-esPSI(as)
  junctionsPIR<-junctionsPIR(as)
  junctionsPSI<-junctionsPSI(as)
  genesDE<-genesDE(du)
  binsDU<-binsDU(du)
  junctionsDU<-junctionsDU(du)
  binJoin<-joint(du)
}

\arguments{
  \item{features}{ An ASpli-features object }
  \item{counts}{ An ASpli-counts object }
  \item{as}{ An ASpli-AS object }
  \item{du}{ An ASpli-DU object } 
}
\details{ Extract dataframes from ASpli objects }
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }
\examples{
library(ASpli)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome<-loadDb(chr14)
features<-binGenome(genome)
GeneRanges<-featuresg(features)
BinRanges<-featuresb(features)
JunctionRanges<-featuresj(features)
}

