\name{ASpliCounts-class}
\docType{class}
\alias{ASpliCounts-class}
\alias{AsDiscover,ASpliCounts-method}
\alias{countsb,ASpliCounts-method}
\alias{countse1i,ASpliCounts-method}
\alias{countsg,ASpliCounts-method}
\alias{countsie2,ASpliCounts-method}
\alias{countsj,ASpliCounts-method}
\alias{DUreport_DEXSeq,ASpliCounts-method}
\alias{DUreport,ASpliCounts-method}
\alias{rdsb,ASpliCounts-method}
\alias{rdsg,ASpliCounts-method}
\alias{rds,ASpliCounts-method}
\alias{writeCounts,ASpliCounts-method}
\alias{writeRds,ASpliCounts-method}
\title{Class \code{"ASpliCounts"}}
\description{ Contains results of read overlaps against all feature levels summarization }
\section{Slots}{
  \describe{
    \item{\code{gene.counts}:}{Object of class \code{"data.frame"}  }
    \item{\code{exon.intron.counts}:}{Object of class \code{"data.frame"}  }
    \item{\code{junction.counts}:}{Object of class \code{"data.frame"} }
    \item{\code{e1i.counts}:}{Object of class \code{"data.frame"} }
    \item{\code{ie2.counts}:}{Object of class \code{"data.frame"}  }
    \item{\code{gene.rd}:}{Object of class \code{"data.frame"} }
    \item{\code{bin.rd}:}{Object of class \code{"data.frame"} }
  }
}
\section{Methods}{
  \describe{
    \item{AsDiscover}{ psi and pir metrics}
    \item{countsb}{bin counts accesor}
    \item{countse1i}{ e1i counts accesor}
    \item{countsg}{ gene counts accesor}
    \item{countsie2}{ie2 counts accesor}
    \item{countsj}{junction counts accesor}
    \item{DUreport_DEXSeq}{ differential expression and usage estimation using DEXSeq }
    \item{DUreport}{differential expression and usage estimation using DEXSeq}
    \item{rdsb}{bin read densities accesor}
    \item{rdsg}{gen read densities acceesor}
    \item{rds}{compute read densities on genes and bins}
    \item{writeCounts}{Export count tables}
    \item{writeRds}{Export read density tables}
	 }
}
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz}
