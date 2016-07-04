\name{ASpliAS-class}
\Rdversion{1.1}
\docType{class}
\alias{ASpliAS-class}
\alias{altPSI,ASpliAS-method}
\alias{esPSI,ASpliAS-method}
\alias{irPIR,ASpliAS-method}
\alias{joint,ASpliAS-method}
\alias{junctionsPIR,ASpliAS-method}
\alias{junctionsPSI,ASpliAS-method}
\alias{writeAS,ASpliAS-method}
\title{Class \code{"ASpliAS"}}
\description{
Results of PSI and PIR using experimental junctions
}
\section{Slots}{
  \describe{
    \item{\code{irPIR}:}{Reports: event, e1i counts (J1), ie1 counts (J2), j_within (J3), PIR by condition. J1, J2, J3 sum of junctions (J1, J2, J3) by condition.}
    \item{\code{altPSI}:}{ Reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum of junctions (J1, J2, J3) by condition.}
    \item{\code{esPSI}:}{Reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum of junctions (J1, J2, J3) by condition.}
    \item{\code{join}:}{It is a combination of irPIR, altPSI and esPSI tables}
    \item{\code{junctionsPIR}:}{PIR metric for each experimental junction using e1i and ie2 counts. Exclusion junction is the junction itself. This output helps to discover new introns as well as new retention events}
    \item{\code{junctionsPSI}:}{Given a junction, it is possible to analyze if it shares start, end or both with another junction. If so, is because there is more than one way for/of splicing. Using strand information it is possible to classify those pair of junctions into Alt5'ss, Alt3'ss or ES. Ratio between them along samples is reported.}
    
  }
}
\author{
Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz
}
\seealso{\code{ 
Methods: \link{AsDiscover}, Accesors: \link{irPIR},\link{esPSI},\link{junctionsPIR},\link{junctionsPSI}
} 
}
