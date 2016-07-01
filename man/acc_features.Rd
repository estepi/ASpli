\name{features accesors}
\alias{featuresb}
\alias{featuresg}
\alias{featuresj}
\title{
Accessors for ASpliFeatures object
}
\description{
Accessors for ASpliFeatures object
}
\usage{
featuresg(x)
featuresb(x)
featuresj(x)
}
\arguments{
  \item{x}{
  An ASpliFeatures object
}
}
\value{
Resturns a GenomicRanges object
}
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }
\examples{
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome <- loadDb(chr14)
features <- binGenome(genome)
featuresg(features)
featuresb(features)
featuresj(features)
}


