\name{write}
\alias{writeCounts}
\alias{writeRds}
\alias{writeAS}
\alias{writeDU}
\alias{writeAll}
\title{
Write results 
}
\description{
Export tab delimited files in structured output
}
\usage{
writeCounts(counts, output.dir="counts")
writeRds(counts, output.dir="rds")
writeDU(du, output.dir="du")
writeAS(as, output.dir="as")
writeAll(counts, du, as, output.dir="output")
}
\arguments{
\item{counts}{
An ASpli-counts object
}
\item{as}{
An ASpli-AS object
}
\item{du}{
An ASpli-DU object
}
\item{output.dir}{
Name of output folder (new or existing)
}
}
\details{
Tab delimited files are exported in a tidy manner into output folder
}
\author{
Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz 
}
