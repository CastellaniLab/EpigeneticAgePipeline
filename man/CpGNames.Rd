\name{CpGNames}
\alias{CpGNames}
\title{CpG Name Array}
\description{
Character array of CpG's, used in generation of test data.
}
\usage{
data("CpGNames", package = "EpigeneticAgePipeline")
}
\format{
An R data object.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("CpGNames", package = "EpigeneticAgePipeline", envir = myEnv)
print(myEnv$CpGNames)
}
