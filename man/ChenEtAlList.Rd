\name{ChenEtAlList}
\alias{ChenEtAlList}
\title{List of Methylation Probes from Chen et al.}
\description{
This dataset includes a list of methylation probes from the study by
\href{https://doi.org/10.4161/epi.23470}{Chen et al., 2013}.
It is used to remove cross reactive probes during IDAT file processing.
}
\usage{
data("ChenEtAlList", package = "EpigeneticAgePipeline")
}
\format{
A data frame where each row represents a methylation probe.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("ChenEtAlList", package = "EpigeneticAgePipeline", envir = myEnv)
head(myEnv$ChenEtAlList)
}
