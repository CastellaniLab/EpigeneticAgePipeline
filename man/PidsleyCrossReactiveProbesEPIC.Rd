\name{PidsleyCrossReactiveProbesEPIC}
\alias{PidsleyCrossReactiveProbesEPIC}
\title{Cross-Reactive Probes List for EPIC Arrays}
\description{
A list of methylation probes. Used to remove cross reactive probes during
IDAT file processing. For EPIC microarray's.
\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/}{Pidsley
et al., 2016}
}
\usage{
data("PidsleyCrossReactiveProbesEPIC", package = "EpigeneticAgePipeline")
}
\format{
A data frame where each row represents a methylation probe.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("PidsleyCrossReactiveProbesEPIC", package = "EpigeneticAgePipeline",
    envir = myEnv)
head(myEnv$PidsleyCrossReactiveProbesEPIC)
}
