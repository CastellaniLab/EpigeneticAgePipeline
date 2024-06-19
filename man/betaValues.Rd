\name{betaValues}
\alias{betaValues}
\title{Beta Values for Methylation Sites}
\description{
This dataset contains the beta values for methylation sites across samples.
It is used when the parameter \code{useBeta} is set to TRUE in the main
function, allowing direct reading of processed beta values
without needing to process IDAT files. Publically available data from
GEO database GSE197646.
}
\usage{
data(betaValues)
}
\format{
A data frame with rows as methylation sites and columns as samples.
}
\examples{
data(betaValues, package = "EpigeneticAgePipeline")
head(betaValues)
}
