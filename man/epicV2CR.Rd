\name{epicV2CR}
\alias{epicV2CR}
\title{Cross Reactive Probes for Illumina EPICv2 Methylation Array}
\description{
A list of methylation probes. Used to remove cross reactive probes during
IDAT file processing.
\href{https://pmc.ncbi.nlm.nih.gov/articles/PMC10916044/#MOESM4}{Peters et al., 2012}
}
\usage{
data("epicV2CR", package = "EpigeneticAgePipeline")
}
\format{
A data frame where each row represents a methylation probe.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("epicV2CR", package = "EpigeneticAgePipeline",
    envir = myEnv)
head(myEnv$epicV2CR)
}

