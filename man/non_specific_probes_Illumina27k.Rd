\name{non_specific_probes_Illumina27k}
\alias{non_specific_probes_Illumina27k}
\title{Cross Reactive Probes for Illumina 27k Methylation Array}
\description{
A list of methylation probes. Used to remove cross reactive probes during
IDAT file processing.
\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484480/}{Chen et al., 2012}
}
\usage{
data("non_specific_probes_Illumina27k", package = "EpigeneticAgePipeline")
}
\format{
A data frame where each row represents a methylation probe.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("non_specific_probes_Illumina27k", package = "EpigeneticAgePipeline",
    envir = myEnv)
head(myEnv$non_specific_probes_Illumina27k)
}

