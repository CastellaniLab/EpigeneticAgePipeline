\name{non-specific-probes-Illumina27k}
\alias{non-specific-probes-Illumina27k}
\title{Non-Specific Probes for Illumina 27k Methylation Array}
\description{
A list of methylation probes. Used to remove cross reactive probes during
IDAT file processing.
\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484480/}{Chen et al., 2012}
}
\usage{
data("non-specific-probes-Illumina27k")
}
\format{
A data frame where each row represents a methylation probe.
}
\examples{
data("non-specific-probes-Illumina27k", package = "EpigeneticAgePipeline")
head(non-specific-probes-Illumina27k)
}

