\name{DunedinPACE}
\alias{DunedinPACE}
\title{Dunedin PACE Clock Data}
\description{
This dataset is utilized for calculating the Dunedin PACE epigenetic age.
It is accessed when the \code{calculateDunedinPACE} function is called.
Part of dnaMethyAge package.
\href{https://link.springer.com/article/10.1007/s11357-023-00871-w}{Wang et
al., 2023}
}
\usage{
data("DunedinPACE", package = "EpigeneticAgePipeline")
}
\format{
An R data object.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("DunedinPACE", package = "EpigeneticAgePipeline", envir = myEnv)
print(myEnv$coefs)
print(myEnv$gold_standard_means)
}
