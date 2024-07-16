\name{PCGrimAge_agesex}
\alias{PCGrimAge_agesex}
\title{GrimAge Data}
\description{
This dataset is utilized for calculating GrimAge epigenetic age acceleration.
Part of dnaMethyAge package.
\href{https://link.springer.com/article/10.1007/s11357-023-00871-w}{Wang et
al., 2023}
}
\usage{
data("PC-clocks", package = "EpigeneticAgePipeline")
}
\format{
An R data object.
}
\examples{
myEnv <- new.env(parent = emptyenv())
data("PC-clocks", package = "EpigeneticAgePipeline", envir = myEnv)
print(myEnv$PCGrimAge_agesex)
}
