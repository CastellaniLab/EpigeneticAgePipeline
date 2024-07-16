\name{coefs}
\alias{coefs}
\title{Coefficient Data}
\description{
Epigenetic age probe data. Part of dnaMethyAge package.
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
print(myEnv$coefs)
}
