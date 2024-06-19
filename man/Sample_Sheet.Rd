\name{Sample_Sheet}
\alias{Sample_Sheet}
\title{Sample Metadata Sheet}
\description{
This CSV file contains covariate data for each sample included in the
analysis. Contains publically available data from GEO database GSE197646.
}
\usage{
data(Sample_Sheet)
}
\format{
A data frame with rows as samples and columns containing covariate data
such as `Sample_ID`, `Age`, `Sex`, and other relevant information.
}
\examples{
data(Sample_Sheet, package = "EpigeneticAgePipeline")
head(Sample_Sheet)
}
