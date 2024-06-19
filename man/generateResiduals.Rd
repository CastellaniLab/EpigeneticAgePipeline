\name{generateResiduals}
\alias{generateResiduals}

\title{generateResiduals}

\description{
A function used to perform principal component analysis and residual generation
post epigenetic age measure generation.
}

\usage{
generateResiduals(directory = getwd(),
useBeta = FALSE,
arrayType = "450K")

}

\arguments{
\item{directory}{
Directory containing input data files (default: current working directory).
}

\item{useBeta}{
Logical. If TRUE,
expect a .csv containing beta values
(scaled between 0 and 1).
If FALSE, process raw intensity data (IDAT).
}

\item{arrayType}{
Type of DNA methylation array used (options: "27K", "450K", or "EPIC").
}

}

\details{
Input Files:

Sample_Sheet.csv:
.csv file containing phenotypic data for each sample
*See note for conventions

IDAT Files:
IDAT files containing methylation data for each sample

betaValues.csv:
If IDAT files are not available, processed beta values can be provided

}

\value{
The function will return the follwing:

Residuals.csv: A .csv file containing sample residual data, with your epigenetic
age measure against your independent variables.

Outliers.csv: A .csv file containing samples that were deemed to be outliers
during PCA.
}


\note{
Naming Conventions for Data in "Sample_Sheet.csv":

When using the generateResiduals function, name the column with epigenetic age
values "EpiAge@@@1".

To include a variable from Sample_Sheet.csv:
If a variable contains numeric type data, append "@@@1" to the column name

If a variable contains factor type data, append "@@@2" to the column name.

Ex. variable "isSmoker" would become "isSmoker@@@2".

To generate GrimAge, chronological age must be included ("Age@@@1"),
as well as sample sex ("Sex@@@2"). Valid values for male sex would be "1",
"M" or "Male". Valid values for female sex would be "2", "F", or "Female".

For a more detailed guide, visit the github page:
https://github.com/CastellaniLab/EpigeneticAgePipeline

}

\examples{
#Example of using generateResiduals function
library(EpigeneticAgePipeline)
directory <- paste0(path.package("EpigeneticAgePipeline"),"/data/")
generateResiduals(
    directory = directory,
    useBeta = TRUE,
    arrayType = "450K",
)
message("Outliers")
knitr::kable(outliersCSV)
message("Residuals")
knitr::kable(residualsCSV)
message("Acceleration Residuals")
knitr::kable(accelResidualsCSV)
}



