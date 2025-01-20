\name{generateResiduals}
\alias{generateResiduals}

\title{Generate Epigenetic Age Residual Data}

\description{
A function used to perform principal component analysis and residual generation
post epigenetic age measure generation.
}

\usage{
generateResiduals(directory = getwd(),
useBeta = FALSE,
formula = NULL,
arrayType = "450K",
ignoreCor = TRUE,
PCs = 5,
threshold = 3,
doParallel = TRUE,
doCellCounts = TRUE,
useAdult = FALSE)

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

\item{formula}{
String. If the user wants to use a custom formula, one may input it as a string.
If using PC's, input in the format of "PC1".
Adult blood cell names: "Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK"
Cord blood cell names: "Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC", "NK"
Example:
"myAge ~ Age + BMI + Smoking_Status + PC1"
}

\item{arrayType}{
Type of DNA methylation array used (options: "27K", "450K", or "EPIC").
}

\item{ignoreCor}{
Logical. If TRUE then removing highly correlated (>=0.6) covariates will not be an
option.
}

\item{PCs}{
Numerical. Number of PC's to generate.
}

\item{threshold}{
Numerical. Threshold for MAD outlier detection.
median +- threshold * MAD
}

\item{doParallel}{
Logical. If TRUE, use parallel processing if reading beta values.
}

\item{doCellCounts}{
Logical. If TRUE include cell counts for analysis.
}

\item{useAdult}{
Logical. If TRUE asssume adult blood sample during cell count generation, else
assuming cord blood.
}


}

\details{
Input Files:

Sample_Sheet.csv:
.csv file containing phenotypic data for each sample
*See note for conventions

IDAT Files:
IDAT files containing methylation data for each sample. Should be contained in
the directory passed.

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
EpigeneticAgePipeline:::loadTestData()
directory <- paste0(path.package("EpigeneticAgePipeline"),"/data/")
generateResiduals(
    directory = directory,
    useBeta = TRUE,
    arrayType = "450K",
    ignoreCor = TRUE
)
message("Outliers")
knitr::kable(read.csv(paste0(directory, "OutlierSamples.csv")))
message("Residuals")
knitr::kable(read.csv(paste0(directory, "Residuals.csv")))
}



