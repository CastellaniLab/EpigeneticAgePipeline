\name{main}
\alias{main}

\title{Epigenetic Age Generation}

\description{
A comprehensive function designed for
processing and analyzing DNA methylation data.
Its goal is to provide a simple solution
for generating epigenetic age and performing covariate analysis.
}

\usage{
main(directory = getwd(),
normalize = TRUE,
useBeta = FALSE,
arrayType = "450K",
useSampleSheet = TRUE,
doParallel = TRUE,
writeBeta = TRUE,
useAdult = FALSE)

}

\arguments{
\item{directory}{
Directory containing input data files (default: current working directory).
}

\item{normalize}{
Logical. Perform normalization of beta values if TRUE. (Horvath's normalization)
}

\item{useBeta}{
Logical. If TRUE,
expect a .csv containing beta values
(scaled between 0 and 1).
If FALSE, process raw intensity data (IDAT).
}

\item{arrayType}{
Type of DNA methylation array used (options: "27K", "450K", "EPIC", "EPICv2"
or "MSA").
}

\item{useSampleSheet}{
Logical. If TRUE use Sample_Sheet.csv.
}

\item{doParallel}{
Logical. If TRUE use parallel processing in .csv file reads.
}

\item{writeBeta }{
Logical. If TRUE write extracted beta values from processed IDAT files.
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
If IDAT files are not available, a .csv containing beta values (CpG names first
column, then sample columns) can be provided.

}

\value{
The function will return the follwing:

output.txt: A .txt file containing
epigenetic age/acceleration estimates,
covariate data and residual data.

matrixplot[Clockname].pdf: A set of .pdf files
illustrating a the correaltions between a specific
epigenetic age estimate and covariates.

epigeneticAge.txt: A .txt
file showing epigenetic age/acceleration estimates,
as well as some covariate data.
This file is better suited for importing
generated data into a spreadsheet program.

plot_[Clockname].png: A set of
.png files showing a line
plot of an epigenetic age
estimate against chronological age.

SampleIDandAge.png:
A .png file containing a grouped bar chart showing each
sample and their assocaited epigenetic age
estimates as well as chronological age.
Note that this file is typically a more useful
analysis tool when using smaller sample sizes.

}


\note{
Naming Conventions for Data in "Sample_Sheet.csv":

To include a variable from Sample_Sheet.csv:
If a variable contains numeric type data, append "@@@1" to the column name

If a variable contains factor type data, append "@@@2" to the column name.

Ex. variable "isSmoker" would become "isSmoker@@@2".

To generate GrimAge, chronological age must be included ("Age@@@1"),
as well as sample sex ("Sex@@@2"). Valid values for male sex would be "1",
"M" or "Male". Valid values for female sex would be "2", "F", or "Female".
Mind capatalization.

For a more detailed guide, visit the github page:
https://github.com/CastellaniLab/EpigeneticAgePipeline

}

\examples{
#Example of using the main function
library(EpigeneticAgePipeline)
EpigeneticAgePipeline:::loadTestData()
directory <- paste0(path.package("EpigeneticAgePipeline"),"/data/")
main(
    directory = directory,
    normalize = TRUE,
    useBeta = TRUE,
    arrayType = "450K",
    useSampleSheet = FALSE
)
knitr::kable(read.table(paste0(directory, "epigeneticAge.txt")))
EpigeneticAgePipeline:::removeTestData()
}



