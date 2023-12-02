\name{main}
\alias{main}

\title{main}

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
generateResiduals = TRUE,
useSampleSheet = TRUE)

}

\arguments{
\item{directory}{
Directory containing input data files (default: current working directory).
}

\item{normalize}{
Logical. Perform normalization of beta values if TRUE.
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

\item{generateResiduals}{
Whether or not to generate residuals.
(options to be specified during runtime)
}

\item{useSampleSheet}{
Whether or not to use Sample_Sheet.csv.
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

ChenEtAllList:
.csv file used for removing cross reactive probes. (450K array platform)

PidsleyCrossReactiveProbesEPIC.csv:
.csv file used for removing cross reactive probes. (EPIC array platform)

non-specific-probes-Illumina27k:
.csv file used for removing cross reactive probes. (27K array platform)

PC-clocks.rda and golden_ref.rda:
Prerequesite data for generating GrimAge age acceleration data.

}

\value{
The function will return the follwing:

output.txt: A .txt file containing
epigenetic age/acceleration estimates,
covariate data and residual data.

matrixplot{Clockname}.pdf: A set of .pdf files
illustrating a the correaltions between a specific
epigenetic age estimate and covariates.

epigeneticAge.txt: A .txt
file showing epigenetic age/acceleration estimates,
as well as some covariate data.
This file is better suited for importing
generated data into a spreadsheet program.

plot_{Clockname}.png: A set of
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

ID:
Variable containing the IDs of each sample.
Note: Although optional in most circumstances
adding an ID variable will significantly
reduce the chances of encountering a runtime error.

Array:
Array plate name in format "RXCY" (X: row number, Y: column number).

Age:
Age values of samples.

Sex:
Gender information of samples (M or F, 0 or 1)

Smoking_Status:
Smoking status of samples.

Batch:
Batch information of samples.

Slide:
Slide information of samples.

Bcell, CD4T, CD8T, Gran, Mono, nRBC:
Optional* cell type proportions if IDAT files are not available

*If data is not following the guidelines listed above,
all non-obliging variables will be treated as custom variables,
whose specification will be determined by the user during runtime.
}

\examples{
#Example of using the main function to generate epigenetic age values from
#IDAT files

message("The IDAT files used in this example were obtained from:\n")
message("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201872\n\n")

library(EpigeneticAgePipeline)
directory <- paste0(path.package("EpigeneticAgePipeline"),"/extdata/")

if (file.exists(paste0(directory, "betaValues.csv")))
{
    EpigeneticAgePipeline::main(
        directory = directory,
        normalize = TRUE,
        useBeta = FALSE,
        arrayType = "450K",
        generateResiduals = FALSE,
        useSampleSheet = FALSE
    )

    knitr::kable(exportDf)
}

}



