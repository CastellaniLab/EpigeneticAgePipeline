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
main(inputDirectory = getwd(),
outputDirectory = inputDirectory,
normalize = TRUE,
useBeta = FALSE,
arrayType = "450K",
useSampleSheet = TRUE,
doParallel = TRUE,
writeBeta = TRUE,
useAdult = FALSE,
useImputation = FALSE)

}

\arguments{
\item{inputDirectory}{
String. Directory containing input data files (default: current working directory).
}

\item{outputDirectory}{
String. Directory to output created files (default: inputDirectory).
}

\item{normalize}{
Logical. Perform normalization of beta values if TRUE. (Horvath's normalization, from methylclock)
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

\item{useImputation}{
Logical. If TRUE use mean values from GSE40279 to fill in missing CpG's during
GrimAge and DunedinPACE calculation.
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

betaValues.csv or betaValues.csv.gz:
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
as well as sample sex ("Sex@@@2"). Valid values for male sex would be "1", or "Male".
Valid values for female sex would be "2", or "Female". Mind capatalization.

For a more detailed guide, visit the github page:
https://github.com/CastellaniLab/EpigeneticAgePipeline

Citations:


}

\examples{
#Example of using the main function
library(EpigeneticAgePipeline)
directory <- paste0(path.package("EpigeneticAgePipeline"),"/extdata/")
main(
    inputDirectory = directory,
    outputDirectory = getwd(),
    normalize = TRUE,
    useBeta = TRUE,
    arrayType = "450K",
    useSampleSheet = FALSE
)
knitr::kable(read.table("epigeneticAge.txt"))
}

\references{
Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA (2014). “Minfi: A flexible and
  comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays.” _Bioinformatics_,
  *30*(10), 1363-1369. doi:10.1093/bioinformatics/btu049 <https://doi.org/10.1093/bioinformatics/btu049>.

Triche TJ, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD (2013). “Low-level processing of Illumina Infinium DNA
  Methylation BeadArrays.” _Nucleic Acids Research_, *41*(7), e90. doi:10.1093/nar/gkt090
  <https://doi.org/10.1093/nar/gkt090>.

Fortin J, Triche TJ, Hansen KD (2017). “Preprocessing, normalization and integration of the Illumina
  HumanMethylationEPIC array with minfi.” _Bioinformatics_, *33*(4). doi:10.1093/bioinformatics/btw691
  <https://doi.org/10.1093/bioinformatics/btw691>.

Dolors Pelegri-Siso, Paula de Prado, Justiina Ronkainen, Mariona Bustamante, Juan R Gonzalez, methylclock: a
  Bioconductor package to estimate DNA methylation age, Bioinformatics, Volume 37, Issue 12, 15 June 2021, Pages
  1759–1760, https://doi.org/10.1093/bioinformatics/btaa825.

Wang Y (2025). _dnaMethyAge: Predict epigenetic age from DNA methylation data and calculate age acceleration._. R
  package version 0.2.0, commit f8ef2c7a781193cef426a28de4d1c6eed4bc7654, <https://github.com/yiluyucheng/dnaMethyAge>.

Fortin J, Hansen KD (2016).
  _IlluminaHumanMethylation27kanno.ilmn12.hg19:
  Annotation for Illumina's 27k methylation arrays_. R
  package version 0.6.0.

Hansen KD (2021).
  _IlluminaHumanMethylation450kanno.ilmn12.hg19:
  Annotation for Illumina's 450k methylation arrays_. R
  package version 0.6.1.

Hansen KD (2017).
  _IlluminaHumanMethylationEPICanno.ilm10b4.hg19:
  Annotation for Illumina's EPIC methylation arrays_. R
  package version 0.6.0,
  <https://bitbucket.com/kasperdanielhansen/Illumina_EPIC>.

Gu Z (2024).
  _IlluminaHumanMethylationEPICv2anno.20a1.hg38:
  Annotation for Illumina's EPIC v2.0 methylation
  arrays_. R package version 1.0.0,
  <https://www.illumina.com/products/by-type/microarray-kits   /infinium-methylation-epic.html>.

MacDonald JW (2024).
  _IlluminaHumanMethylationMSAanno.ilm10a1.hg38:
  Annotation for Illumina's MSA methylation arrays_. R
  package version 0.1.0,
  <https://github.com/jmacdon/IlluminaHumanMethylationMSAanno.ilm10a1.hg38>.

Salas LA, Gervin K, Jones MC (2024).
  _FlowSorted.CordBloodCombined.450k: Illumina 450k/EPIC
  data on FACS and MACS umbilical blood cells_.
  doi:10.18129/B9.bioc.FlowSorted.CordBloodCombined.450k
  <https://doi.org/10.18129/B9.bioc.FlowSorted.CordBloodCombined.450k>,
  R package version 1.22.0,
  <https://bioconductor.org/packages/FlowSorted.CordBloodCombined.450k>.

Salas LA, Koestler DC (2024). _FlowSorted.Blood.EPIC:
  Illumina EPIC data on immunomagnetic sorted peripheral
  adult blood cells_.
  doi:10.18129/B9.bioc.FlowSorted.Blood.EPIC
  <https://doi.org/10.18129/B9.bioc.FlowSorted.Blood.EPIC>,
  R package version 2.10.0,
  <https://bioconductor.org/packages/FlowSorted.Blood.EPIC>.
}

