\name{main}
\alias{main}

\title{Epigenetic Age Generation}

\description{
A comprehensive function for processing and analyzing DNA methylation data.
It generates epigenetic age measures, optional accelerations, cell-type
estimates, and covariate analyses with plotting and export utilities.
}

\usage{
main(
 inputDirectory = getwd(),
 outputDirectory = inputDirectory,
 normalize = TRUE,
 useBeta = FALSE,
 arrayType = "450K",
 useSampleSheet = TRUE,
 sampleSheetFile = "Sample_Sheet.csv",
 columnTypes = NULL,
 doParallel = TRUE,
 writeBeta = TRUE,
 tissue = "bloodAdult",
 cellDeconvMethod = "CP",
 placentaTrimester = "third",
 useImputation = FALSE,
 detPSampleCutoff = 0.05,
 detPProbeCutoff = 0.01,
 minBeads = 3,
 beadSampleMaxFailFrac = 0.05,
 beadProbeMaxFailFrac  = 0.05
)
}

\arguments{
\item{inputDirectory}{
String. Directory containing input data files (default: current working directory).
}

\item{outputDirectory}{
String. Directory to write outputs (default: \code{inputDirectory}).
}

\item{normalize}{
Logical. If \code{TRUE}, apply Horvath-style normalization (as implemented in \pkg{methylclock}) when computing CpG-based clocks.
}

\item{useBeta}{
Logical. If \code{TRUE}, expect a \code{betaValues.csv} (or \code{.csv.gz}) with beta values
(first column CpG IDs, remaining columns sample IDs). If \code{FALSE}, process raw IDATs.
}

\item{arrayType}{
String. DNA methylation array platform: \code{"27K"}, \code{"450K"}, \code{"EPIC"}, \code{"EPICv2"}, or \code{"MSA"}.
}

\item{useSampleSheet}{
Logical. If \code{TRUE}, read phenotypic/covariate data from \code{sampleSheetFile}.
}

\item{sampleSheetFile}{
String. File name of the sample sheet (default: \code{"Sample_Sheet.csv"}).
}

\item{columnTypes}{
Named integer vector or \code{NULL}. Declares types for columns in the sample sheet: \code{1 = factor}, \code{2 = numeric}.\\
Example: \code{c("Age"=2, "Sex"=1, "Batch"=1, "BMI"=2, "Array"=1)}.\\
If an \code{Array} column is present in \code{RXCY} format (e.g., \code{R01C05}), \code{Row} and \code{Column} are derived automatically.
When \code{useSampleSheet = TRUE}, this argument is required.
}

\item{doParallel}{
Logical. If \code{TRUE}, use parallelized CSV reading for large \code{betaValues.csv}.
}

\item{writeBeta}{
Logical. If \code{TRUE}, write \code{extractedBetaValues.csv} after IDAT processing (ignored when \code{useBeta = TRUE}).
}

\item{tissue}{
String. Reference panel for cell deconvolution: \code{"bloodAdult"}, \code{"bloodCord"}, \code{"saliva"}, \code{"placenta"}, \code{"brainDLPFC"} or \code{"buccal"}.
}

\item{cellDeconvMethod}{
String. Cell proportion method: \code{"CP"} (constrained projection) or \code{"RPC"} (robust partial correlations; \pkg{EpiDISH}).
}

\item{placentaTrimester}{
String. When \code{tissue = "placenta"}, choose \code{"first"} or \code{"third"} trimester reference (via \pkg{planet}).
}

\item{useImputation}{
Logical. If \code{TRUE}, perform mean-imputation of missing CpGs (from GSE40279) for PC clocks and DunedinPACE using packaged reference means.
}

\item{detPSampleCutoff}{
Numeric between 0 and 1. Sample-level detection-p filter.
For each sample, the mean detection-p across probes is computed.
Samples with mean detection-p greater than or equal to this cutoff are removed.
Default: 0.05.
}

\item{detPProbeCutoff}{
Numeric between 0 and 1. Probe-level detection-p filter.
After sample filtering, a probe is retained only if its detection-p is less than this cutoff
in all remaining samples. Probes failing in any sample are dropped.
Default: 0.01.
}

\item{minBeads}{
Integer greater than or equal to 0. Minimum bead count required for a probe–sample
measurement to be considered reliable. Observations with bead count less than this
value are marked as low-bead failures for downstream fraction tests.
Default: 3.
}

\item{beadSampleMaxFailFrac}{
Numeric between 0 and 1. Maximum fraction of low-bead observations permitted per sample.
Samples with a proportion of low-bead probe calls greater than this value are removed.
Default: 0.05.
}

\item{beadProbeMaxFailFrac}{
Numeric between 0 and 1. Maximum fraction of low-bead observations permitted per probe
across the (post-filtered) samples. Probes with a proportion of low-bead calls greater
than this value are removed. Default: 0.05.
}
}

\details{
\strong{Input Files}

\emph{Sample_Sheet.csv:} A CSV with one row per sample and one column per variable.
Types are declared via \code{columnTypes} (see \samp{columnTypes} argument). If an \code{Array} column
is supplied in \code{RXCY} format, \code{Row} and \code{Column} are derived automatically.

\emph{IDAT files:} Raw intensity files (RG) for each sample in \code{inputDirectory}.

\emph{betaValues.csv / betaValues.csv.gz:} Optional precomputed beta matrix; first column CpG IDs, remaining columns samples.
When both exist, the \code{.gz} version is preferred.
}

\value{
The function writes files to \code{outputDirectory} and assigns selected objects to \code{.GlobalEnv}:

\itemize{
  \item \code{output.txt}: Text summary of covariate--clock correlations across processed age measures.
  \item \code{epigeneticAge.txt}: Numeric table of epigenetic age estimates per sample (spreadsheet-friendly).
  \item \code{results.md}: Markdown table summarizing epigenetic age results.
  \item \code{extractedBetaValues.csv} (conditional): Written if \code{useBeta = FALSE} and \code{writeBeta = TRUE}.
  \item \code{GroupedAge.png}: Grouped bar chart of age measures per sample (includes chronological age when available).
  \item \code{matrixplot\{ClockName\}.png}: Correlation matrix per clock versus covariates.
  \item \code{plot_\{ClockName\}.png}: Scatter plots of DNAm age vs. chronological age (per clock when \code{Age} is present).
  \item \code{\{ClockName\}SampleData.csv}: Per-clock covariate frames used in plotting/correlation.
  \item \code{EpiAgeResultsDf} (R object): Data frame of numeric results assigned to \code{.GlobalEnv}.
  \item \code{CellCountsDf} (R object, if applicable): Estimated cell-type proportions assigned to \code{.GlobalEnv}.
  \item \code{myEnv} (R environment): Working environment storing intermediate/final objects (\code{bVals}, \code{pdataSVs}, etc.).
}
}

\note{
\strong{Sample sheet typing via \code{columnTypes}}

Declare variable types with a named integer vector: \code{1 = factor}, \code{2 = numeric}.
Required for \code{useSampleSheet = TRUE}.
\emph{GrimAge requirements:} Include \code{Age} (numeric) and \code{Sex} (factor). Valid sex encodings include \code{1}/\code{Male} and \code{2}/\code{Female} (case-insensitive).
If an \code{Array} column is supplied in \code{RXCY} format, \code{Row} and \code{Column} are auto-derived for mixed-effects modeling.

For additional guidance, see the package repository.
}

\examples{
# Example of using the main function
library(EpigeneticAgePipeline)

directory <- paste0(path.package("EpigeneticAgePipeline"), "/extdata/")

main(
  inputDirectory   = directory,
  outputDirectory  = getwd(),
  normalize        = TRUE,
  useBeta          = TRUE,
  arrayType        = "450K",
  useSampleSheet   = FALSE,
  sampleSheetFile  = "Sample_Sheet.csv",
  columnTypes      = ct,
  doParallel       = TRUE,
  writeBeta        = TRUE,
  tissue           = "bloodAdult",
  cellDeconvMethod = "CP",
  placentaTrimester= "third",
  useImputation    = FALSE
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

