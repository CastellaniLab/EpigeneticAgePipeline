\name{generateResiduals}
\alias{generateResiduals}

\title{Generate Epigenetic Age Residual Data}

\description{
Performs principal component analysis (PCA) and residual generation
after epigenetic age measures are available. Optionally includes
cell-type proportions, derives array layout covariates, and exports
diagnostic plots and tables.
}

\usage{
generateResiduals(
  inputDirectory = getwd(),
  outputDirectory = inputDirectory,
  sampleSheetFile = "Sample_Sheet.csv",
  columnTypes = NULL,
  useBeta = FALSE,
  formula = NULL,
  arrayType = "450K",
  ignoreCor = TRUE,
  PCs = 5,
  threshold = 3,
  doParallel = TRUE,
  doCellCounts = TRUE,
  tissue = "bloodAdult",
  cellDeconvMethod = "CP",
  placentaTrimester = "third"
)
}

\arguments{
\item{inputDirectory}{
String. Directory containing input data files (default: current working directory).
}

\item{outputDirectory}{
String. Directory to write outputs (default: \code{inputDirectory}).
}

\item{sampleSheetFile}{
String. File name of the sample sheet (default: \code{"Sample_Sheet.csv"}).
}

\item{columnTypes}{
Named integer vector or \code{NULL}. Declares sample-sheet column types:
\code{1 = factor}, \code{2 = numeric}.\\
Example: \code{c("EpiAge"=2,"Age"=2,"Sex"=1,"Batch"=1,"Array"=1)}.\\
If \code{Array} is in \code{RXCY} format (e.g., \code{R01C05}), \code{Row} and \code{Column} are derived automatically.
\emph{Required} when reading a sample sheet.
}

\item{useBeta}{
Logical. If \code{TRUE}, expect \code{betaValues.csv} (or \code{.csv.gz}) with beta values
(first column CpG IDs, remaining columns samples). If \code{FALSE}, process raw IDATs.
}

\item{formula}{
String or \code{NULL}. Custom model formula (e.g., \code{"EpiAge ~ Age + BMI + Smoking_Status + PC1"}).\\
PCs are referenced as \code{PC1}, \code{PC2}, ...\\
Adult blood cells: \code{Bcell, CD4T, CD8T, Mono, Neu, NK}.\\
Cord blood cells: \code{Bcell, CD4T, CD8T, Gran, Mono, nRBC, NK}.\\
Saliva cells: \code{Epithelial, Immune}.\\
Placenta cells: \code{Trophoblasts, Stromal, Hofbauer, Endothelial, nRBC, Syncytiotrophoblast}.\\
Frontal cortex cells: \code{NeuN_neg, NeuN_pos}.\\
Buccal cells: \code{NeuN_neg, NeuN_pos}.\\
If \code{NULL}, a formula is constructed dynamically from available covariates.
}

\item{arrayType}{
String. DNA methylation array platform: \code{"27K"}, \code{"450K"}, \code{"EPIC"}, \code{"EPICv2"}, or \code{"MSA"}.
}

\item{ignoreCor}{
Logical. If \code{TRUE}, skip pruning of highly correlated covariates (|r| > 0.6); if \code{FALSE}, an interactive removal step is offered.
}

\item{PCs}{
Integer. Number of principal components to compute/include (default: \code{5}).
}

\item{threshold}{
Numeric. MAD-based outlier threshold; flags samples outside \code{median ± threshold × MAD} (default: \code{3}).
}

\item{doParallel}{
Logical. If \code{TRUE}, use parallelized CSV reading for large \code{betaValues.csv}.
}

\item{doCellCounts}{
Logical. If \code{TRUE}, estimate and include cell-type proportions in the model (requires IDATs).
}

\item{tissue}{
String. Reference panel for cell deconvolution: \code{"bloodAdult"}, \code{"bloodCord"}, \code{"saliva"}, \code{"placenta"}, \code{"brainDLPFC"} or \code{"buccal"}.
}

\item{cellDeconvMethod}{
String. Cell proportion method: \code{"CP"} (constrained projection) or \code{"RPC"} (robust partial correlations; \pkg{EpiDISH}).
}

\item{placentaTrimester}{
String. When \code{tissue = "placenta"}, choose \code{"first"} or \code{"third"} trimester reference.
}
}

\details{
\strong{Input Files}

\emph{Sample_Sheet.csv:} One row per sample, one column per variable.
Variable types are declared via \code{columnTypes}. If an \code{Array} column in \code{RXCY}
format is present, \code{Row} and \code{Column} are derived automatically for mixed-effects terms.

\emph{IDAT files:} Raw intensity files for each sample in \code{inputDirectory}.

\emph{betaValues.csv / betaValues.csv.gz:} Optional precomputed beta matrix; first column CpG IDs, remaining columns samples.
}

\value{
Writes results to \code{outputDirectory} and assigns selected objects to \code{.GlobalEnv}:

\itemize{
  \item \code{Residuals.csv}: Residuals from the fitted (G)LMM regressing the chosen age measure on covariates/PCs/cell counts.
  \item \code{OutlierSamples.csv}: Samples flagged by MAD-based outlier detection across selected PCs.
  \item \code{\{AgeType\}SampleData.csv}: Covariate frame used in modeling (e.g., \code{EpiAgeSampleData.csv}).
  \item \code{matrixplot\{AgeType\}.png}: Correlation matrix of the target age measure versus covariates (e.g., \code{matrixplotEpiAge.png}).
  \item \code{CellCountsDf} (R object, if applicable): Estimated cell-type proportions assigned to \code{.GlobalEnv}.
  \item \code{myEnv} (R environment): Stores intermediate/final objects (\code{bVals}, \code{pdataSVs}, \code{outliersCSV}, \code{residualsCSV}, ...).
}
}

\note{
\strong{Sample sheet typing via \code{columnTypes}} — declare variable types with a named integer vector: \code{1 = factor}, \code{2 = numeric}.
For GrimAge/PCGrimAge, include \code{Age} (numeric) and \code{Sex} (factor); valid sex encodings include \code{1}/\code{Male} and \code{2}/\code{Female} (case-insensitive).
If an \code{Array} column is supplied in \code{RXCY} format (e.g., \code{R01C05}), \code{Row} and \code{Column} are auto-derived for random-effects construction.
}

\examples{
# Example of using generateResiduals
library(EpigeneticAgePipeline)

directory <- paste0(path.package("EpigeneticAgePipeline"), "/extdata/")

ct <- c(
  "EpiAge" = 2,
  "Age"    = 2,
  "Sex"    = 1
)

generateResiduals(
  inputDirectory   = directory,
  outputDirectory  = getwd(),
  sampleSheetFile  = "Sample_Sheet.csv",
  columnTypes      = ct,
  useBeta          = TRUE,
  formula          = NULL,
  arrayType        = "450K",
  ignoreCor        = TRUE,
  PCs              = 5,
  threshold        = 3,
  doParallel       = TRUE,
  doCellCounts     = TRUE,
  tissue           = "bloodAdult",
  cellDeconvMethod = "CP",
  placentaTrimester= "third"
)

message("Outliers")
knitr::kable(read.csv("OutlierSamples.csv"))
message("Residuals")
knitr::kable(read.csv("Residuals.csv"))
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


