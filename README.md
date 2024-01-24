# Epigenetic Age Pipeline

'EpigeneticAgePipeline' is a comprehensive package designed for processing and analyzing DNA methylation data. The package provides a variety of epigenetic age measures, epigenetic age acceleration measures, residual generation, cell-count generation and provides a set of plots and tables used for further analysis.

## Installation and Use

1. Install the following packages which act as dependencies for the pipeline.

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages_to_install <- c(
  "minfi",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICmanifest",
  "methylclock",
  "methylclockData",
  "FlowSorted.CordBlood.450k",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation27kmanifest",
  "IlluminaHumanMethylation27kanno.ilmn12.hg19"
)

for (package in packages_to_install) {
    if (!requireNamespace(package, quietly = TRUE)) {
        BiocManager::install(package)
    }
}

if (!requireNamespace("S4Vectors", quietly = TRUE)) {
	BiocManager::install("S4Vectors", force = TRUE)
}

remotes::install_github('CastellaniLab/EpigeneticAgePipeline')
library(EpigeneticAgePipeline)
```

2. Run this code separately from the code above.

```
downloadDirectory <- "EpigeneticAgePipelineDataset-main"
downloadURL <- paste0("https://github.com/StanRaye/",
                        "EpigeneticAgePipelineDataset",
                        "/archive/refs/heads/main.zip")
installDirectory <- paste0(path.package("EpigeneticAgePipeline"),"/extdata/")
setwd(installDirectory)
utils::download.file(url = downloadURL, destfile = "asdf.zip")
utils::unzip("asdf.zip", exdir = installDirectory)
extractedFiles <- list.files(paste0(installDirectory, downloadDirectory),
                                    full.names = TRUE)
file.copy(from = extractedFiles, to = installDirectory, overwrite = TRUE,
            recursive = TRUE)
```

3. Once installed, the package is used by invoking the main function with the needed parameters.

```
library(EpigeneticAgePipeline)
main(
	directory = directory, #directory containing IDAT/beta values/supporting data
    normalize = TRUE, #should normalize beta values?
    useBeta = FALSE, #should use beta values contained within a betaValues.csv file?
    arrayType = "450K", #specification of array type used to geneerate idat files ("450K", "27K", "EPIC)
    generateResiduals = FALSE, #should generate residuals?
    useSampleSheet = FALSE #should use phenotypic data found within a Sample_Sheet.csv file?
)
```

## Description:

### Epigenetic Age and Acceleration Measures Provided:

**Horvath Clock:**

- Description: Provides DNAm age.
- CpG Sites: It uses a set of 353 CpG sites.
- Array Type: It was originally trained on Illumina 27K and 450K arrays.
- Cell Type: This clock is designed to work with a variety of tissues and cell types.

**Horvath2/skinHorvath:**

- Description: Provides DNAm age.
- CpG Sites: It uses a set of 391 CpG sites.
- Array Type: It was originally trained on Illumina 450K arrays.
- Cell Type: Primarly designed for skin/blood tissue.

**Hannum Clock:**

- Description: Provides DNAm age.
- CpG Sites: It relies on a set of 71 CpG sites.
- Array Type: This clock was trained on Illumina 450K arrays.
- Cell Type: Designed primarly for blood samples.

**Levine/PhenoAge Clock:**

- Description: Provides DNAm age.
- CpG Sites: Uses 513 CpG sites.
- Array Type: This clock was trained on Illumina 27K, 450K and EPIC array types.
- Cell Type: Designed primarly for blood samples.

**DunedinPACE:**

- Description: Provides the estimated rate of biological aging.
- CpG Sites: Uses 173 CpG sites.
- Array Type: Trained using Illumina 450K and EPIC array types.
- Cell Type: Designed primarly for blood samples.

**GrimAge:**

- Description: Provides a measure of epigenetic age acceleration.
- CpG Sites: Uses 1030 CpG sites.
- Array Type: Trained using Illumina 450K and EPIC array types.
- Cell Type: Designed primarly for blood samples.

### Cell Counts

If IDAT files are provided, methylation data can be used to determine cell counts of the following cell types.

- B Cells
- CD4T Cells
- CD8T Cells
- Granulocytes
- Monocytes
- Nucleated Red Blood cells

### Introduction to Residual Generation

This section describes the proccess of residual generation. The function constructs a formula string for the linear model based on the presence of the variables “Row”, “Column”, “Slide”, and “Batch”.

#### Dynamic Formula Construction for Linear Models
Based on the available data variables, the function dynamically constructs the formula for the linear model. Below are the possible formulae:

| Variables Present             | Formula                                                   |
|-------------------------------|-----------------------------------------------------------|
| Without “Column”              | `EpigeneticAgeMeasure ~ Xi + (Row\|Slide) + (1\|Batch)`   |
| Without “Slide”               | `EpigeneticAgeMeasure ~ Xi + (Row & Column) + (1\|Batch)` |
| Both “Column” and “Slide”     | `EpigeneticAgeMeasure ~ Xi + (1\|Slide) + (Row + Column\|Slide) + (1\|Batch)` |
| Without “Row” or “Batch”      | `EpigeneticAgeMeasure ~ Xi`                               |

*Note: In these formulae, `Xi` represents the independent variables.*


#### Implementation of the Linear Model Generation

Once the formula is constructed, the linear model is generated using the Gaussian method. The function glmmTMB from the glmmTMB package is used to fit the model. The maximum number of iterations and evaluations for the optimizer are set to 10000 to ensure convergence. The user can specify to remove highly correlated explanatory variables (greater than 0.6) during runtime.

### Output

**output.txt:**  
A .txt file containing epigenetic age/acceleration estimates, covariate data and residual data.

**matrixplot{Clockname}.pdf:**  
A set of .pdf files illustrating a the correaltions between a specific epigenetic age estimate and covariates.

**epigeneticAge.txt:**  
A .txt file showing epigenetic age/acceleration estimates. This file is better suited for importing into a spreadsheet program than output.txt.

**plot{Clockname}.pdf:**  
A set of .png files showing a line plot of an epigenetic age estimate against chronological age.

**SampleIDandAge.png:**  
A .png file containing a grouped bar chart showing each sample and their associated epigenetic age estimates as well as chronological age. Note that this file is typically a more useful analysis tool when using 
smaller sample sizes.

### Common Erros
1. If using the package Compute Canada, and are getting a errors such as:
```
Error in ExperimentHub::ExperimentHub() :
  DEFUNCT: As of ExperimentHub (>1.17.2), default caching location has changed.
  Problematic cache: /users/ccastell/.cache/ExperimentHub
  See https://bioconductor.org/packages/devel/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#default-caching-location-update
```
Try running:
```
library(ExperimentHub)
oldcache = path.expand(rappdirs::user_cache_dir(appname="ExperimentHub"))
setExperimentHubOption("CACHE", oldcache)
eh = ExperimentHub(localHub=TRUE)

## removes old location and all resources
removeCache(eh, ask=FALSE)

## create the new default caching location
newcache = tools::R_user_dir("ExperimentHub", which="cache")
setExperimentHubOption("CACHE", newcache)
eh = ExperimentHub()
```



### References

Citation for dnaMethyAge Package:
Wang, Y., Grant, O. A., Zhai, X., McDonald-Maier, K. D., & Schalkwyk, L. C. (2023). Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock. *GeroScience*, 1–18. [Link to article](https://link.springer.com/article/10.1007/s11357-023-00871-w)
