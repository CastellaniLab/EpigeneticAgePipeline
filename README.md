# Epigenetic Age Pipeline

'EpigeneticAgePipeline' is a comprehensive package designed for processing and analyzing DNA methylation data. The package provides a variety of epigenetic age measures, epigenetic age acceleration measures, residual generation, cell-count generation and provides a set of plots and tables used for further analysis.

## Installation 

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
  "IlluminaHumanMethylation27kanno.ilmn12.hg19",
  "AnnotationHub",
  "AnnotationHubData",
  "base64",
  "beanplot",
  "ggplot2",
  "dplyr",
  "tidyr",
  "annotate",
  "magick",
  "pdftools",
  "S4Vectors",
  "preprocessCore",
  "reshape2",
  "glmmTMB",
  "SparseM",
  "nloptr"
)
for (package in packages_to_install) {
    if (!requireNamespace(package, quietly = TRUE)) {
        BiocManager::install(package)
    }
}
remotes::install_github('CastellaniLab/EpigeneticAgePipeline')
```

2. Once installed, the package is used by invoking the main function with the needed parameters.

```
library(EpigeneticAgePipeline)
main(directory = directory, #directory containing IDAT/beta values/supporting data
    normalize = TRUE, #should normalize beta values?
    useBeta = FALSE, #should use beta values contained within a betaValues.csv file?
    arrayType = "450K", #specification of array type used to geneerate idat files ("450K", "27K", "EPIC)
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

This section describes the proccess of residual generation. The function constructs a formula string for the linear model based on the presence of the variables “Row”, “Column”, “Slide”, and “Batch” in Sample_Sheet.csv. If you have access to these variables in your Sample_Sheet.csv, make sure to append "@@@2" to each variable, and mind capitalization. More info on using the "Array" variable is in Usage Guidlines.    

#### Dynamic Formula Construction for Linear Models
Based on the available data variables, the function dynamically constructs the formula for the linear model. Below are the possible formulae:

```mermaid
graph TD
    A["Start"]
    B{"Available Variables"}
    C["formula_string = 'EpiAge ~ Xi + (Row|Slide) + (1|Batch)'"]
    D["formula_string = 'EpiAge ~ Xi + (Row&Column) + (1|Batch)'"]
    E["formula_string = 'EpiAge ~ Xi + (1|Slide) + (Row+Column|Slide) + (1|Batch)'"]
    F["formula_string = 'EpiAge ~ Xi'"]
    G["End"]

    A --> B
    B -->|Not 'Column', Yes 'Row', 'Slide', 'Batch'| C
    B -->|Not 'Slide', Yes 'Row', 'Column', 'Batch'| D
    B -->|Yes 'Row', 'Column', 'Slide', 'Batch'| E
    B -->|Else| F
    C --> G
    D --> G
    E --> G
    F --> G
```

*Note: In these formulae, `Xi` represents the independent variables.*


#### Implementation of the Linear Model Generation

Once the formula is constructed, the linear model is generated using the Gaussian method. The function glmmTMB from the glmmTMB package is used to fit the model. The maximum number of iterations and evaluations for the optimizer are set to 10000 to ensure convergence. The user can specify to remove highly correlated explanatory variables (greater than 0.6) during runtime.

## Usage Guidelines 
### Using the main Function
```
main(directory = getwd(),
normalize = TRUE,
useBeta = FALSE,
arrayType = "450K",
useSampleSheet = TRUE)
```
**directory** argument:  
Directory containing input data files (default: current working directory).

**normalize** argument:  
Logical. Perform normalization of beta values if TRUE.

**useBeta** argument:  
Logical. If TRUE, will expect a betaValues.csv file containing beta values (scaled between 0 and 1). If FALSE, process raw intensity data (IDAT).

**arrayType** argument:  
Type of DNA methylation array used (options: "27K", "450K", or "EPIC").

**useSampleSheet** argument:  
Logical. If TRUE, will expect a Sample_Sheet.csv containing phenotypic data.

#### Output  

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

### Using the generateResiduals function
```
generateResiduals(directory = getwd(),
useBeta = FALSE,
arrayType = "450K"
``` 
**directory** argument:  
Directory containing input data files (default: current working directory). 

**useBeta** argument:  
Logical. If TRUE, will expect a betaValues.csv file containing beta values (scaled between 0 and 1). If FALSE, process raw intensity data (IDAT).  

**arrayType** argument:  
Type of DNA methylation array used (options: "27K", "450K", or "EPIC").  

#### Output  

**Residuals.csv:**  
A .csv file containing residuals from the linear model.

**ResidualsAcceleartion.csv:**  
A .csv file containing age acceleration residuals from the linear model.

### Description of Client-Side Input Files
**Sample_Sheet.csv**  
.csv file containing phenotypic data for each sample.
*Guidlines listed below

**IDAT Files**  
IDAT files containing methylation data for each sample.

**betaValues.csv**  
If IDAT files are not available, processed beta values can be provided. First **column** should contain CpG names. First **row** should contain sample names. 

### Guidelines for Sample_Sheet.csv
**To include a variable from Sample_Sheet.csv:**  
If a variable contains **numeric** type data, append "@@@1" to the column name  

If a variable contains **factor** type data, append "@@@2" to the column name.  

Ex. variable "isSmoker" would become "isSmoker@@@2".  

To generate **GrimAge**, chronological age must be included ("Age@@@1"), as well as sample sex ("Sex@@@2"). Valid values for male sex would be "1", "M" or "Male". Valid values for female sex would be "2", "F", or "Female".  

If using the **generateResiduals** function, name the column with epigenetic age values as "EpiAge@@@1".  

If you want to use random effects in the **generateResiduals** function, view the guidelines found in "Introduction to Residual Generation" above.  

**Specification for Using 'Array' Variable**  
If using an Array variable to store row and column information, please make sure it follows the format "RXCY" (X: row number, Y: column number).  
Make sure to append @@@2 to the column name  


### Common Erros
1. If using the package Compute Canada, and are getting errors such as:
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
