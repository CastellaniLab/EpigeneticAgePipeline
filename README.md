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

3. Once installed, you can use the package in the following ways:

```
library(EpigeneticAgePipeline)
main(...)
```

```
EpigeneticAgePipeline::main(...)
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

### Residual Generation Linear Model

This section describes the proccess of residual generation. The function constructs a formula string for the linear model based on the presence of the variables “Row”, “Column”, “Slide”, and “Batch”.

#### Possible Combinations of Random Effects

The function checks for the presence of these variables in the data and constructs the formula accordingly. Here are the possible combinations of random effects:

If “Column” is not present, the random effects include “Row” nested within “Slide” and “Batch” as a random intercept.

_EpigeneticAgeMeasure ~ X<sub>i</sub> + (Row|Slide) + (1|Batch)_

If “Slide” is not present, the random effects include an interaction of “Row” and “Column” and “Batch” as a random intercept.

_EpigeneticAgeMeasure ~ X<sub>i</sub> + (Row & Column) + (1|Batch)_

If both “Column” and “Slide” are present, the random effects include “Slide” as a random intercept, an interaction of “Row” and “Column” nested within “Slide”, and “Batch” as a random intercept.

_EpigeneticAgeMeasure ~ X<sub>i</sub> + (1|Slide) + (Row + Column|Slide) + (1|Batch)_

If “Row” or “Batch” is not present, no random effects are included.

_EpigeneticAgeMeasure ~ X<sub>i</sub>_

#### Linear Model Generation

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
A .png file containing a grouped bar chart showing each sample and their assocaited epigenetic age estimates as well as chronological age. Note that this file is typically a more useful analysis tool when using smaller sample sizes.

### References

Citation for dnaMethyAge Package:
Wang, Y., Grant, O. A., Zhai, X., McDonald-Maier, K. D., & Schalkwyk, L. C. (2023). Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock. *GeroScience*, 1–18. [Link to article](https://link.springer.com/article/10.1007/s11357-023-00871-w)
