# Epigenetic Age Pipeline

'EpigeneticAgePipeline' is a comprehensive package designed for processing and analyzing DNA methylation data. The package provides a variety of epigenetic age measures, epigenetic age acceleration measures, residual generation, cell-count generation and provides a set of plots and tables used for further analysis.

## Installation and Use

1. Install the following packages which act as dependencies for the pipeline.

```
install.packages('remotes')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("minfi", quietly = TRUE))
    BiocManager::install("minfi")

if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE))
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse")

if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

if (!requireNamespace("ggpubr", quietly = TRUE))
    install.packages("ggpubr")

if (!requireNamespace("umap", quietly = TRUE))
    install.packages("umap")

 if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE))
    BiocManager::install("IlluminaHumanMethylationEPICmanifest")

if (!requireNamespace("methylclock", quietly = TRUE))
    BiocManager::install("methylclock")

if (!requireNamespace("DunedinPACE", quietly = TRUE))
    remotes::install_github("danbelsky/DunedinPACE")  # Use 'remotes' for devtools::install_github()

if (!requireNamespace("FlowSorted.CordBlood.450k", quietly = TRUE))
    BiocManager::install("FlowSorted.CordBlood.450k")

if (!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE))
    BiocManager::install("IlluminaHumanMethylation450kmanifest")

if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
    BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

if (!requireNamespace("lme4", quietly = TRUE))
    install.packages("lme4")

if (!requireNamespace("dnaMethyAge", quietly = TRUE))
    BiocManager::install("yiluyucheng/dnaMethyAge")
```

2. Use the 'remotes' package congruently with 'install_github', or load using library().

```
remotes::install_github('CastellaniLab/myEpigeneticAgePipeline')
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
