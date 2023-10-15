# Epigenetic Age Pipeline

'EpigeneticAgePipeline' is a comprehensive package designed for processing and analyzing DNA methylation data. The package provides a variety of epigenetic age measures, epigenetic age acceleration measures, residual generation, cell-count generation and provides a set of plots and tables used for further analysis.

## Installation and Use

1. Install 'remotes' package in order to install R pacakges hosted on remote repositories, as well as Bioconducter which will support some associated packages.

```
install.packages('remotes')
install.packages("BiocManager")
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

**_ Horvath Clock: _**

- Description: Provides DNAm age.
- CpG Sites: It uses a set of 353 CpG sites.
- Array Type: It was originally trained on Illumina 27K and 450K arrays.
- Cell Type: This clock is designed to work with a variety of tissues and cell types.

**_ Horvath2/skinHorvath: _**

- Description: Provides DNAm age.
- CpG Sites: It uses a set of 391 CpG sites.
- Array Type: It was originally trained on Illumina 450K arrays.
- Cell Type: Primarly designed for skin/blood tissue.

**_ Hannum Clock: _**

- Description: Provides DNAm age.
- CpG Sites: It relies on a set of 71 CpG sites.
- Array Type: This clock was trained on Illumina 450K arrays.
- Cell Type: Designed primarly for blood samples.

**_ Levine/PhenoAge Clock: _**

- Description: Provides DNAm age.
- CpG Sites: Uses 513 CpG sites.
- Array Type: This clock was trained on Illumina 27K, 450K and EPIC array types.
- Cell Type: Designed primarly for blood samples.

**_ DunedinPACE: _**

- Description: Provides the estimated rate of biological aging.
- CpG Sites: Uses 173 CpG sites.
- Array Type: Trained using Illumina 450K and EPIC array types.
- Cell Type: Designed primarly for blood samples.

**_ GrimAge: _**

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
