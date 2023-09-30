# Epigenetic Age Pipeline 

## Installation 

1. Install "remotes" package in order to install R pacakges hosted on remote repositories.
```
install.packages("remotes")
```
2. Use the "remotes" package congruently with "install_github".
```
remotes::install_github("StanRaye/myEpigeneticAgePipeline")
```

## Function

### main(normalize = TRUE, useBeta = FALSE, directory = getwd(), arrayType = "450K")

## Description

The DNAm-Age Pipeline is a comprehensive function designed for processing and analyzing DNA methylation data. It performs data normalization and analysis using specified parameters and input files.

## Function Parameters

- **normalize:** Logical. Perform data normalization if TRUE.
- **useBeta:** Logical. If TRUE, expect input data as beta values (scaled between 0 and 1). If FALSE, process raw intensity data.
- **directory:** Directory containing input data files (default: current working directory).
- **arrayType:** Type of DNA methylation array used (options: "27K", "450K", or "EPIC").

## Input Files

The following input files are required:

- **Sample_Sheet.csv:** A CSV file containing phenotypic data for each sample.
- **IDAT Files:** IDAT files containing methylation data for each sample.
- **betaValues.csv:** If IDAT files are not available, processed beta values can be provided.

## Naming Conventions for Data in "Sample_Sheet.csv"

Ensure that the "Sample_Sheet.csv" file follows these naming conventions:

- "Array": Array plate name in format "RXCY" (X: row number, Y: column number).
- "Age": Age values of samples.
- "Sex": Gender information of samples (M or F, 0 or 1).
- "Smoking_Status": Smoking status of samples.
- "Batch": Batch information of samples.
- "Slide": Slide information of samples.
- "Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC": Optional\* cell type proportions if IDAT files are not available.

## Output

The function generates "output.txt" and analysis plots. Be sure to install the necessary R packages. Input CSV files must follow the provided naming conventions.

## Example Usage

```R
main(normalize = TRUE, useBeta = FALSE, directory = "path/to/data", arrayType = "450K")
```
