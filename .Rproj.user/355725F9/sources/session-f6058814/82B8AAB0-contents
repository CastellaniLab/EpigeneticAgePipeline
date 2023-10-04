\name{main}
\alias{main}

\title{main}

\description{
  The `DNAm-Age Pipeline` is a comprehensive function designed for processing and analyzing DNA methylation data. It performs data normalization and analysis using specified parameters and input files.
}

\usage{
  main(directory = getwd(), normalize = TRUE, useBeta = FALSE, arrayType = "450k")
}

\arguments{
  \item{directory}{
    Directory containing input data files (default: current working directory).
  }

  \item{normalize}{
    Logical. Perform data normalization if `TRUE`.
  }

  \item{useBeta}{
    Logical. If `TRUE`, expect input data as beta values (scaled between 0 and 1). If `FALSE`, process raw intensity data.
  }

  \item{arrayType}{
    Type of DNA methylation array used (options: "27K", "450K", or "EPIC").
  }
}

\details{
  Input Files:

  `Sample_Sheet.csv`:
  .csv file containing phenotypic data for each sample

  `IDAT Files`:
  IDAT files containing methylation data for each sample

  `betaValues.csv`:
  If IDAT files are not available, processed beta values can be provided

  Naming Conventions for Data in "Sample_Sheet.csv":

  `Array`:
  Array plate name in format "RXCY" (X: row number, Y: column number).

  "Array" column should be in "RXCY" format for row and column data extraction.

  `Age`:
  Age values of samples.

  `Sex`:
  Gender information of samples (M or F, 0 or 1)

  `Smoking_Status`:
  Smoking status of samples.

  `Batch`:
  Batch information of samples.

  `Slide`:
  Slide information of samples.

  `Bcell`, `CD4T`, `CD8T`, `Gran`, `Mono`, `nRBC`:
  Optional* cell type proportions if IDAT files are not available
}

\value{
  The function returns results and analysis plots.
}

\references{
  % Put references to the literature/web site here
}

\author{
  % Author information
}

\note{
  Example Usage:

  main(directory = "path/to/data", normalize = TRUE, useBeta = FALSE, arrayType = "450K")
}

\seealso{
  % Objects to See Also as \code{\link{help}}
}

\examples{

}
