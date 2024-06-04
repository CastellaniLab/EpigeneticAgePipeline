library("testthat")

test_that("Testing IDAT file processing", {
    library(EpigeneticAgePipeline)
    directory <- paste0(path.package("EpigeneticAgePipeline"),"/extdata/")
    if (file.exists(paste0(directory, "betaValues.csv")))
    {
        main(
            directory = directory,
            normalize = TRUE,
            useBeta = FALSE,
            arrayType = "450K",
            useSampleSheet = FALSE
        )

        expect_true(nrow(exportDf) > 0,
            "Number of rows should be greater than 0")
    }
})

test_that("Testing betaValues file processing", {
    library(EpigeneticAgePipeline)
    directory <- paste0(path.package("EpigeneticAgePipeline"),"/extdata/")
    if (file.exists(paste0(directory, "betaValues.csv")))
    {
        main(
            directory = directory,
            normalize = TRUE,
            useBeta = TRUE,
            arrayType = "450K",
            useSampleSheet = FALSE
        )

        expect_true(nrow(exportDf) > 0,
                    "Number of rows should be greater than 0")
    }
})
