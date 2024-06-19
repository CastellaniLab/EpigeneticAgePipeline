library("testthat")
test_that("Testing betaValues file processing", {
    library(EpigeneticAgePipeline)
    directory <- paste0(path.package("EpigeneticAgePipeline"),"/data/")
    main(
        directory = directory,
        normalize = TRUE,
        useBeta = TRUE,
        arrayType = "450K",
        useSampleSheet = FALSE
    )
    testthat::expect_true(nrow(exportDf) > 0,
        "Number of rows should be greater than 0")
})
