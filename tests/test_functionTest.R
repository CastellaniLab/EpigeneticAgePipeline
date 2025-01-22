library("testthat")
test_that("Testing betaValues file processing", {
    library(EpigeneticAgePipeline)
    EpigeneticAgePipeline:::loadTestData()
    directory <- paste0(path.package("EpigeneticAgePipeline"),"/data/")
    setwd(directory)
    main(
        directory = directory,
        normalize = TRUE,
        useBeta = TRUE,
        arrayType = "450K",
        useSampleSheet = FALSE
    )
    testthat::expect_true(file.exists(paste0(directory, "epigeneticAge.txt")),
        "epigeneticAge.txt should exist")
    EpigeneticAgePipeline:::removeTestData()
})
