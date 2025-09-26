# Generating regression formula
formulaGeneration <- function(columnsUsed) {
    formula_string <- ""
    if (!("Column" %in% columnsUsed) && "Row" %in% columnsUsed && "Slide"
        %in% columnsUsed && "Batch" %in% columnsUsed) {
        columnsUsed <- setdiff(columnsUsed, c("Row", "Slide", "Batch"))
        string <- paste0(columnsUsed, collapse = " + ")
        formula_string <- paste0(
            "EpiAge",
            " ~ ",
            string,
            " + ",
            "(Row|Slide)",
            " + ",
            "(1|Batch)"
        )
    } else if ("Row" %in% columnsUsed && "Column" %in% columnsUsed && "Slide"
               %in% columnsUsed && "Batch" %in% columnsUsed) {
        columnsUsed <- setdiff(columnsUsed, c("Row", "Column", "Batch", "Slide"))
        string <- paste0(columnsUsed, collapse = " + ")
        formula_string <- paste0(
            "EpiAge",
            " ~ ",
            string,
            " + ",
            "(Row + Column|Slide)",
            " + ",
            "(1|Batch)"
        )
    } else if (!("Row" %in% columnsUsed) && "Column" %in% columnsUsed && "Slide"
               %in% columnsUsed && "Batch" %in% columnsUsed) {
        columnsUsed <- setdiff(columnsUsed, c("Column", "Batch", "Slide"))
        string <- paste0(columnsUsed, collapse = " + ")
        formula_string <- paste0(
            "EpiAge",
            " ~ ",
            string,
            " + ",
            "(Column|Slide)",
            " + ",
            "(1|Batch)"
        )
    } else if ("Batch" %in% columnsUsed) {
        columnsUsed <- setdiff(columnsUsed, c("Batch"))
        string <- paste0(columnsUsed, collapse = " + ")
        formula_string <- paste0(
            "EpiAge",
            " ~ ",
            string,
            " + ",
            "(1|Batch)"
        )
    } else {
        string <- paste0(columnsUsed, collapse = " + ")
        formula_string <- paste("EpiAge", "~", string)
    }
    return(formula_string)
}

# Definition for function used for residual generation
residGeneration <- function(pdata, formula) {
    columns <- colnames(pdata)
    columnsUsed <- columns[columns != "EpiAge"]
    if (is.null(formula)) {
        formula_string <- formulaGeneration(columnsUsed)
    } else {
        formula_string <- formula
    }
    runlme <- function(formula) {
        lme1 <- glmmTMB(formula,
                        data = pdata,
                        family = "gaussian",
                        control = glmmTMBControl(optCtrl = list(
                            iter.max = 10000,
                            eval.max = 10000
                        ))
        )
        smodel <- lme1
        return(smodel)
    }
    lme_formula <- formula_string
    message(lme_formula)
    lme_formula <- as.formula(lme_formula)
    lme_summary <- try(runlme(lme_formula), silent = FALSE)
    resids <- residuals(lme_summary)
    return(resids)
}

# Definition for function specifically used for generating pca's
pcaGeneration <- function(PCs, threshold) {
    myEnv$bVals <- na.omit(myEnv$bVals)
    bValst <- t(myEnv$bVals)
    bpca <- prcomp(bValst, center = TRUE, scale = FALSE)
    pca_scores <- as.data.frame(bpca$x)
    constant <- threshold
    sample_outliers <- c()
    alloutliers <- c()


    if (ncol(pca_scores) < PCs) {
        loopNum <- ncol(pca_scores)
    } else {
        loopNum <- PCs
    }


    for (i in seq.default(from = 1, to = loopNum))
    {
        median <- median(bpca$x[, i])
        MAD <- stats::mad(bpca$x[,i])
        upper <- median + constant * MAD
        lower <- median - constant * MAD

        a <- subset(
            rownames(bpca$x),
            bpca$x[, i] > upper)
        b <- subset(
            rownames(bpca$x),
            bpca$x[, i] < lower)

        out <- c(a, b)
        sample_outliers <- c(sample_outliers, out)
        alloutliers <- c(alloutliers, sample_outliers)
        sample_outliers <- c()
    }

    myEnv$outliersCSV <- unique(alloutliers)
    bpca <- prcomp(bValst, center = TRUE, scale = FALSE)
    pca_scores <- as.data.frame(bpca$x)
    return(pca_scores)
}

# Residual and PCA Generation function
generateResiduals <- function(inputDirectory = getwd(),
                              outputDirectory = inputDirectory,
                              useBeta = FALSE,
                              formula = NULL,
                              arrayType = "450K",
                              sampleSheetFile = "Sample_Sheet.csv",
                              columnTypes = NULL,
                              ignoreCor = TRUE,
                              PCs = 5,
                              threshold = 3,
                              doParallel = TRUE,
                              doCellCounts = TRUE,
                              tissue = "bloodAdult",
                              cellDeconvMethod = "CP",
                              placentaTrimester = "third",
                              detPSampleCutoff = 0.05,
                              detPProbeCutoff = 0.01,
                              minBeads = 3,
                              beadSampleMaxFailFrac = 0.05,
                              beadProbeMaxFailFrac  = 0.05) {

    if ((missing(columnTypes) || is.null(columnTypes) || length(columnTypes) == 0)) {
        stop("columnTypes is required to specify columns and must be a named integer vector: 1=factor, 2=numeric.\n
             eg. c(\"Age\"=2, \"Sex\"=1, \"Batch\"=1, \"BMI\"=2, \"Row\"=1, \"Column\"=1)")
    }
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    startup(outputDirectory)
    setwd(inputDirectory)

    if (useBeta) {
        file <- sort(list.files(pattern = "^betaValues\\.csv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE), decreasing = TRUE)[1]
        if (doParallel) {
            message(paste0("Reading ", file, ", utilizing parallel processing..."))
            myEnv$bVals <- as.data.frame(data.table::fread(file))
            rownames(myEnv$bVals) <- myEnv$bVals$V1
            myEnv$bVals$V1 <- NULL
            myEnv$origSampleNames <- colnames(myEnv$bVals)
        } else {
            message(paste0("Reading ", file, "..."))
            myEnv$bVals <- read.csv(file, row.names = 1)
            myEnv$origSampleNames <- colnames(myEnv$bVals)
        }
    } else {
        message("Processing IDAT files...")
        processIDAT(inputDirectory,
                    arrayType, useSampleSheet = TRUE,
                    sampleSheetFile,
                    detPSampleCutoff,
                    detPProbeCutoff,
                    minBeads,
                    beadSampleMaxFailFrac,
                    beadProbeMaxFailFrac)
    }

    if (PCs != 0) {
        pca_scores <- pcaGeneration(PCs, threshold)
    }

    myEnv$pdataSVs <- data.frame(row.names = myEnv$origSampleNames)
    createAnalysisDF(inputDirectory, sampleSheetFile = sampleSheetFile, columnTypes = columnTypes)

    if (PCs != 0) {
        myEnv$pdataSVs <- cbind(myEnv$pdataSVs, pca_scores[, seq(from = 1, to = PCs)])
    }

    #if (!is.numeric(myEnv$rgSet) && arrayType != "27K" && doCellCounts) {
        CC <- estimateCellCountsFlexible(
            myEnv$rgSet,
            arrayType = arrayType,
            tissue = tissue,
            method = cellDeconvMethod,
            placentaTrimester = placentaTrimester
        )
        assign("CellCountsDf", CC, envir = .GlobalEnv)
        addCellCountsToPdataSvs(CC)
    #}

    if (is.null(formula)) {
        processAgeType(myEnv$pdataSVs, "EpiAge", " ")
    } else {
        processAgeType(myEnv$pdataSVs, gsub(" ", "", sub("~.*", "", formula)), " ")
    }

    validColumns <- c()
    for (i in colnames(myEnv$pdataSVs)) {
        if (length(unique(myEnv$pdataSVs[[i]])) != 1) {
            validColumns <- c(validColumns, i)
        }
    }
    x <- corCovariates(" ", validColumns)
    if (!ignoreCor && is.null(formula)) removeCovariates()
    myEnv$corsToRemove <- c()

    myEnv$residualsCSV <- residGeneration(myEnv$pdataSVs, formula)
    write.csv(myEnv$outliersCSV, paste0(myEnv$baseDirectory, "/OutlierSamples.csv"))
    write.csv(myEnv$residualsCSV, paste0(myEnv$baseDirectory, "/Residuals.csv"))
}

