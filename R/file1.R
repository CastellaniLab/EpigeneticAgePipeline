globalVariables("myEnv")
myEnv <- new.env(parent = emptyenv())
# Main function starts here
main <- function(inputDirectory = getwd(),
                 outputDirectory = inputDirectory,
                 normalize = TRUE,
                 useBeta = FALSE,
                 arrayType = "450K",
                 useSampleSheet = TRUE,
                 sampleSheetFile = "Sample_Sheet.csv",
                 columnTypes = NULL,
                 doParallel = TRUE,
                 writeBeta = TRUE,
                 tissue = "bloodAdult",
                 cellDeconvMethod = "CP",
                 placentaTrimester = "third",
                 useImputation = FALSE,
                 detPSampleCutoff = 0.05,
                 detPProbeCutoff = 0.01,
                 minBeads = 3,
                 beadSampleMaxFailFrac = 0.05,
                 beadProbeMaxFailFrac  = 0.05) {

    if ((missing(columnTypes) || is.null(columnTypes) || length(columnTypes) == 0) && useSampleSheet == TRUE) {
        stop("columnTypes is required if useSampleSheet = TRUE and must be a named integer vector: 1=factor, 2=numeric.\n
             eg. c(\"Age\"=2, \"Sex\"=1, \"Batch\"=1, \"BMI\"=2, \"Row\"=1, \"Column\"=1)")
    }
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    startup(outputDirectory, useImputation)
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
                    arrayType, useSampleSheet,
                    sampleSheetFile,
                    detPSampleCutoff,
                    detPProbeCutoff,
                    minBeads,
                    beadSampleMaxFailFrac,
                    beadProbeMaxFailFrac)
    }
    preparePdataSVs(myEnv$bVals, useSampleSheet, sampleSheetFile, columnTypes)
    #CC <- NULL
    #if (!is.numeric(myEnv$rgSet) & arrayType != "27K") {
        CC <- estimateCellCountsFlexible(
            myEnv$rgSet,
            arrayType = arrayType,
            tissue = tissue,
            method = cellDeconvMethod,
            placentaTrimester = placentaTrimester
        )
        addCellCountsToPdataSvs(CC)
        assign("CellCountsDf", CC, envir = .GlobalEnv)
    #}

    if (!useBeta && writeBeta) {
        message("Writing extracted beta values...")
        write.csv(myEnv$bVals, file = paste0(myEnv$baseDirectory, "/extractedBetaValues.csv"))
    }

    message("Generating epigenetic age...")
    results <- calculateDNAmAge(myEnv$bVals, myEnv$pdataSVs, normalize)
    if ("age" %in% names(results)) {
        results$age <- NULL
    }
    results$DunedinPACE <- calculateDunedinPACE()
    names(results) <- gsub("skinHorvath", "Horvath_SkinBlood", names(results))

    if ("Age" %in% colnames(myEnv$pdataSVs) && "Sex" %in% colnames(myEnv$pdataSVs)) {
        results <- calculateGrimAge(results)
    }
    results <- calculatePC(results)
    results <- recomputeAccels(results)

    finalOutput <- " "
    if (ncol(myEnv$pdataSVs) != 0 && nrow(myEnv$pdataSVs) > 1) {
        finalOutput <- processAllAgeTypes(results)
    } else {
        message("Skipping covariate correlation analysis: only one observation or variable present.")
    }
    exportResults(results, myEnv$bVals, finalOutput)
    clockSummary()
}

# Function for loading tools and setting variables
startup <- function(outputDirectory, useImputation = FALSE) {
    message("Loading dependencies...")

    if ("myEnv" %in% search()) detach("myEnv")
    rm(list = ls(envir = myEnv), envir = myEnv)

    assign("baseDirectory", sub("/$", "", outputDirectory), envir = myEnv)

    assign("bVals", 0, envir = myEnv)
    assign("rgSet", 0, envir = myEnv)
    assign("listofCors", c(), envir = myEnv)
    assign("corsToRemove", c(), envir = myEnv)
    assign("pdataSVs", 0, envir = myEnv)
    assign("exportDf", 0, envir = myEnv)
    assign("outliersCSV", 0, envir = myEnv)
    assign("residualsCSV", 0, envir = myEnv)
    assign("useImputation", useImputation, envir = myEnv)
    assign("PC_Coefs", 0, envir = myEnv)
    assign("DunedinPACE_Coefs", 0, envir = myEnv)
    assign("origSampleNames", 0, envir = myEnv)
    assign("failedSamples", 0, envir = myEnv)

    data("PC-clocks", envir = myEnv, package = "EpigeneticAgePipeline")
    data("DunedinPACE", envir = myEnv, package = "EpigeneticAgePipeline")
    data(list = "GrimAge_V1", envir = myEnv, package = "EpigeneticAgePipeline")
    data(list = "GrimAge_V2", envir = myEnv, package = "EpigeneticAgePipeline")
    attach(myEnv, name = "myEnv")
}

# Creating pdataSVs dataframe
preparePdataSVs <- function(bVals, useSampleSheet, sampleSheetFile, columnTypes) {
    myEnv$pdataSVs <- data.frame(row.names = myEnv$origSampleNames)
    if (isTRUE(useSampleSheet)) {
        createAnalysisDF(getwd(), sampleSheetFile = sampleSheetFile, columnTypes = columnTypes)
    }
}

# Calculating epigenetic age
calculateDNAmAge <- function(bVals, pdataSVs, shouldNormalize) {
    betaValues <- as.matrix(bVals)
    if ("Age" %in% colnames(pdataSVs)) {
        results <- methylclock::DNAmAge(betaValues, normalize = shouldNormalize,
                                        age = pdataSVs$Age)
    } else {
        results <- methylclock::DNAmAge(betaValues, normalize = shouldNormalize)
    }
    return(results)
}

# Calculating dunedinPACE
calculateDunedinPACE <- function() {
    message("Calculating DunedinPACE...")
    clockname <- "DunedinPACE"
    dunedinPACEDf <-
        methyAge(betas = myEnv$bVals, clock = clockname)
    return(as.data.frame(dunedinPACEDf)[, 2])
}

# Calculating GrimAge
calculateGrimAge <- function(df) {
    message("Calculating GrimAge...")
    sampleIds <- colnames(myEnv$bVals)
    grimDf <- data.frame(Sample = colnames(myEnv$bVals),
        Age = myEnv$pdataSVs$Age,
        Sex = as.character(myEnv$pdataSVs$Sex))
    for (i in seq(from = 1, to = length(grimDf$Sex))) {
        if (grimDf$Sex[i] == 1) {
            grimDf$Sex[i] <- "Male"
        } else if (grimDf$Sex[i] == 2) {
            grimDf$Sex[i] <- "Female"
        }
    }

    getAccelVec <- function(age, dnamAge) {
        keep <- !is.na(age) & !is.na(dnamAge)
        out  <- rep(NA_real_, length(dnamAge))
        if (sum(keep) >= 2) out[keep] <- resid(lm(dnamAge[keep] ~ age[keep]))
        out
    }

    grimDf$Sex <- as.factor(grimDf$Sex)
    grimAgePcRes <- methyAge(betas = myEnv$bVals, clock = "PCGrimAge", age_info = grimDf)
    idxPc <- match(sampleIds, grimAgePcRes$Sample)
    df$GrimAge_PC      <- grimAgePcRes$mAge[idxPc]
    #df$GrimAgeAccel_PC <- grimAgePcRes$PCGrimAge_Acceleration[idxPc]

    parsedV1 <- parseGrimAgeTable(myEnv$GrimAge_V1)
    grimAgeV1Res <- computeGrimAgeV1(
        betas   = myEnv$bVals,
        ageInfo = grimDf,
        comps   = parsedV1$components,
        hazard  = parsedV1$hazard,
        cal     = parsedV1$cal
    )
    idxV1 <- match(sampleIds, grimAgeV1Res$Sample)
    grimAgeV1 <- grimAgeV1Res$mAge[idxV1]
    df$GrimAge_V1      <- grimAgeV1
    #df$GrimAgeAccel_V1 <- getAccelVec(grimDf$Age, grimAgeV1)

    parsedV2 <- parseGrimAgeTable(myEnv$GrimAge_V2)

    grimAgeV2Res <- computeGrimAgeV2(
        betas   = myEnv$bVals,
        ageInfo = grimDf,
        comps   = parsedV2$components,
        hazard  = parsedV2$hazard,
        cal     = parsedV2$cal
    )
    idxV2 <- match(sampleIds, grimAgeV2Res$Sample)
    grimAgeV2 <- grimAgeV2Res$mAge[idxV2]
    df$GrimAge_V2      <- grimAgeV2
    #df$GrimAgeAccel_V2 <- getAccelVec(grimDf$Age, grimAgeV2)

    return(df)
}

calculatePC <- function (results) {
    Horvath_PC <- methyAge(betas = myEnv$bVals, clock = "PCHorvathS2013")
    Hannum_PC <- methyAge(betas = myEnv$bVals, clock = "PCHannumG2013")
    Levine_PC <- methyAge(betas = myEnv$bVals, clock = "PCPhenoAge")
    Horvath_SkinBlood_PC <- methyAge(betas = myEnv$bVals, clock = "PCHorvathS2018")
    results$Horvath_PC <- Horvath_PC$mAge
    results$Hannum_PC <- Hannum_PC$mAge
    results$Levine_PC <- Levine_PC$mAge
    results$Horvath_SkinBlood_PC <- Horvath_SkinBlood_PC$mAge
    return(results)
}

# Exporting results
exportResults <- function(results, bVals, finalOutput) {
    # Build exportDf with all numeric columns (clocks + accelerations + age if present)
    numCols <- vapply(results, is.numeric, logical(1))
    exportDf <- results[, numCols, drop = FALSE]
    if ("id" %in% colnames(results)) exportDf <- cbind(id = results$id, exportDf)

    plotCols <- pickPlotClocks(results)

    # Grouped bar chart: Age + clocks if Age exists
    if (length(plotCols) > 0 && "Age" %in% colnames(myEnv$pdataSVs)) {
        plotDf <- data.frame(
            sample = colnames(myEnv$bVals),
            age = myEnv$pdataSVs$Age,
            results[, plotCols, drop = FALSE],
            check.names = FALSE
        )
        createGroupedBarChart(
            data = plotDf,
            x = "sample",
            y = "value",
            fill = "Age_Measure",
            title = "Sample ID and Type of Age Measure"
        )
    } else if (length(plotCols) > 0) {
        plotDf <- data.frame(
            sample = colnames(myEnv$bVals),
            results[, plotCols, drop = FALSE],
            check.names = FALSE
        )
        createGroupedBarChart(
            data = plotDf,
            x = "sample",
            y = "value",
            fill = "Age_Measure",
            title = "Sample ID and Type of Age Measure"
        )
    }

    myEnv$exportDf <- exportDf
    writeResults(finalOutput, myEnv$exportDf, results)
}

# Writing out results
writeResults <- function(finalOutput, exportDf, results) {
    formattedResults <- kable(exportDf, format = "markdown")
    exportDf <- as.data.frame(exportDf)
    assign("EpiAgeResultsDf", exportDf, envir = .GlobalEnv)
    write.table(exportDf, file = paste0(myEnv$baseDirectory,"/epigeneticAge.txt"))
    write_file(finalOutput, file = paste0(myEnv$baseDirectory,"/output.txt"))
    write(formattedResults, file = paste0(myEnv$baseDirectory,"/results.md"))
}

# Definition for initiation function for process age type
processAllAgeTypes <- function(results) {
    finalOutput <- ""
    ageTypes <- getClockColumns(results)
    ageTypes <- c(ageTypes, grep("accel|ageAcc", colnames(results), value = TRUE, ignore.case = TRUE))

    ageTypes <- setdiff(ageTypes, "age")  # age is a covariate, not a clock

    for (ageType in ageTypes) {
        finalOutput <- processAgeType(results, ageType, finalOutput)
        if (ageType %in% names(myEnv$pdataSVs)) myEnv$pdataSVs[[ageType]] <- NULL
    }
    return(finalOutput)
}

# Definition for function used for processing epigenetic age measures
processAgeType <- function(data, ageType, output) {
    myEnv$pdataSVs[[ageType]] <- data[[ageType]]
    colnames <- names(myEnv$pdataSVs)
    myEnv$pdataSVs <- myEnv$pdataSVs[c(
        ageType,
        colnames[colnames != ageType]
    )]

    write.csv(myEnv$pdataSVs, file = paste0(myEnv$baseDirectory, "/", ageType,"SampleData.csv"))
    validColumns <- c()
    for (i in colnames(myEnv$pdataSVs)) {
        if (length(unique(myEnv$pdataSVs[[i]])) != 1) {
            validColumns <- c(validColumns, i)
        }
    }
    diag.labels <- validColumns
    pdataColumns <- setdiff(validColumns, ageType)
    plot.formula <- as.formula(paste(
        ageType, "~",
        paste(pdataColumns,
              collapse = "+"
        )
    ))
    generateMatrixPlot(plot.formula, diag.labels, ageType)
    finalOutput <- paste(output, "\n", ageType, "Covariates\n")
    finalOutput <- corCovariates(finalOutput, validColumns)
    if (ageType != "EpiAge") {
        myEnv$corsToRemove <- c()
    }
    return(finalOutput)
}

# Definition for function for creating the df used for analyses
createAnalysisDF <- function(directory,
                             sampleSheetFile = "Sample_Sheet.csv",
                             columnTypes) {
    if (is.null(names(columnTypes)) || any(names(columnTypes) == "")) {
        stop("columnTypes must be a *named* integer vector. Names are column names in the sample sheet.")
    }
    if (!all(columnTypes %in% c(1L, 2L))) {
        stop("columnTypes values must be 1 (factor) or 2 (numeric).")
    }

    setwd(directory)
    if (!file.exists(sampleSheetFile)) {
        stop("Sample sheet file not found: ", file.path(directory, sampleSheetFile))
    }
    message("Reading ", sampleSheetFile)
    sampleData <- read.csv(sampleSheetFile, header = TRUE, check.names = FALSE)
    sampleData <- as.data.frame(sampleData)

    # Validate requested columns (allow Row/Column derived from 'Array')
    reqCols <- setdiff(names(columnTypes), c("Row","Column"))
    missingCols <- setdiff(reqCols, colnames(sampleData))
    if (length(missingCols) > 0) {
        stop("Missing required column(s) in sample sheet: ", paste(missingCols, collapse = ", "))
    }

    # Build pdataSVs strictly from requested columns
    for (nm in names(columnTypes)) {
        flag <- as.integer(columnTypes[[nm]])
        if (nm == "Array") {
            arr <- as.character(sampleData[["Array"]])
            rowVals <- as.factor(gsub("R(\\d+).*", "\\1", arr))
            colVals <- as.factor(gsub(".*C(\\d+)", "\\1", arr))
            myEnv$pdataSVs$Row <- rowVals
            myEnv$pdataSVs$Column <- colVals
            next
        }
        vec <- sampleData[[nm]]
        if (flag == 2L) {
            message("Reading variable as numeric: ", nm)
            myEnv$pdataSVs[[nm]] <- suppressWarnings(as.numeric(vec))
        } else { # 1L
            message("Reading variable as factor: ", nm)
            myEnv$pdataSVs[[nm]] <- as.factor(vec)
        }
    }

    # Warn on one-level covariates (still allowed)
    for (nm in names(myEnv$pdataSVs)) {
        if (length(unique(myEnv$pdataSVs[[nm]])) == 1) {
            message("\nCovariate with only 1 unique level detected (", nm, "); consider excluding.\n")
        }
    }

    if (!is.null(myEnv$failedSamples) && length(myEnv$failedSamples) > 0) {
        keep_rows <- !(rownames(myEnv$pdataSVs) %in% myEnv$failedSamples)
        removed_n <- sum(!keep_rows)
        myEnv$pdataSVs <- myEnv$pdataSVs[keep_rows, , drop = FALSE]
        message("Removed ", removed_n, " failed sample(s) from pdataSVs.")
    }

}

getClockColumns <- function(results) {
    # keep numeric columns that have at least one non-NA
    isNum <- vapply(results, is.numeric, logical(1))
    keep <- names(results)[isNum]
    # Drop obvious non-clock columns
    drop <- c("id")  # 'age' stays for plotting alongside clocks
    setdiff(keep, drop)
}

pickPlotClocks <- function(results) {
    # 1) numeric columns only
    is_num <- vapply(results, is.numeric, logical(1))
    cols <- names(results)[is_num]

    # 2) drop acceleration columns of any type
    accel_pat <- "(^|\\.)ageAcc(2|3)?\\.|accel"
    cols <- cols[!grepl(accel_pat, cols, ignore.case = TRUE)]

    # 3) drop obvious non-clocks
    cols <- setdiff(cols, c("Age", "age", "id", "DunedinPACE"))

    # keep original order
    cols
}

# Definition for function used for user specified covariate removal
removeCovariates <- function() {
    if (length(myEnv$corsToRemove) == 0) {
        message("No significant correlations")
        return()
    }
    message(
        "The following covariates ",
        "were found to be highly correlated: \n"
    )
    for (i in myEnv$corsToRemove)
    {
        message(i)
    }
    message(
        "\nTo remove ",
        "one of the covariates or several, ",
        "enter 1 ",
        "to  remove, ",
        "0   to  keep"
    )
    userInput <- scan(file = "", n = length(myEnv$corsToRemove))
    userInput <- as.numeric(userInput)
    for (i in seq.default(from = 1, to = length(userInput))) {
        if (userInput[i] == 1) {
            column <- myEnv$corsToRemove[i]
            message(column)
            myEnv$pdataSVs <- myEnv$pdataSVs[
                ,
                !(names(myEnv$pdataSVs) %in% column)
            ]
        }
        message(i)
    }
}

# Definition for function used to find highly correlated covariates
corCovariates <- function(reportText, validColumns) {
    # Subset pdataSVs to validColumns for correlation
    corDf <- myEnv$pdataSVs[, validColumns, drop = FALSE]

    # Compute correlation matrix directly, safely handling factors and NAs
    corDfNumeric <- data.frame(
        lapply(corDf, function(col) {
            if (is.factor(col) || is.character(col)) {
                as.numeric(as.factor(col))
            } else {
                as.numeric(col)
            }
        }),
        check.names = FALSE,
        row.names = rownames(corDf)
    )

    corMatrix <- cor(corDfNumeric, use = "pairwise.complete.obs")

    # Find and report high correlations
    for (i in seq_len(ncol(corMatrix) - 1)) {
        for (j in seq(i + 1, ncol(corMatrix))) {
            corrValue <- corMatrix[i, j]
            if (!is.na(corrValue) && abs(corrValue) > 0.6) {
                covariate1 <- colnames(corMatrix)[i]
                covariate2 <- colnames(corMatrix)[j]
                reportText <- paste0(
                    reportText, "\n", covariate1, " and ", covariate2,
                    " are highly correlated: ", round(corrValue, 3), "\n"
                )
                myEnv$corsToRemove <- union(myEnv$corsToRemove, c(covariate1, covariate2))
            }
        }
    }
    return(reportText)
}

recomputeAccels <- function(results) {
    # need Age to compute any acceleration
    if (!"Age" %in% names(myEnv$pdataSVs)) return(results)

    # find cell counts (prefer myEnv, fall back to .GlobalEnv)
    CC <- if (exists("CellCountsDf", envir = myEnv, inherits = FALSE)) {
        get("CellCountsDf", envir = myEnv)
    } else if (exists("CellCountsDf", envir = .GlobalEnv, inherits = FALSE)) {
        get("CellCountsDf", envir = .GlobalEnv)
    } else {
        warning("CellCountsDf not found; leaving accelerations unchanged.")
        return(results)
    }

    # meffil-style tiny-variance filter
    ok <- which(apply(CC, 2, stats::IQR, na.rm = TRUE) > 1e-6)
    if (length(ok)) CC <- CC[, ok, drop = FALSE]

    sids <- as.character(results$id)

    # Build the design data.frame expected by methylclock:::ageAcc2
    df <- data.frame(
        age = myEnv$pdataSVs$Age,
        CC,
        check.names = FALSE
    )

    is_num <- vapply(results, is.numeric, logical(1))
    candidates <- names(results)[is_num]
    candidates <- setdiff(candidates, c("Age", "DunedinPACE"))
    candidates <- candidates[!grepl("accel|ageAcc", candidates, ignore.case = TRUE)]
    labs <- candidates

    # Always skip DunedinPACE even if passed in explicitly
    labs <- setdiff(labs, "DunedinPACE")

    for (lab in labs) {
        if (!lab %in% names(results)) next
        x <- data.frame(Sample = sids, mAge = results[[lab]])
        a <- methylclock:::ageAcc2(x, df = df, lab = lab)

        for (nm in paste0(c("ageAcc", "ageAcc2", "ageAcc3"), ".", lab)) {
            results[[nm]] <- a[[nm]]
        }
    }

    results
}


clockSummary <- function() {
    data("PC-clocks", envir = myEnv, package = "EpigeneticAgePipeline")

    cpg.names <- rownames(myEnv$bVals)

    coefHannum_PC <- data.frame(CpGmarker = myEnv$coefs$Probe,
                                PCHannum2013 = myEnv$coefs$PCHannumG2013)
    coefHorvath2013_PC <- data.frame(CpGmarker = myEnv$coefs$Probe,
                                     PCHorvaths2013 = myEnv$coefs$PCHorvathS2013)
    coefHorvath2018_PC <- data.frame(CpGmarker = myEnv$coefs$Probe,
                                     PCHorvaths2018 = myEnv$coefs$PCHorvathS2018)
    coefPhenoAge_PC <- data.frame(CpGmarker = myEnv$coefs$Probe,
                                  PCHPhenoAge = myEnv$coefs$PCPhenoAge)
    coefGrimAge_PC <- data.frame(CpGmarker = myEnv$coefs$Probe,
                                 PCGrimAge = myEnv$coefs$PCGrimAge)

    data("DunedinPACE", envir = myEnv, package = "EpigeneticAgePipeline")

    coefDunedinPACE <- data.frame(CpGmarker = myEnv$coefs$Probe,
                                 PCGrimAge = myEnv$coefs$Coefficient)


    checkHorvath <- coefHorvath$CpGmarker[-1][!coefHorvath$CpGmarker[-1] %in% cpg.names]
    checkHannum <- coefHannum$CpGmarker[!coefHannum$CpGmarker %in% cpg.names]
    checkLevine <- coefLevine$CpGmarker[-1][!coefLevine$CpGmarker[-1] %in% cpg.names]
    checkSkin <- coefSkin$CpGmarker[-1][!coefSkin$CpGmarker[-1] %in% cpg.names]
    checkPedBE <- coefPedBE$CpGmarker[-1][!coefPedBE$CpGmarker[-1] %in% cpg.names]
    checkWu <- coefWu$CpGmarker[-1][!coefWu$CpGmarker[-1] %in% cpg.names]
    checkTL <- coefTL$CpGmarker[-1][!coefTL$CpGmarker[-1] %in% cpg.names]
    checkBLUP <- coefBLUP$CpGmarker[-1][!coefBLUP$CpGmarker[-1] %in% cpg.names]
    checkEN <- coefEN$CpGmarker[-1][!coefEN$CpGmarker[-1] %in% cpg.names]
    checkDunedin <- coefDunedinPACE$CpGmarker[-1][!coefDunedinPACE$CpGmarker[-1] %in% cpg.names]
    checkGrimAge_PC <-coefGrimAge_PC$CpGmarker[-1][!coefGrimAge_PC$CpGmarker[-1] %in% cpg.names]
    checkHorvath_PC <- coefHorvath2013_PC$CpGmarker[-1][!coefHorvath2013_PC$CpGmarker[-1] %in% cpg.names]
    checkHannum_PC <- coefHannum_PC$CpGmarker[-1][!coefHannum_PC$CpGmarker[-1] %in% cpg.names]
    checkLevine_PC <- coefPhenoAge_PC$CpGmarker[-1][!coefPhenoAge_PC$CpGmarker[-1] %in% cpg.names]
    checkSkin_PC <- coefHorvath2018_PC$CpGmarker[-1][!coefHorvath2018_PC$CpGmarker[-1] %in% cpg.names]

    parsedV1 <- parseGrimAgeTable(myEnv$GrimAge_V1)
    parsedV2 <- parseGrimAgeTable(myEnv$GrimAge_V2)
    covListV1 <- lapply(parsedV1$components, coverageReport, betas = myEnv$bVals)
    covDfV1 <- data.frame(
        Surrogate = names(covListV1),
        present   = sapply(covListV1, `[[`, "present"),
        total     = sapply(covListV1, `[[`, "total"),
        pct       = sapply(covListV1, `[[`, "pct"),
        row.names = NULL
    )
    covDfV1$Surrogate <- paste0(covDfV1$Surrogate, "_GrimAgeV1")
    covDfV1$pct <- 1 - covDfV1$pct
    covListV2 <- lapply(parsedV2$components, coverageReport, betas = myEnv$bVals)
    covDfV2 <- data.frame(
        Surrogate = names(covListV2),
        present   = sapply(covListV2, `[[`, "present"),
        total     = sapply(covListV2, `[[`, "total"),
        pct       = sapply(covListV2, `[[`, "pct"),
        row.names = NULL
    )
    covDfV2$Surrogate <- paste0(covDfV2$Surrogate, "_GrimAgeV2")
    covDfV2$pct <- 1 - covDfV2$pct


    sizes <- lengths(list(
        checkHorvath,           checkHorvath_PC,
        checkHannum,            checkHannum_PC,
        checkLevine,            checkLevine_PC,
        checkSkin,              checkSkin_PC,
        checkPedBE,
        checkWu,
        checkTL,
        checkBLUP,
        checkEN,
        checkDunedin,
        checkGrimAge_PC
    ))

    n <- c(
        nrow(coefHorvath[-1, ]),        nrow(coefHorvath2013_PC[-1, ]),
        nrow(coefHannum),               nrow(coefHannum_PC[-1, ]),
        nrow(coefLevine[-1, ]),         nrow(coefPhenoAge_PC[-1, ]),
        nrow(coefSkin[-1, ]),           nrow(coefHorvath2018_PC[-1, ]),
        nrow(coefPedBE[-1, ]),
        nrow(coefWu[-1, ]),
        nrow(coefTL[-1, ]),
        nrow(coefBLUP[-1, ]),
        nrow(coefEN[-1, ]),
        nrow(coefDunedinPACE[-1, ]),
        nrow(coefGrimAge_PC[-1, ])
    )

    df <- data.frame(
        clock = c(
            "Horvath", "Horvath_PC",
            "Hannum",  "Hannum_PC",
            "Levine",  "Levine_PC",
            "Horvath_SkinBlood", "Horvath_SkinBlood_PC",
            "PedBE",
            "Wu",
            "TL",
            "BLUP",
            "EN",
            "DunedinPACE",
            "GrimAge_PC"
        ),
        Cpgs_in_clock = n,
        missing_CpGs = sizes,
        percentage = round((sizes / n) * 100, 1)
    )

    # Reformat covDfV1 to match df
    covDfV1_renamed <- data.frame(
        clock = covDfV1$Surrogate,
        Cpgs_in_clock = covDfV1$total,
        missing_CpGs  = covDfV1$total - covDfV1$present,
        percentage    = round((1 - (covDfV1$present/covDfV1$total)) * 100, 1),
        stringsAsFactors = FALSE
    )

    # Do the same for V2
    covDfV2_renamed <- data.frame(
        clock = covDfV2$Surrogate,
        Cpgs_in_clock = covDfV2$total,
        missing_CpGs  = covDfV2$total - covDfV2$present,
        percentage    = round((1 - (covDfV2$present/covDfV2$total)) * 100, 1),
        stringsAsFactors = FALSE
    )

    # Bind everything together
    df_all <- rbind(df, covDfV1_renamed, covDfV2_renamed)

    assign("ClockStats", df_all, envir = myEnv)

    print(df_all)
}






