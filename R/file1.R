globalVariables("myEnv")
myEnv <- new.env(parent = emptyenv())
# Main function starts here
main <- function(directory = getwd(),
                 normalize = TRUE,
                 useBeta = FALSE,
                 arrayType = "450K",
                 useSampleSheet = TRUE) {
    baseDirectory <- getwd()
    setwd(directory)
    startup()
    if (useBeta) {
        message("Reading betaValues.csv...")
        myEnv$bVals <- read.csv("betaValues.csv", row.names = 1)
    } else {
        message("Processing IDAT files...")
        processIDAT(directory, useSampleSheet, arrayType)
    }
    CC <- NULL
    if (!is.numeric(myEnv$rgSet) & arrayType != "27K") {
        CC <- estimateCellCounts(myEnv$rgSet)
    }
    if (useBeta == FALSE) {
        message("Writing extracted beta values...")
        write.csv(myEnv$bVals, file = "extractedBetaValues.csv")
    }
    preparePdataSVs(myEnv$bVals, useSampleSheet, CC, arrayType)
    message("Generating epigenetic age...")
    results <- calculateDNAmAge(myEnv$bVals, myEnv$pdataSVs,
                                normalize)
    clockname <- "DunedinPACE"
    results$DunedinPACE <- calculateDunedinPACE()
    if ("Age" %in% colnames(myEnv$pdataSVs) &&
                            "Sex" %in% colnames(myEnv$pdataSVs)) {
        results <- calculateGrimAge(results)
    } else {
        results$GrimAge <- NA
        results$GrimAgeAccel <- NA
    }
    finalOutput <- ""
    if (useSampleSheet | useBeta == FALSE) {
        finalOutput <- processAllAgeTypes(results)
    }
    exportResults(results, myEnv$bVals, finalOutput)
    setwd(baseDirectory)
}

# Function for loading tools and setting variables
startup <- function() {
    message("Loading dependencies...")
    installDirectory <- paste0(
        path.package("EpigeneticAgePipeline"),
        "/data/"
    )
    assign("bVals", 0, envir = myEnv)
    assign("rgSet", 0, envir = myEnv)
    assign("listofCors", c(), envir = myEnv)
    assign("corsToRemove", c(), envir = myEnv)
    assign("pdataSVs", 0, envir = myEnv)
    assign("exportDf", 0, envir = myEnv)
    assign("outliersCSV", 0, envir = myEnv)
    assign("residualsCSV", 0, envir = myEnv)
    data("PC-clocks", envir = myEnv, package = "EpigeneticAgePipeline")
    data("DunedinPACE", envir = myEnv, package = "EpigeneticAgePipeline")
}

# Function for getting cell counts
estimateCellCounts <- function(rgSet) {
    message("Generating cell counts...")
    FlowSorted.CordBlood.450k::FlowSorted.CordBlood.450k
    CC <- minfi::estimateCellCounts(
        rgSet, compositeCellType = "CordBlood", processMethod = "auto",
        probeSelect = "auto", cellTypes = c("Bcell", "CD4T", "CD8T", "Gran",
                                            "Mono", "nRBC"),
        referencePlatform = c("IlluminaHumanMethylation450k"),
        returnAll = FALSE, meanPlot = FALSE, verbose = TRUE
    )
    return(CC)
}

# Creating pdataSVs dataframe
preparePdataSVs <- function(bVals, useSampleSheet, CC, arrayType) {
    myEnv$pdataSVs <- data.frame(row.names = colnames(bVals))
    if (useSampleSheet) {
        createAnalysisDF(getwd())
        if ("EpiAge" %in% colnames(myEnv$pdataSVs)) {
            myEnv$pdataSVs$EpiAge <- NULL
        }
    }
    if (!is.numeric(myEnv$rgSet) & arrayType != "27K") {
        addCellCountsToPdataSVs(CC)
    }
}

# Adding cell counts to pdataSVs
addCellCountsToPdataSVs <- function(CC) {
    cellTypes <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC")
    for (cellType in cellTypes) {
        myEnv$pdataSVs[[cellType]] <- as.numeric(CC[, cellType])
    }
}

# Calculating epigenetic age
calculateDNAmAge <- function(bVals, pdataSVs, shouldNormalize) {
    results <- NULL
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
    grimDf <- data.frame(Sample = colnames(myEnv$bVals),
        Age = myEnv$pdataSVs$Age,
        Sex = as.character(myEnv$pdataSVs$Sex))
    for (i in seq(from = 1, to = length(grimDf$Sex))) {
        if (grimDf$Sex[i] == "M" | grimDf$Sex[i] == 1) {
            grimDf$Sex[i] <- "Male"
        } else if (grimDf$Sex[i] == "F" | grimDf$Sex[i] == 2) {
            grimDf$Sex[i] <- "Female"
        }
    }
    grimDf$Sex <- as.factor(grimDf$Sex)
    grimage <-
        methyAge(betas = myEnv$bVals, clock = "PCGrimAge",
            age_info = grimDf)

    df$GrimAgeAccel <- grimage$Age_Acceleration
    df$GrimAge <- grimage$mAge
    return(df)
}

# Exporting results
exportResults <- function(results, bVals, finalOutput) {
    if ("Age" %in% colnames(myEnv$pdataSVs) &&
        "Sex" %in% colnames(myEnv$pdataSVs)) {
        plotDf <- preparePlotDf(results, myEnv$pdataSVs)
        plotDf$age <- NULL
        plotDf$grimage <- results$GrimAge
        plotDf$age <- myEnv$pdataSVs$Age
        createGroupedBarChart(plotDf, "sample", "value", "Age_Measure",
            "Sample ID and Type of Age Measure")
        myEnv$exportDf <- results[, c("id", "Horvath", "Hannum",
            "Levine", "GrimAge", "skinHorvath", "DunedinPACE", "GrimAgeAccel", "age")]
    } else if ("Age" %in% colnames(myEnv$pdataSVs)) {
        plotDf <- preparePlotDf(results, myEnv$pdataSVs)
        createGroupedBarChart(plotDf, "sample", "value", "Age_Measure",
                              "Sample ID and Type of Age Measure")
        myEnv$exportDf <- results[, c("id", "Horvath", "Hannum",
                                      "Levine", "skinHorvath", "DunedinPACE", "age")]
    } else {
        myEnv$exportDf <- results[, c("id", "Horvath", "Hannum", "Levine",
                                "skinHorvath", "DunedinPACE")]
    }
    writeResults(finalOutput, myEnv$exportDf, results)
}

# Making data frame for plotting
preparePlotDf <- function(results, pdataSVs) {
    data.frame(
        sample = colnames(myEnv$bVals),
        horvath = results$Horvath,
        skinhorvath = results$skinHorvath,
        hannum = results$Hannum,
        levine = results$Levine,
        age = pdataSVs$Age
    )
}

# Writing out results
writeResults <- function(finalOutput, exportDf, results) {
    formattedResults <- kable(results[,c("id", "Horvath", "skinHorvath",
        "Hannum", "Levine", "GrimAge", "DunedinPACE", "GrimAgeAccel")], format = "markdown")
    exportDf <- as.data.frame(exportDf)
    write.table(exportDf, file = "epigeneticAge.txt")
    write_file(finalOutput, file = "output.txt")
    write(formattedResults, file = "results.md")
}

# Definition for initiation function for process age type
processAllAgeTypes <- function(results) {
    finalOutput <- ""
    ageTypes <- c("Horvath", "skinHorvath", "Hannum", "Levine", "DunedinPACE")
    for (ageType in ageTypes) {
        finalOutput <- processAgeType(results, ageType, finalOutput)
        myEnv$pdataSVs[[ageType]] <- NULL
    }
    if ("Age" %in% colnames(myEnv$pdataSVs) &&
       "Sex" %in% colnames(myEnv$pdataSVs)) {
        finalOutput <- processAgeType(results, "GrimAge", finalOutput)
        myEnv$pdataSVs$GrimAge <- NULL
        finalOutput <- processAgeType(results, "GrimAgeAccel", finalOutput)
        myEnv$pdataSVs$GrimAgeAccel <- NULL
    }
    return(finalOutput)
}

# Definition for function used for generating the grouped bar chart
createGroupedBarChart <- function(data, x, y, fill, title) {
    meltedDf <- melt(data, id.vars = x, variable.name = fill)

    customPalette <- c( "age" = "red", "horvath" = "#66c2a5",
                        "skinhorvath" = "#3288bd", "hannum" = "#5e4fa2", "levine" = "#3288dd"
    )

    if ("Age" %in% colnames(myEnv$pdataSVs) &&
        "Sex" %in% colnames(myEnv$pdataSVs)) {
        customPalette <- c( "age" = "red", "horvath" = "#66c2a5",
                            "skinhorvath" = "#3288bd", "hannum" = "#5e4fa2",
                            "levine" = "#3288dd", "grimage" = "#3798de")
    }

    plot <- ggplot(
        data = meltedDf,
        aes_string(
            x = x,
            y = y,
            fill = fill
        )
    ) +
        geom_bar(stat = "identity", width = 1, position = position_dodge(width = 0.6)) +
        labs(x = x, y = "Age", title = title) +
        scale_fill_manual(values = customPalette) +
        theme_minimal()

    plot <- plot +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 65, hjust = 1) # Adjust x-axis text for better readability
        )

    ggsave("SampleIDandAge.png", plot = plot, width = 35, height = 10, units = "in", dpi = 300)

    for (i in 2:(ncol(data) - 1))
    {
        plot <- ggscatter(data,
            x = "age", y = colnames(data)[i],
            add = "reg.line", #
            add.params = list(color = "blue", fill = "lightgray")
        )
        plot <- plot + stat_cor(method = "pearson")
        ggsave(
            filename = paste0("plot_", colnames(data)[i], ".png"),
            plot = plot,
            width = 1000, height = 1000,
            units = "px"
        )
    }
}

# Definition of processIDAT function starts here
processIDAT <- function(directory, useSampleSheet, arrayType) {
    installDirectory <- paste0(
        path.package("EpigeneticAgePipeline"),
        "/data/"
    )
    dataDirectory <- directory
    myEnv$rgSet <- read.metharray.exp(dataDirectory, force = TRUE)
    #    Calculate    the    detection    p-values
    detP <- detectionP(myEnv$rgSet)
    samples_before <- dim(myEnv$rgSet)[2]
    keep <- colMeans(detP) < 0.05
    myEnv$rgSet <- myEnv$rgSet[, keep]
    samples_removed <- samples_before - dim(detP)[2]
    message(
        "-----    ",
        samples_removed,
        " sample(s) removed due to poor quality"
    )

    mSetSq <- myEnv$rgSet
    mSetSq <- preprocessRaw(mSetSq)
    mSetSq <- mapToGenome(mSetSq)
    mSetSq <- ratioConvert(mSetSq)
    detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
    probes_before <- dim(mSetSq)[1]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    mSetSqFlt <- mSetSq[keep, ]
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message(
        "-----    ",
        probes_removed,
        " probe(s) removed for failing in",
        "one or more samples"
    )
    probes_before <- dim(mSetSqFlt)[1]
    #    Remove    probes    with    SNPs    at    CpG    site
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message(
        "-----    ",
        probes_removed,
        " probe(s)  removed",
        "for having SNPs at CpG site"
    )
    probes_before <- dim(mSetSqFlt)[1]
    #    Exclude    cross    reactive    probes
    if (arrayType == "450K") {
        data("ChenEtAlList", envir = myEnv, package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$ChenEtAlList
    } else if (arrayType == "27K") {
        data("non-specific-probes-Illumina27k", envir = myEnv,
            package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$non_specific_probes_Illumina27k
    } else {
        data("PidsleyCrossReactiveProbesEPIC", envir = myEnv,
            package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$PidsleyCrossReactiveProbesEPIC
    }
    keep <-
        !(featureNames(mSetSqFlt)
        %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep, ]
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message( "-----    ", probes_removed,
            " probe(s) removed for being cross reactive"
    )
    probes_before <- dim(mSetSqFlt)[1]
    #    Remove    Sex    Probes
    if (arrayType == "EPIC") {
        ann <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest
    } else if (arrayType == "450K") {
        ann <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Manifest
    } else {
        ann <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Manifest
    }
    sexProbes <- ann[which(ann$chr %in% c("chrX", "chrY")), ]
    keep <- !(featureNames(mSetSqFlt) %in% sexProbes$Name)
    mSetSqFlt <- mSetSqFlt[keep, ]
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("-----    ",
        probes_removed, " probe(s)  removed", "for  being  on  sex  chromosomes"
    )
    #    Print    out    the    number    of    probes    remaining
    message(
        "-----    ",
        dim(mSetSqFlt)[1],
        " probe(s) remaining for analysis"
    )
    #    Calculate    methylation    beta    values
    myEnv$bVals <- getBeta(mSetSqFlt)
    myEnv$bVals <- myEnv$bVals[,order(colnames(myEnv$bVals))]
}

# Definition for function used in matrix generation
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 2, ...) {
    usr <- par("usr")
    on.exit(par(usr = usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "complete.obs"))
    myEnv$listofCors <- append(myEnv$listofCors, r)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if (missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor)
}

# Definition for function used in matrix generation
mydiag.panel <- function(x, labels, ...) {
    ll <- par("usr")
    rect(ll[1], ll[3], ll[2], ll[4], col = "#CC7178")
}

# Definition for function used in matrix generation
twolines <- function(x, y) {
    points(x, y, pch = 20)
    abline(lm(y ~ x), col = "#CC7178")
}

# Definition for function used for processing epigenetic age measures
processAgeType <- function(data, ageType, output) {
    myEnv$pdataSVs[[ageType]] <- data[[ageType]]
    colnames <- names(myEnv$pdataSVs)
    myEnv$pdataSVs <- myEnv$pdataSVs[c(
        ageType,
        colnames[colnames != ageType]
    )]
    if ("V1" %in% names(myEnv$pdataSVs)) {
        myEnv$pdataSVs$V1 <- NULL
    }
    for (i in colnames(myEnv$pdataSVs)) {
        if (length(unique(myEnv$pdataSVs[[i]])) == 1) {
            myEnv$pdataSVs[[i]] <- NULL
        }
    }
    diag.labels <- colnames(myEnv$pdataSVs)
    pdataColumns <-
        names(myEnv$pdataSVs)[names(myEnv$pdataSVs) != ageType]
    plot.formula <- as.formula(paste(
        ageType, "~",
        paste(pdataColumns,
              collapse = "+"
        )
    ))
    generateMatrixPlot(plot.formula, diag.labels, ageType)
    finalOutput <- paste(output, "\n", ageType, "Covariates\n")
    finalOutput <- corCovariates(finalOutput)
    covariate_data <- myEnv$pdataSVs
    if (ageType != "EpiAge") {
        myEnv$listofCors <- c()
        myEnv$corsToRemove <- c()
    }
    return(finalOutput)
}

# Definition for function used for generating the correlation matrices
generateMatrixPlot <- function(plotFormula, diagLabels, ageType) {
    cairo_pdf(paste("matrixplot", ageType, ".pdf", sep = ""),
            width = 14, height = 14, fallback_resolution = 1000)
    pairs(plotFormula, data = myEnv$pdataSVs, upper.panel = twolines,
            labels = diagLabels, diag.panel = mydiag.panel,
            lower.panel = panel.cor, label.pos = 0.5, main = "")
    dev.off()
}

# Definition for function for creating the df used for analyses
createAnalysisDF <- function(directory) {
    message("Reading Sample_Sheet.csv...")
    setwd(directory)
    sampleData <- read.csv("Sample_Sheet.csv", header = TRUE)
    sampleData <- as.data.frame(sampleData)
    for (i in colnames(sampleData)) {
        newVarName <- gsub("\\...[1-2]$", "", i)
        if (grepl("...1", i)) {
            myEnv$pdataSVs[[newVarName]] <- as.numeric(sampleData[[i]])
        } else if (grepl("...2", i)) {
            if (i == "Array...2") {
                row <- as.factor(gsub(
                    "R(\\d+).*",
                    "\\1",
                    sampleData[[i]]
                ))
                column <- as.factor(gsub(
                    ".*C(\\d+)",
                    "\\1",
                    sampleData[[i]]
                ))
                myEnv$pdataSVs$Row <- row
                myEnv$pdataSVs$Column <- column
            } else {
                myEnv$pdataSVs[[newVarName]] <- as.factor(sampleData[[i]])
            }
        }
    }
    for (i in colnames(myEnv$pdataSVs)) {
        if (length(unique(myEnv$pdataSVs[[i]])) == 1 & i != "V1") {
            message("\nCovariate with only 1 unique level detected,",
                    "consider excluding\n")
        }
    }
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
corCovariates <- function(x) {
    corDf <- myEnv$pdataSVs
    corDf <- corDf[seq.default(from = 1, to = (ncol(myEnv$pdataSVs))), ]
    row.names(corDf) <- colnames(corDf)
    corDf[, ] <- 0
    counter <- 1
    for (i in seq.default(from = 1, to = (ncol(corDf))))
    {
        for (j in seq.default(from = i, to = (ncol(corDf)) - 1))
        {
            corDf[j + 1, i] <- myEnv$listofCors[counter]
            counter <- counter + 1
        }
    }
    corDf <- corDf[-nrow(corDf), ]
    for (i in seq.default(from = 1, to = (ncol(corDf)))) {
        for (j in seq.default(from = i, to = (ncol(corDf)) - 1)) {
            if (!is.na(corDf[j + 1, i])) {
                if (corDf[j + 1, i] > 0.6) {
                    covariate1 <- rownames(corDf)[j + 1]
                    covariate2 <- colnames(corDf)[i]
                    x <- paste(
                        x, "\n", covariate1, " and ", covariate2,
                        " are highly correlated: ", corDf[j + 1, i], "\n"
                    )
                    message("\n")
                    message( covariate1, " and ", covariate2,
                        " are highly correlated: ", corDf[j + 1, i]
                    )
                    message("\n")
                    if (!(covariate1 %in% myEnv$corsToRemove)) {
                        myEnv$corsToRemove <- append(
                            myEnv$corsToRemove,
                            covariate1
                        )
                    }
                    if (!(covariate2 %in% myEnv$corsToRemove)) {
                        myEnv$corsToRemove <- append(
                            myEnv$corsToRemove,
                            covariate2
                        )
                    }
                }
            }
        }
    }
    return(x)
}

# Generating regression formula
formulaGeneration <- function(string, columnsUsed) {
    formula_string <- ""
    if (!("Column" %in% columnsUsed) && "Row" %in% columnsUsed && "Slide"
        %in% columnsUsed && "Batch" %in% columnsUsed) {
        formula_string <- paste(
            "EpiAge",
            " ~ ",
            string,
            " + ",
            "(Row|Slide)",
            " + ",
            "(1|Batch)"
        )
    } else if (!("Slide" %in% columnsUsed) && "Row" %in% columnsUsed && "Column"
               %in% columnsUsed && "Batch" %in% columnsUsed) {
        formula_string <- paste(
            "EpiAge",
            "~",
            string,
            " + ",
            "Row + Column",
            " + ",
            "(1|Batch)"
        )
    } else if ("Row" %in% columnsUsed && "Column" %in% columnsUsed && "Slide"
               %in% columnsUsed && "Batch" %in% columnsUsed) {
        formula_string <- paste(
            "EpiAge",
            "~",
            string,
            " + ",
            "((Row+Column)|Slide)",
            " + ",
            "(1|Batch)"
        )
    } else if (!("Row" %in% columnsUsed) && "Column" %in% columnsUsed && "Slide"
        %in% columnsUsed && "Batch" %in% columnsUsed) {
        formula_string <- paste(
            "EpiAge",
            " ~ ",
            string,
            " + ",
            "(Column|Slide)",
            " + ",
            "(1|Batch)"
        )
    } else {
        formula_string <- paste("EpiAge", "~", string)
    }
    return(formula_string)
}

# Definition for function used for residual generation
residGeneration <- function(pdata) {
    columns <- colnames(pdata)
    columnsUsed <- columns[columns != "EpiAge"]
    string <- paste(columnsUsed, collapse = " + ")
    formula_string <- formulaGeneration(string, columnsUsed)
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
pcaGeneration <- function(PCs) {
    myEnv$bVals <- na.omit(myEnv$bVals)
    bValst <- t(myEnv$bVals)
    bpca <- prcomp(bValst, center = TRUE, scale = FALSE)
    pca_scores <- as.data.frame(bpca$x)
    constant <- 3
    sample_outliers <- c()
    alloutliers <- c()
    if (ncol(pca_scores) < PCs) {
        loopNum <- ncol(pca_scores)
    } else {
        loopNum <- PCs
    }
    for (i in seq.default(from = 1, to = loopNum))
    {
        a <- subset(
            rownames(bpca$x),
            bpca$x[, i] > (mean(bpca$x[, i]) + constant * sd(bpca$x[, i]))
        )
        b <- subset(
            rownames(bpca$x),
            bpca$x[, i] < (mean(bpca$x[, i]) - constant * sd(bpca$x[, i]))
        )
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
generateResiduals <- function(directory = getwd(), useBeta = FALSE,
                            arrayType = "450K", ignoreCor = FALSE, PCs = 5) {
    baseDirectory <- getwd()
    setwd(directory)
    startup()
    if (useBeta == TRUE) {
        message("Reading betaValues.csv...")
        myEnv$bVals <- read.csv("betaValues.csv", row.names = 1)
    } else {
        message("Processing IDAT files...")
        processIDAT(directory, useSampleSheet = TRUE, arrayType)
    }
    if (PCs != 0) {
        pca_scores <- pcaGeneration(PCs = PCs)
    }
    # Processing and Writing Residuals ####
    myEnv$pdataSVs <- as.data.frame(matrix(NA,
        nrow = ncol(myEnv$bVals),
        ncol = 1
    ))
    rownames(myEnv$pdataSVs) <- colnames(myEnv$bVals)
    createAnalysisDF(directory)
    if (PCs != 0) {
        myEnv$pdataSVs <-
            cbind(myEnv$pdataSVs, pca_scores[, seq(from = 1, to = PCs)])
    }
    if (!"EpiAge" %in% colnames(myEnv$pdataSVs)) {
        warning(
            "You did not specify a column called EpiAge@@@1 in your",
            "Sample_Sheet.csv"
        )
        return()
    }
    if (!is.numeric(myEnv$rgSet) & arrayType != "27K") {
        CC <- estimateCellCounts(myEnv$rgSet)
        addCellCountsToPdataSVs(CC)
    }
    processAgeType(myEnv$pdataSVs, "EpiAge", " ")
    x <- corCovariates(" ")
    if (!ignoreCor) {
        removeCovariates()
    }
    myEnv$listofCors <- c()
    myEnv$corsToRemove <- c()
    myEnv$residualsCSV <- residGeneration(myEnv$pdataSVs)
    write.csv(myEnv$outliersCSV, "OutlierSamples.csv")
    write.csv(myEnv$residualsCSV, "Residuals.csv")
    setwd(baseDirectory)
}

#FUNCTION FROM:
#https://github.com/yiluyucheng/dnaMethyAge
methyAge <- function(betas, clock, age_info=NA) {
    if(grepl('^PC[A-Z]', clock)){
        data(list="PC-clocks", envir=myEnv, package = "EpigeneticAgePipeline")
        myEnv$coefs$Coefficient <- myEnv$coefs[[clock]]
    }else{
        data(list="DunedinPACE", envir=myEnv, package = "EpigeneticAgePipeline")
    }
    is_beta <- TRUE
    if(clock == 'DunedinPACE'){
        betas <- preprocessDunedinPACE(betas,
            ref_means=myEnv$gold_standard_means)
    }
    if(is_beta){
        myEnv$coefs <- setNames(myEnv$coefs$Coefficient, myEnv$coefs$Probe)
        betas <- rbind(betas, Intercept=1)
        betas <- betas[rownames(betas) %in% names(myEnv$coefs), ]
        betas[is.na(betas)] <- 0
        m_age <- t(betas) %*% matrix(data=myEnv$coefs[rownames(betas)])
    }
    m_age <- data.frame(Sample=rownames(m_age), mAge=as.numeric(m_age))
    warning_message <- "\n'age_info' should be a dataframe which contains
    sample ID and age information,
    like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge acceleration
    will not be calculated."
    if (is(age_info, "data.frame")) {
        if (all(c('Sample', 'Age') %in% colnames(age_info))){
            m_age <- merge(age_info, m_age, by='Sample', sort=FALSE)
            if (nrow(m_age) < 1){
                stop(message("Colnames of the input beta dataframe do not match
                    any of the values of the 'Sample' column in
                    age_info!"))
            }
            if(clock == 'PCGrimAge'){
                if('Sex' %in% colnames(age_info)){
                    m_age$is_Female <- gsub('^F.*', 1, m_age$Sex)
                    m_age$is_Female <- gsub('^M.*', 0, m_age$is_Female)
                    m_age$is_Female <- as.numeric(m_age$is_Female)
                    m_age$mAge <- m_age$mAge +
                        as.matrix(m_age[, c('is_Female', 'Age')]) %*%
                        myEnv$PCGrimAge_agesex$PCGrimAge
                }else{
                    stop(message("\nTo calculate 'PCGrimage', 'age_info'
                        should include a 'Sex' column that contains
                        binary sex annotation, i.e. either Female or Male."))
                }
            }
            m_age$Age_Acceleration <- getAccel(m_age$Age, m_age$mAge)
        }else{
            warning(message("\nThe colnames of age_info should include both
                'Sample' and 'Age',
                like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge\nAge
                acceleration will not be calculated."))
        }
    } else if (is.na(age_info[1])){
        if(clock == 'PCGrimAge'){
            stop(message("\nTo calculate 'PCGrimage': \n'age_info' should be
                a dataframe which contains sample ID, age,
                sex information,
                like:\nSample\tAge\tSex\nname1\t30\tFemale\nname2\t60\tMale
                \nname3\t40\tFemale\n"))
        }
    } else {
        warning(message(warning_message))
    }
    return(m_age)
}

#FUNCTION FROM:
#https://github.com/yiluyucheng/dnaMethyAge
preprocessDunedinPACE <- function(betas, ref_means, least_proportion=0){
    common_p <- intersect(rownames(betas), names(ref_means))
    back_p <- length(common_p) / length(ref_means)
    if(back_p > least_proportion){
        betas <- betas[common_p, ]
        betas[,] <- normalize.quantiles.use.target(as.matrix(betas),
            target=ref_means[common_p])
        missing_p <- setdiff(names(ref_means), common_p)
        #### replace NA with reference mean
        betas <- data.frame(betas, check.names = FALSE)
        betas[missing_p, ] <- NA
        ref_means <- ref_means[rownames(betas)]
        for(c in 1:ncol(betas)){
            na_col <- which(is.na(betas[, c]))
            betas[na_col, c] <- ref_means[na_col]
        }
        return(betas)
    }else{
        stop(sprintf("Missing too many probes. Only %.2f%s of
                     the required probes have been found!", back_p * 100, "%"))
    }
}

#FUNCTION FROM:
#https://github.com/yiluyucheng/dnaMethyAge
getAccel <- function(c_age, m_age){
    message("Age acceleration is calculated as the residual resulting",
        " from a linear regression model which DNAm age is regressed",
        " on chronological age.") ## copied
    fit_model <- lm(m_age ~ c_age, na.action = na.exclude)
    accel <- residuals(fit_model)
    return(accel)
}

loadTestData <- function() {
    myEnv <- new.env(parent = emptyenv())
    data("CpGNames", package = "EpigeneticAgePipeline", envir = myEnv)
    CpGNames <- myEnv$CpGNames
    df <- data.frame(matrix(runif(length(CpGNames) * 9),
        nrow = length(CpGNames), ncol = 9))
    colnames(df) <- paste0("Sample", 1:9)
    rownames(df) <- CpGNames
    SampleSheet <- data.frame(
        Sex...2 = ifelse(runif(9) < 0.3, "Male", "Female"),
        Age...1 = sample(0:100, 9, replace = TRUE),
        EpiAge...1 = sample(0:100, 9, replace = TRUE)
    )
    write.csv(SampleSheet, paste0(system.file(package =
        "EpigeneticAgePipeline"),"/data/Sample_Sheet.csv"))
    write.csv(df, paste0(system.file(package = "EpigeneticAgePipeline"),
        "/data/betaValues.csv"))
}

removeTestData <- function() {
    filePath <- paste0(system.file(package = "EpigeneticAgePipeline"),
        "/data/betaValues.csv")
    filePath2 <- paste0(system.file(package = "EpigeneticAgePipeline"),
        "/data/Sample_Sheet.csv")
    if (file.exists(filePath)) {
        unlink(filePath)
        unlink(filePath2)
        message("Test data deleted")
    } else {
        message("Test data not yet loaded")
    }
}




