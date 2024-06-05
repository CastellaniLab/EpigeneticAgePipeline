utils::globalVariables(c(
    "bVals", "rgSet", "listofCors", "corsToRemove",
    "methyAge", "pdataSVs", "exportDf", "outliersCSV",
    "residualsCSV", "accelResidualsCSV"
))

# Main function starts here
main <- function(directory = getwd(),
                 normalize = TRUE,
                 useBeta = FALSE,
                 arrayType = "450K",
                 useSampleSheet = TRUE) {
    setwd(directory)
    startup()

    if (useBeta) {
        bVals <- getBetaValues()
    } else {
        processIDAT(directory, useSampleSheet, arrayType)
    }
    CC <- NULL
    if (!is.numeric(.GlobalEnv$rgSet) & arrayType != "27K") {
        CC <- estimateCellCounts(.GlobalEnv$rgSet)
    }

    write.csv(.GlobalEnv$bVals, file = "extractedBetaValues.csv")

    preparePdataSVs(.GlobalEnv$bVals, useSampleSheet, CC, arrayType)

    results <- calculateDNAmAge(.GlobalEnv$bVals, .GlobalEnv$pdataSVs,
                                                                    normalize)
    clockname <- "DunedinPACE"
    results$DunedinPACE <- calculateDunedinPACE(.GlobalEnv$bVals)

    results$GrimAge <- if ("Age" %in% colnames(.GlobalEnv$pdataSVs) &&
                            "Sex" %in% colnames(.GlobalEnv$pdataSVs)) {
        calculateGrimAge(.GlobalEnv$bVals, .GlobalEnv$pdataSVs)
    } else {
        NA
    }
    finalOutput <- ""
    if (useSampleSheet | useBeta == FALSE) {
        finalOutput <- processAllAgeTypes(results)
    }

    exportResults(results, .GlobalEnv$bVals, finalOutput)
}

#Function for loading tools and setting variables
startup <- function() {
    installDirectory <- paste0(
        path.package("EpigeneticAgePipeline"),
        "/extdata/"
    )
    devtools::load_all(paste0(installDirectory, "dnaMethyAge-main"))
    if (!exists("bVals")) {
        base::assign("bVals", 0, envir = .GlobalEnv)
    }
    base::assign("rgSet", 0, envir = .GlobalEnv)
    base::assign("listofCors", c(), envir = .GlobalEnv)
    base::assign("corsToRemove", c(), envir = .GlobalEnv)
    base::assign("pdataSVs", 0, envir = .GlobalEnv)
    base::assign("exportDf", 0, envir = .GlobalEnv)
    base::assign("outliersCSV", 0, envir = .GlobalEnv)
    base::assign("residualsCSV", 0, envir = .GlobalEnv)
    load(paste0(installDirectory, "PC-clocks.rda"), envir = .GlobalEnv)
    load(paste0(installDirectory, "golden_ref.rda"), envir = .GlobalEnv)
}

#Getting beta values from .csv file
getBetaValues <- function() {
    if (typeof(.GlobalEnv$bVals) != "list") {
        .GlobalEnv$bVals <- read.csv("betaValues.csv", row.names = 1)
    }
    return(.GlobalEnv$bVals)
}

#Function for getting cell counts
estimateCellCounts <- function(rgSet) {
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

#Creating pdataSVs dataframe
preparePdataSVs <- function(bVals, useSampleSheet, CC, arrayType) {
    .GlobalEnv$pdataSVs <- data.frame(row.names = colnames(bVals))
    if (useSampleSheet) {
        createAnalysisDF(getwd())
        if ("EpiAge" %in% colnames(.GlobalEnv$pdataSVs)) {
            .GlobalEnv$pdataSVs$EpiAge <- NULL
        }
    }
    if (!is.numeric(.GlobalEnv$rgSet) & arrayType != "27K") {
        addCellCountsToPdataSVs(CC)
    }
}

#Adding cell counts to pdataSVs
addCellCountsToPdataSVs <- function(CC) {
    cellTypes <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC")
    for (cellType in cellTypes) {
        .GlobalEnv$pdataSVs[[cellType]] <- as.numeric(CC[, cellType])
    }
}

#Calculating epigenetic age
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

calculateDunedinPACE <- function(bVals) {
    clockname <- "DunedinPACE"
    dunedinPACEDf <- methyAge(betas = bVals, clock = clockname, do_plot = FALSE)
    return(as.data.frame(dunedinPACEDf)[, 2])
}

calculateGrimAge <- function(bVals, pdataSVs) {
    grimDf <- data.frame(Sample = colnames(bVals), Age = pdataSVs$Age,
        Sex = as.character(pdataSVs$Sex))
    for (i in seq(from = 1, to = length(grimDf$Sex))) {
        if (grimDf$Sex[i] == "M" | grimDf$Sex[i] == 1) {
            grimDf$Sex[i] <- "Male"
        } else if (grimDf$Sex[i] == "F" | grimDf$Sex[i] == 2) {
            grimDf$Sex[i] <- "Female"
        }
    }
    grimDf$Sex <- as.factor(grimDf$Sex)
    clockname <- "PCGrimAge"
    grimage <- methyAge(betas = bVals, clock = clockname, age_info = grimDf,
                do_plot = FALSE)
    return(grimage$Age_Acceleration)
}

exportResults <- function(results, bVals, finalOutput) {
    if ("Age" %in% colnames(.GlobalEnv$pdataSVs)) {
        plotDf <- preparePlotDf(results, .GlobalEnv$pdataSVs)
        createGroupedBarChart(plotDf, "sample", "value", "Age_Measure",
            "Sample ID and Type of Age Measure")
        .GlobalEnv$exportDf <- results[, c("id", "Horvath", "Hannum",
            "Levine", "skinHorvath", "DunedinPACE", "GrimAge", "age")]
    } else {
        .GlobalEnv$exportDf <- results[, c("id", "Horvath", "Hannum", "Levine",
                                "skinHorvath", "DunedinPACE", "GrimAge")]
    }
    writeResults(finalOutput, .GlobalEnv$exportDf, results)
}

preparePlotDf <- function(results, pdataSVs) {
    data.frame(
        sample = colnames(.GlobalEnv$bVals),
        horvath = results$Horvath,
        skinhorvath = results$skinHorvath,
        hannum = results$Hannum,
        levine = results$Levine,
        age = pdataSVs$Age
    )
}

writeResults <- function(finalOutput, exportDf, results) {
    formattedResults <- knitr::kable(results[,c("id", "Horvath", "skinHorvath",
        "Hannum", "Levine", "DunedinPACE", "GrimAge")], format = "markdown")
    exportDf <- as.data.frame(exportDf)
    write.table(exportDf, file = "epigeneticAge.txt")
    readr::write_file(finalOutput, file = "output.txt")
    write(formattedResults, file = "results.md")
}

processAllAgeTypes <- function(results) {
    finalOutput <- ""
    ageTypes <- c("Horvath", "skinHorvath", "Hannum", "Levine", "DunedinPACE")
    for (ageType in ageTypes) {
        finalOutput <- processAgeType(results, ageType, finalOutput)
        print(finalOutput)
        .GlobalEnv$pdataSVs[[ageType]] <- NULL
    }
    if ("Age" %in% colnames(.GlobalEnv$pdataSVs) &&
       "Sex" %in% colnames(.GlobalEnv$pdataSVs)) {
        finalOutput <- processAgeType(results, "GrimAge", finalOutput)
        print(finalOutput)
        .GlobalEnv$pdataSVs$GrimAge <- NULL
    }
    return(finalOutput)
}

# Definition for function used for generating the grouped bar chart
createGroupedBarChart <- function(data, x, y, fill, title) {
    melted_df <- reshape2::melt(data, id.vars = x, variable.name = fill)
    custom_palette <- c( "age" = "red", "horvath" = "#66c2a5",
        "skinhorvath" = "#3288bd", "hannum" = "#5e4fa2", "levine" = "#3288dd"
    )
    plot <- ggplot2::ggplot(
        data = melted_df,
        ggplot2::aes_string(
            x = x,
            y = y,
            fill = fill
        )
    ) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::labs(x = x, y = "Age", title = title) +
        ggplot2::scale_fill_manual(values = custom_palette) +
        ggplot2::theme_minimal()

    plot <- plot +
        ggplot2::theme(
            plot.background =
                ggplot2::element_rect(fill = "white")
        )

    ggplot2::ggsave("SampleIDandAge.png", plot = plot, width = 10000,
                    height = 1000, units = "px", dpi = 300
    )

    for (i in 2:(ncol(data) - 1))
    {
        plot <- ggpubr::ggscatter(data,
            x = "age", y = colnames(data)[i],
            add = "reg.line", #
            add.params = list(color = "blue", fill = "lightgray")
        )
        plot <- plot + ggpubr::stat_cor(method = "pearson")
        ggplot2::ggsave(
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
        "/extdata/"
    )
    dataDirectory <- directory

    rgSet <- minfi::read.metharray.exp(dataDirectory, force = TRUE)

    #    Calculate    the    detection    p-values
    detP <- minfi::detectionP(rgSet)

    samples_before <- dim(rgSet)[2]

    keep <- colMeans(detP) < 0.05
    rgSet <- rgSet[, keep]

    samples_removed <- samples_before - dim(detP)[2]
    message(
        "-----    ",
        samples_removed,
        " sample(s) removed due to poor quality"
    )

    mSetSq <- rgSet

    mSetSq <- minfi::preprocessRaw(mSetSq)

    mSetSq <- minfi::mapToGenome(mSetSq)

    mSetSq <- minfi::ratioConvert(mSetSq)

    detP <- detP[match(minfi::featureNames(mSetSq), rownames(detP)), ]

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
    mSetSqFlt <- minfi::dropLociWithSnps(mSetSqFlt)

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
        xReactiveProbes <- read.csv(file = paste(installDirectory,
            "ChenEtAlList.csv",
            sep = ""
        ))
    } else if (arrayType == "27K") {
        xReactiveProbes <- read.csv(file = paste(installDirectory,
            "non-specific-probes-Illumina27k.csv",
            sep = ""
        ))
    } else {
        xReactiveProbes <- read.csv(file = paste(installDirectory,
            "PidsleyCrossReactiveProbesEPIC.csv",
            sep = ""
        ))
    }

    keep <-
        !(minfi::featureNames(mSetSqFlt)
        %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep, ]

    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message(
        "-----    ",
        probes_removed,
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
    keep <- !(minfi::featureNames(mSetSqFlt) %in% sexProbes$Name)
    mSetSqFlt <- mSetSqFlt[keep, ]

    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message(
        "-----    ",
        probes_removed,
        " probe(s)  removed",
        "for  being  on  sex  chromosomes"
    )

    #    Print    out    the    number    of    probes    remaining
    message(
        "-----    ",
        dim(mSetSqFlt)[1],
        " probe(s) remaining for analysis"
    )

    #    Calculate    methylation    beta    values
    .GlobalEnv$rgSet <- rgSet
    .GlobalEnv$bVals <- minfi::getBeta(mSetSqFlt)
    .GlobalEnv$bVals <- bVals[,order(colnames(bVals))]
}

# Definition for function used in matrix generation
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 2, ...) {
    usr <- par("usr")
    on.exit(par(usr = usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "complete.obs"))
    .GlobalEnv$listofCors <- append(listofCors, r)
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
    .GlobalEnv$pdataSVs[[ageType]] <- data[[ageType]]
    colnames <- names(.GlobalEnv$pdataSVs)
    .GlobalEnv$pdataSVs <- .GlobalEnv$pdataSVs[c(
        ageType,
        colnames[colnames != ageType]
    )]
    if ("V1" %in% names(pdataSVs)) {
        .GlobalEnv$pdataSVs$V1 <- NULL
    }
    for (i in colnames(.GlobalEnv$pdataSVs)) {
        if (length(unique(.GlobalEnv$pdataSVs[[i]])) == 1) {
            .GlobalEnv$pdataSVs[[i]] <- NULL
        }
    }
    diag.labels <- colnames(pdataSVs)
    pdataColumns <- names(pdataSVs)[names(pdataSVs) != ageType]
    plot.formula <- as.formula(paste(
        ageType, "~",
        paste(pdataColumns,
              collapse = "+"
        )
    ))
    generateMatrixPlot(plot.formula, diag.labels, ageType)
    finalOutput <- paste(output, "\n", ageType, "Covariates\n")
    finalOutput <- corCovariates(finalOutput)
    covariate_data <- .GlobalEnv$pdataSVs
    if (ageType != "EpiAge") {
        .GlobalEnv$listofCors <- c()
        .GlobalEnv$corsToRemove <- c()
    }
    return(finalOutput)
}

# Definition for function used for generating the correlation matrices
generateMatrixPlot <- function(plotFormula, diagLabels, ageType) {
    grDevices::cairo_pdf(paste("matrixplot", ageType, ".pdf", sep = ""),
            width = 14, height = 14, fallback_resolution = 1000)
    pairs(plotFormula, data = .GlobalEnv$pdataSVs, upper.panel = twolines,
            labels = diagLabels, diag.panel = mydiag.panel,
            lower.panel = panel.cor, label.pos = 0.5, main = "")
    grDevices::dev.off()
}

# Definition for function for creating the df used for analyses
createAnalysisDF <- function(directory) {
    setwd(directory)
    sampleData <- read.csv("Sample_Sheet.csv", header = TRUE)
    sampleData <- as.data.frame(sampleData)
    for (i in colnames(sampleData)) {
        newVarName <- gsub("\\...[1-2]$", "", i)
        if (grepl("...1", i)) {
            .GlobalEnv$pdataSVs[[newVarName]] <- as.numeric(sampleData[[i]])
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
                .GlobalEnv$pdataSVs$Row <- row
                .GlobalEnv$pdataSVs$Column <- column
            } else {
                .GlobalEnv$pdataSVs[[newVarName]] <- as.factor(sampleData[[i]])
            }
        }
    }
    for (i in colnames(.GlobalEnv$pdataSVs)) {
        if (length(unique(.GlobalEnv$pdataSVs[[i]])) == 1) {
            message("Covariate with only 1 unique level detected,
                    consider excluding")
        }
    }
}

# Definition for function used for user specified covariate removal
removeCovariates <- function() {
    if (length(corsToRemove) == 0) {
        message("No significant correlations")
        return()
    }
    message(
        "The following covariates",
        "were found to be highly correlated: \n"
    )
    for (i in corsToRemove)
    {
        message(i)
    }
    message(
        "\nTo remove",
        "one of the covariates or several,",
        "enter 1 for",
        "each index  want to  remove,",
        "0   to  keep"
    )
    userInput <- scan(file = "", n = length(corsToRemove))
    userInput <- as.numeric(userInput)
    for (i in seq.default(from = 1, to = length(userInput))) {
        if (userInput[i] == 1) {
            column <- corsToRemove[i]
            message(column)
            .GlobalEnv$pdataSVs <- .GlobalEnv$pdataSVs[
                ,
                !(names(.GlobalEnv$pdataSVs) %in% column)
            ]
        }
        message(i)
    }
}

# Definition for function used to find highly correlated covariates
corCovariates <- function(x) {
    corDf <- pdataSVs
    corDf <- corDf[seq.default(from = 1, to = (ncol(pdataSVs))), ]
    row.names(corDf) <- colnames(corDf)
    corDf[, ] <- 0
    counter <- 1
    for (i in seq.default(from = 1, to = (ncol(corDf))))
    {
        for (j in seq.default(from = i, to = (ncol(corDf)) - 1))
        {
            corDf[j + 1, i] <- listofCors[counter]
            counter <- counter + 1
        }
    }
    corDf <- corDf[-nrow(corDf), ]
    print(.GlobalEnv$listofCors)
    print(corDf)
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
                    if (!(covariate1 %in% corsToRemove)) {
                        .GlobalEnv$corsToRemove <- append(
                            corsToRemove,
                            covariate1
                        )
                    }
                    if (!(covariate2 %in% corsToRemove)) {
                        .GlobalEnv$corsToRemove <- append(
                            corsToRemove,
                            covariate2
                        )
                    }
                }
            }
        }
    }
    return(x)
}

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
            "(Row&Column)",
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
            "(1|Slide)",
            " + ",
            "(Row+Column|Slide)",
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
        lme1 <- glmmTMB::glmmTMB(formula,
                            data = pdata,
                            family = "gaussian",
                            control = glmmTMB::glmmTMBControl(optCtrl = list(
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
pcaGeneration <- function() {
    bVals <- na.omit(bVals)
    bValst <- t(bVals)
    bpca <- prcomp(bValst, center = TRUE, scale = FALSE)
    pca_scores <- as.data.frame(bpca$x)
    constant <- 3
    sample_outliers <- c()
    alloutliers <- c()
    if (ncol(pca_scores) < 5) {
        loopNum <- ncol(pca_scores)
    } else {
        loopNum <- 5
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
    .GlobalEnv$outliersCSV <- unique(alloutliers)
    bValst <- t(bVals)
    bpca <- prcomp(bValst, center = TRUE, scale = FALSE)
    pca_scores <- as.data.frame(bpca$x)
    return(pca_scores)
}

# Residual and PCA Generation function
generateResiduals <- function(directory = getwd(), useBeta = FALSE,
                            arrayType = "450K") {
    startup()
    if (useBeta == TRUE) {
        if (typeof(bVals) != "list") {
            bVals <- read.csv("betaValues.csv", row.names = 1)
        }
    } else {
        processIDAT(directory, useSampleSheet = TRUE, arrayType)
    }
    pca_scores <- pcaGeneration()
    # Processing and Writing Residuals ####
    .GlobalEnv$pdataSVs <- as.data.frame(matrix(NA,
        nrow = ncol(bVals),
        ncol = 1
    ))
    rownames(.GlobalEnv$pdataSVs) <- colnames(bVals)
    createAnalysisDF(directory)
    .GlobalEnv$pdataSVs <- cbind(pdataSVs, pca_scores[, seq(from = 1, to = 5)])
    if (!"EpiAge" %in% colnames(.GlobalEnv$pdataSVs)) {
        warning(
            "You did not specify a column called EpiAge@@@1 in your",
            "Sample_Sheet.csv"
        )
        return()
    }
    processAgeType(.GlobalEnv$pdataSVs, "EpiAge", " ")
    message(corCovariates(" "))
    removeCovariates()
    .GlobalEnv$listofCors <- c()
    .GlobalEnv$corsToRemove <- c()
    .GlobalEnv$residualsCSV <- residGeneration(.GlobalEnv$pdataSVs)
    .GlobalEnv$pdataSVs$EpiAge <- .GlobalEnv$residualsCSV
    .GlobalEnv$accelResidualsCSV <- residGeneration(.GlobalEnv$pdataSVs)
    write.csv(.GlobalEnv$outliersCSV, "OutlierSamples.csv")
    write.csv(.GlobalEnv$residualsCSV, "Residuals.csv")
    write.csv(accelResidualsCSV, "ResidualsAcceleration.csv")
}
