globalVariables("myEnv")
myEnv <- new.env(parent = emptyenv())
# Main function starts here
main <- function(inputDirectory = getwd(),
                 outputDirectory = inputDirectory,
                 normalize = TRUE,
                 useBeta = FALSE,
                 arrayType = "450K",
                 useSampleSheet = TRUE,
                 doParallel = TRUE,
                 writeBeta = TRUE,
                 useAdult = FALSE,
                 useImputation = FALSE) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    startup(outputDirectory, useImputation)
    setwd(inputDirectory)
    if (useBeta) {
        if (doParallel) {
            file <- sort(list.files(pattern = "^betaValues\\.csv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE), decreasing = TRUE)[1]
            message(paste0("Reading ", file, ", utilizing parallel processing..."))
            myEnv$bVals <- as.data.frame(data.table::fread(file))
            rownames(myEnv$bVals) <- myEnv$bVals$V1
            myEnv$bVals$V1 <- NULL
        } else {
            file <- sort(list.files(pattern = "^betaValues\\.csv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE), decreasing = TRUE)[1]
            message(paste0("Reading ", file, "..."))
            myEnv$bVals <- read.csv(file, row.names = 1)
        }
    } else {
        message("Processing IDAT files...")
        processIDAT(inputDirectory, arrayType, useSampleSheet)
    }
    CC <- NULL
    if (!is.numeric(myEnv$rgSet) & arrayType != "27K") {
        CC <- estimateCellCounts(myEnv$rgSet, arrayType, useAdult)
        assign("CellCountsDf", CC, envir = .GlobalEnv)
    }
    if (useBeta == FALSE) {
        if (writeBeta ) {
            message("Writing extracted beta values...")
            write.csv(myEnv$bVals, file = paste0(myEnv$baseDirectory, "/extractedBetaValues.csv"))
        }
    }
    preparePdataSVs(myEnv$bVals, useSampleSheet, CC, arrayType, useAdult)
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
    finalOutput <- " "
    if (ncol(myEnv$pdataSVs) != 0 && nrow(myEnv$pdataSVs) > 1) {
        finalOutput <- processAllAgeTypes(results)
    } else {
        message("Skipping covariate correlation analysis: only one observation or variable present.")
    }
    exportResults(results, myEnv$bVals, finalOutput)
}

# Function for loading tools and setting variables
startup <- function(outputDirectory, useImputation = FALSE) {
    message("Loading dependencies...")

    if ("myEnv" %in% search()) {
        detach("myEnv")
    }
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
    data("PC-clocks", envir = myEnv, package = "EpigeneticAgePipeline")
    data("DunedinPACE", envir = myEnv, package = "EpigeneticAgePipeline")
    attach(myEnv, name = "myEnv")
}

# Function for getting cell counts
estimateCellCounts <- function(rgSet, arrayType, useAdult) {
    message("Generating cell counts...")
    if (useAdult) {
        if (arrayType == "MSA") {
            if (arrayType == "MSA" &&
                (minfi::annotation(rgSet)["array"] != "IlluminaHumanMethylationMSA" ||
                 minfi::annotation(rgSet)["annotation"] != "ilm10a1.hg38")) {
                minfi::annotation(rgSet) <- c(array = "IlluminaHumanMethylationMSA",
                                              annotation = "ilm10a1.hg38")
            }
            Betas <- getBeta(minfi::preprocessNoob(rgSet))
            IDOLOptimizedCpGsBlood <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs[
                which(FlowSorted.Blood.EPIC::IDOLOptimizedCpGs %in% rownames(Betas))
            ]
            CC <- FlowSorted.Blood.EPIC::projectCellType_CP(
                Betas[IDOLOptimizedCpGsBlood, ],
                FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBlood, ],
                contrastWBC = NULL,
                nonnegative = TRUE,
                lessThanOne = FALSE
            )
        } else if (arrayType == "EPICv2") {
            Betas <- getBeta(minfi::preprocessNoob(rgSet))
            Betas <- sesame::betasCollapseToPfx(Betas)
            IDOLOptimizedCpGsBlood <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs[
                which(FlowSorted.Blood.EPIC::IDOLOptimizedCpGs %in% rownames(Betas))
            ]
            CC <- FlowSorted.Blood.EPIC::projectCellType_CP(
                Betas[IDOLOptimizedCpGsBlood, ],
                FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBlood, ],
                contrastWBC = NULL,
                nonnegative = TRUE,
                lessThanOne = FALSE
            )
        } else {
            FlowSorted.Blood.450k::FlowSorted.Blood.450k
            CC <- minfi::estimateCellCounts(rgSet)
        }
        return(CC)
    } else {
        if (arrayType == "MSA") {
            if (arrayType == "MSA" &&
                (minfi::annotation(rgSet)["array"] != "IlluminaHumanMethylationMSA" ||
                 minfi::annotation(rgSet)["annotation"] != "ilm10a1.hg38")) {
                minfi::annotation(rgSet) <- c(array = "IlluminaHumanMethylationMSA",
                                              annotation = "ilm10a1.hg38")
            }
            Betas <- getBeta(minfi::preprocessNoob(rgSet))
            IDOLOptimizedCpGsBlood <- FlowSorted.CordBloodCombined.450k::IDOLOptimizedCpGsCordBlood[
                which(FlowSorted.CordBloodCombined.450k::IDOLOptimizedCpGsCordBlood %in% rownames(Betas))
            ]
            CC <- FlowSorted.Blood.EPIC::projectCellType_CP(
                Betas[IDOLOptimizedCpGsBlood, ],
                FlowSorted.CordBloodCombined.450k::FlowSorted.CordBloodCombined.450k.compTable[IDOLOptimizedCpGsBlood, ],
                contrastWBC = NULL,
                nonnegative = TRUE,
                lessThanOne = FALSE
            )
        } else if (arrayType == "EPICv2") {
            Betas <- getBeta(minfi::preprocessNoob(rgSet))
            Betas <- sesame::betasCollapseToPfx(Betas)
            IDOLOptimizedCpGsBlood <- FlowSorted.CordBloodCombined.450k::IDOLOptimizedCpGsCordBlood[
                which(FlowSorted.CordBloodCombined.450k::IDOLOptimizedCpGsCordBlood %in% rownames(Betas))
            ]
            CC <- FlowSorted.Blood.EPIC::projectCellType_CP(
                Betas[IDOLOptimizedCpGsBlood, ],
                FlowSorted.CordBloodCombined.450k::FlowSorted.CordBloodCombined.450k.compTable[IDOLOptimizedCpGsBlood, ],
                contrastWBC = NULL,
                nonnegative = TRUE,
                lessThanOne = FALSE
            )
        } else {
            FlowSorted.CordBlood.450k::FlowSorted.CordBlood.450k
            CC <- minfi::estimateCellCounts(rgSet, compositeCellType = "CordBlood", cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC"))
        }
        return(CC)
    }
}

# Creating pdataSVs dataframe
preparePdataSVs <- function(bVals, useSampleSheet, CC, arrayType, useAdult) {
    myEnv$pdataSVs <- data.frame(row.names = colnames(bVals))
    if (useSampleSheet) {
        createAnalysisDF(getwd())
        #if ("EpiAge" %in% colnames(myEnv$pdataSVs)) {
            #myEnv$pdataSVs$EpiAge <- NULL
        #}
    }
    if (!is.numeric(myEnv$rgSet) & arrayType != "27K") {
        addCellCountsToPdataSVs(CC, arrayType, useAdult)
    }
}

# Adding cell counts to pdataSVs
addCellCountsToPdataSVs <- function(CC, arrayType, useAdult) {
    if (useAdult && (arrayType == "EPICv2" || arrayType == "MSA")) {
        cellTypes <- c("Bcell", "CD4T", "CD8T", "Mono", "Neu", "NK")
    } else if (useAdult) {
        cellTypes <- c("CD8T","CD4T", "NK","Bcell","Mono","Gran")
    } else if (useAdult == FALSE && (arrayType == "EPICv2" || arrayType == "MSA")) {
        cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
    } else {
        cellTypes <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC")
    }
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
        if (grimDf$Sex[i] == 1) {
            grimDf$Sex[i] <- "Male"
        } else if (grimDf$Sex[i] == 2) {
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
        myEnv$exportDf <- results[, c("id", "Horvath", "ageAcc.Horvath",
                                      "ageAcc2.Horvath", "ageAcc3.Horvath",
                                      "Hannum", "ageAcc.Levine", "ageAcc2.Levine",
                                      "ageAcc3.Levine","Levine", "GrimAge", "skinHorvath",
                                      "DunedinPACE", "GrimAgeAccel", "age")]
    } else if ("Age" %in% colnames(myEnv$pdataSVs)) {
        plotDf <- preparePlotDf(results, myEnv$pdataSVs)
        createGroupedBarChart(plotDf, "sample", "value", "Age_Measure",
                              "Sample ID and Type of Age Measure")
        myEnv$exportDf <- results[, c("id", "Horvath", "ageAcc.Horvath",
                                      "ageAcc2.Horvath", "ageAcc3.Horvath",
                                      "Hannum", "ageAcc.Levine", "ageAcc2.Levine",
                                      "ageAcc3.Levine","Levine", "skinHorvath",
                                      "DunedinPACE", "age")]
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
    assign("EpiAgeResultsDf", exportDf, envir = .GlobalEnv)
    write.table(exportDf, file = paste0(myEnv$baseDirectory,"/epigeneticAge.txt"))
    write_file(finalOutput, file = paste0(myEnv$baseDirectory,"/output.txt"))
    write(formattedResults, file = paste0(myEnv$baseDirectory,"/results.md"))
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
        geom_bar(stat = "identity", width = 1, position = ggplot2::position_dodge(width = 0.6)) +
        labs(x = x, y = "Age", title = title) +
        scale_fill_manual(values = customPalette) +
        theme_minimal()

    plot <- plot +
        ggplot2::theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = ggplot2::element_text(angle = 65, hjust = 1) # Adjust x-axis text for better readability
        )

    ggsave(paste0(myEnv$baseDirectory,"/SampleIDandAge.png"), plot = plot, width = 35, height = 10, units = "in", dpi = 300)

    for (i in 2:(ncol(data) - 1))
    {
        plot <- ggscatter(data,
            x = "age", y = colnames(data)[i],
            add = "reg.line", #
            add.params = list(color = "blue", fill = "lightgray")
        )
        plot <- plot + stat_cor(method = "pearson")
        ggsave(
            filename = paste0(myEnv$baseDirectory, "/plot_", colnames(data)[i], ".png"),
            plot = plot,
            width = 1000, height = 1000,
            units = "px"
        )
    }
}

# Definition of processIDAT function starts here
processIDAT <- function(directory, arrayType, useSampleSheet) {
    dataDirectory <- directory
    myEnv$rgSet <- read.metharray.exp(dataDirectory, force = TRUE)
    if (useSampleSheet) {
        message("Checking dimensionality")
        sampleData <- read.csv("Sample_Sheet.csv", header = TRUE)
        message("Number of samples in Sample_Sheet:", nrow(sampleData))
        message("Number of samples from IDAT Files:", length(list.files(path = dataDirectory, pattern = "\\.idat$")) / 2)
        if ((length(list.files(path = dataDirectory, pattern = "\\.idat$")) / 2) != nrow(sampleData)) {
            message("The number of samples in Sample_Sheet is not equal to the number of samples")
            stop()
        }
    }

    if (arrayType == "MSA" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylationMSA" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilm10a1.hg38")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylationMSA",
                                      annotation = "ilm10a1.hg38")
    }

    if (arrayType == "EPICv2" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylationEPICv2" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "20a1.hg38")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylationEPICv2",
                                            annotation = "20a1.hg38")
    }

    if (arrayType == "EPIC" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylationEPIC" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilm10b4.hg19")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylationEPIC",
                                            annotation = "ilm10b4.hg19")
    }

    if (arrayType == "450K" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylation450k" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilmn12.hg19")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylation450k",
                                            annotation = "ilmn12.hg19")
    }

    if (arrayType == "27K" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylation27k" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilmn12.hg19")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylation450k",
                                            annotation = "ilmn12.hg19")
    }

    #    Calculate    the    detection    p-values
    detP <- detectionP(myEnv$rgSet)
    detP <- detP[complete.cases(detP),]
    samples_before <- dim(myEnv$rgSet)[2]
    keep <- colMeans(detP) < 0.05
    myEnv$rgSet <- myEnv$rgSet[, keep]
    samples_removed <- samples_before - dim(detP)[2]
    if (samples_removed != 0) {
        message("The following samples were found to be too poor quality, please remove from analysis:")
        for (name in names(keep[keep == FALSE])) {
            message(name)
        }
        stop()
    }
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
        " probe(s) removed for failing in ",
        "one or more samples"
    )
    probes_before <- dim(mSetSqFlt)[1]
    #    Remove    probes    with    SNPs    at    CpG    site
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message(
        "-----    ",
        probes_removed,
        " probe(s)  removed ",
        "for having SNPs at CpG site"
    )
    probes_before <- dim(mSetSqFlt)[1]
    #    Exclude    cross    reactive    probes
    if (arrayType == "450K") {
        data("ChenEtAlList", envir = myEnv, package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$ChenEtAlList
    } else if (arrayType == "27K") {
        data("non_specific_probes_Illumina27k", envir = myEnv,
            package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$non_specific_probes_Illumina27k
    } else if (arrayType == "EPIC") {
        data("PidsleyCrossReactiveProbesEPIC", envir = myEnv,
            package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$PidsleyCrossReactiveProbesEPIC
    } else if (arrayType == "EPICv2" || arrayType == "MSA") {
        data("epicV2CR", envir = myEnv,
             package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$epicV2CR
    }

    keep <-
        !(sub("_.*", "", featureNames(mSetSqFlt))
        %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep, ]
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message( "-----    ", probes_removed,
            " probe(s) removed for being cross reactive"
    )
    probes_before <- dim(mSetSqFlt)[1]
    #    Remove    Sex    Probes
    if (arrayType == "EPIC") {
        ann <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
    } else if (arrayType == "EPICv2") {
        ann <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations
    } else if (arrayType == "450K") {
        ann <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
    } else if (arrayType == "MSA") {
        ann <- IlluminaHumanMethylationMSAanno.ilm10a1.hg38::Locations
    } else {
        ann <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Locations
    }
    sexProbes <- ann[which(ann$chr %in% c("chrX", "chrY")), ]
    keep <- !(featureNames(mSetSqFlt) %in% rownames(sexProbes))
    mSetSqFlt <- mSetSqFlt[keep, ]
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("-----    ",
        probes_removed, " probe(s)  removed ", "for  being  on  sex  chromosomes"
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
    if (arrayType == "EPICv2") {
        suffix <- rownames(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other)
        wSuffixEPIC <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other$EPICv1_Loci
        wSuffix450 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other$Methyl450_Loci
        wSuffix27 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other$Methyl27_Loci
        matchedIndices <- match(rownames(myEnv$bVals), suffix)
        rownames(myEnv$bVals) <- wSuffixEPIC[matchedIndices]
        emptyIndices <- which(rownames(myEnv$bVals) == "")
        rownames(myEnv$bVals)[emptyIndices] <- wSuffix450[matchedIndices][emptyIndices]
        emptyIndices <- which(rownames(myEnv$bVals) == "")
        rownames(myEnv$bVals)[emptyIndices] <- wSuffix27[matchedIndices][emptyIndices]
        myEnv$bVals <- myEnv$bVals[rownames(myEnv$bVals) != "", ]
    }
}

# Definition for function used in matrix generation
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 2, ...) {
    usr <- par("usr")
    on.exit(par(usr = usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "complete.obs"))
    #myEnv$listofCors <- append(myEnv$listofCors, r)
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
    #if ("V1" %in% names(myEnv$pdataSVs)) {
        #myEnv$pdataSVs$V1 <- NULL
    #}
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
        #myEnv$listofCors <- c()
        myEnv$corsToRemove <- c()
    }
    return(finalOutput)
}

# Definition for function used for generating the correlation matrices
generateMatrixPlot <- function(plotFormula, diagLabels, ageType) {
    cairo_pdf(paste0(myEnv$baseDirectory, "/matrixplot", ageType, ".pdf"),
            width = 14, height = 14, fallback_resolution = 1000)
    pairs(plotFormula, data = myEnv$pdataSVs, upper.panel = twolines,
            labels = diagLabels, diag.panel = mydiag.panel,
            lower.panel = panel.cor, label.pos = 0.5, main = "")
    dev.off()
}

# Definition for function for creating the df used for analyses
createAnalysisDF <- function(directory) {
    message("Reading Sample_Sheet.csv (assuming data is sorted A-Z on sample id's)...")
    setwd(directory)
    sampleData <- read.csv("Sample_Sheet.csv", header = TRUE)
    sampleData <- as.data.frame(sampleData)
    for (i in colnames(sampleData)) {
        newVarName <- gsub("\\...[1-2]$", "", i)
        if (grepl("\\.\\.\\.1", i)) {
            message(paste0("Reading variable as numeric, ", i))
            myEnv$pdataSVs[[newVarName]] <- as.numeric(sampleData[[i]])
        } else if (grepl("\\.\\.\\.2", i)) {
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
                message(paste0("Reading Row and Column data from Array"))
                myEnv$pdataSVs$Row <- row
                myEnv$pdataSVs$Column <- column
            } else {
                message(paste0("Reading variable as factor, ", i))
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
                message("\n", covariate1, " and ", covariate2,
                        " are highly correlated: ", round(corrValue, 3), "\n")
                myEnv$corsToRemove <- union(myEnv$corsToRemove, c(covariate1, covariate2))
            }
        }
    }
    return(reportText)
}


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
                            outputDirectory = inputDirectory, useBeta = FALSE, formula = NULL,
                            arrayType = "450K", ignoreCor = TRUE, PCs = 5, threshold = 3,
                            doParallel = TRUE, doCellCounts = TRUE, useAdult = FALSE) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    startup(outputDirectory)
    setwd(inputDirectory)
    if (useBeta == TRUE) {
        if (doParallel) {
            file <- sort(list.files(pattern = "^betaValues\\.csv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE), decreasing = TRUE)[1]
            message(paste0("Reading ", file, ", utilizing parallel processing..."))
            myEnv$bVals <- as.data.frame(data.table::fread(file))
            rownames(myEnv$bVals) <- myEnv$bVals$V1
            myEnv$bVals$V1 <- NULL
        } else {
            file <- sort(list.files(pattern = "^betaValues\\.csv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE), decreasing = TRUE)[1]
            message(paste0("Reading ", file, "..."))
            myEnv$bVals <- read.csv(file, row.names = 1)
        }
    } else {
        message("Processing IDAT files...")
        processIDAT(inputDirectory, arrayType, useSampleSheet = TRUE)
    }
    if (PCs != 0) {
        pca_scores <- pcaGeneration(PCs, threshold)
    }
    # Processing and Writing Residuals ####
    myEnv$pdataSVs <- data.frame(row.names = colnames(myEnv$bVals))
    createAnalysisDF(inputDirectory)
    if (PCs != 0) {
        myEnv$pdataSVs <-
            cbind(myEnv$pdataSVs, pca_scores[, seq(from = 1, to = PCs)])
    }
    if (!"EpiAge" %in% colnames(myEnv$pdataSVs) & !is.null(formula)) {
        warning(
            "You did not specify a column called EpiAge@@@1 in your",
            "Sample_Sheet.csv", ", rename or use custom formula"
        )
        return()
    }
    if (!is.numeric(myEnv$rgSet) && arrayType != "27K" && doCellCounts) {
        CC <- estimateCellCounts(myEnv$rgSet, arrayType, useAdult)
        assign("CellCountsDf", CC, envir = .GlobalEnv)
        addCellCountsToPdataSVs(CC, arrayType, useAdult)
    }
    if (is.null(formula)) {
        processAgeType(myEnv$pdataSVs, "EpiAge", " ")
    } else {
        processAgeType(myEnv$pdataSVs, gsub(" ", "", sub("~.*", "", formula)), " ")
    }
    x <- corCovariates(" ")
    if (!ignoreCor && is.null(formula)) {
        removeCovariates()
    }
    #myEnv$listofCors <- c()
    myEnv$corsToRemove <- c()
    myEnv$residualsCSV <- residGeneration(myEnv$pdataSVs, formula)
    write.csv(myEnv$outliersCSV, paste0(myEnv$baseDirectory, "/OutlierSamples.csv"))
    write.csv(myEnv$residualsCSV, paste0(myEnv$baseDirectory, "/Residuals.csv"))
}

#FUNCTION FROM:
#https://github.com/yiluyucheng/dnaMethyAge
# All borrowed functions were carefully curated and modified to suit specific project requirements.
# Certain parts of the original implementations were excluded to focus on essential functionality.
# Changes include trimming unused features and restructuring for improved integration.
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

        if (myEnv$useImputation) {
            message("Performing imputation")
            data(list='golden_ref', envir=myEnv, package = "EpigeneticAgePipeline")
            ref_mean <- setNames(myEnv$golden_ref$Mean, rownames(myEnv$golden_ref))
            ref_mean <- ref_mean[names(ref_mean) %in% names(myEnv$coefs)]
            betas <- meanImputation(mt=betas, ref=ref_mean, only_ref_rows=FALSE)
        } else {
            betas[is.na(betas)] <- 0
        }

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

#FUNCTION FROM:
#https://github.com/yiluyucheng/dnaMethyAge
meanImputation <- function(mt, ref, cut_off=0.9, only_ref_rows=TRUE){
    mt <- as.matrix(mt)
    if(only_ref_rows){
        mt <- mt[rownames(mt) %in% names(ref), ]
    }
    t_mt <- rowSums(mt)
    row_mean <- rowMeans(mt[is.na(t_mt),], na.rm=TRUE)
    na_row <- intersect(rownames(mt)[is.na(t_mt)], names(ref))
    for (p in na_row){
        n_miss <- sum(is.na(mt[p,]))
        if(n_miss / nrow(mt) > cut_off) {
            #mt[p,][is.na(mt[p,])] <- ref[p]
            mt[p,] <- ref[p]
        }else{
            mt[p,][is.na(mt[p,])] <- row_mean[p]
        }

    }
    ## replace missing rows with reference values
    miss_row <- setdiff(names(ref), rownames(mt))
    mt <- rbind(mt, matrix(rep(ref[miss_row], ncol(mt)), ncol=ncol(mt),
                           dimnames = list(miss_row, colnames(mt))))

    return(mt)
}

#loadTestData <- function() {
    #myEnv <- new.env(parent = emptyenv())
    #data("CpGNames", package = "EpigeneticAgePipeline", envir = myEnv)
    #CpGNames <- myEnv$CpGNames
    #df <- data.frame(matrix(runif(length(CpGNames) * 9),
        #nrow = length(CpGNames), ncol = 9))
    #colnames(df) <- paste0("Sample", 1:9)
    #rownames(df) <- CpGNames
    #SampleSheet <- data.frame(
        #Sex...2 = ifelse(runif(9) < 0.3, "Male", "Female"),
        #Age...1 = sample(0:100, 9, replace = TRUE),
        #EpiAge...1 = sample(0:100, 9, replace = TRUE)
    #)
    #write.csv(SampleSheet, paste0(system.file(package =
        #"EpigeneticAgePipeline"),"/data/Sample_Sheet.csv"))
    #write.csv(df, paste0(system.file(package = "EpigeneticAgePipeline"),
        #"/data/betaValues.csv"))
#}

#removeTestData <- function() {
    #filePath <- paste0(system.file(package = "EpigeneticAgePipeline"),
        #"/data/betaValues.csv")
    #filePath2 <- paste0(system.file(package = "EpigeneticAgePipeline"),
        #"/data/Sample_Sheet.csv")
    #if (file.exists(filePath)) {
        #unlink(filePath)
        #unlink(filePath2)
        #message("Test data deleted")
    #} else {
        #message("Test data not yet loaded")
    #}
#}




