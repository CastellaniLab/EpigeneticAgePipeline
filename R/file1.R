    utils::globalVariables(c("bVals", "rgSet", "listofCors", "corsToRemove",
                            "methyAge", "pdataSVs", "exportDf"))

    # Main function starts here
    main <- function(directory = getwd(),
                    normalize = TRUE,
                    useBeta = FALSE,
                    arrayType = "450K",
                    generateResiduals = TRUE,
                    useSampleSheet = TRUE) {
        setwd(directory)
        if (!file.exists("output.txt")) {
            file.create("output.txt")
        }

        # Checking the condition for normalization
        if (normalize == TRUE) {
            shouldNormalize <- FALSE
        } else {
            shouldNormalize <- TRUE
        }

        installDirectory <- paste0(path.package("EpigeneticAgePipeline"),
                                                "/extdata/")

        devtools::load_all(paste0(installDirectory, "dnaMethyAge-main"))

        if (!exists("bVals")) {
            base::assign("bVals", 0, envir = .GlobalEnv)
        }

        # More assignments to global environment variables...

        base::assign("rgSet", 0, envir = .GlobalEnv)
        base::assign("listofCors", c(), envir = .GlobalEnv)
        base::assign("corsToRemove", c(), envir = .GlobalEnv)
        base::assign("pdataSVs", 0, envir = .GlobalEnv)
        base::assign("exportDf", 0, envir = .GlobalEnv)

        # Loading necessary data files
        load(paste0(installDirectory, "PC-clocks.rda"), envir = .GlobalEnv)
        load(paste0(installDirectory, "golden_ref.rda"), envir = .GlobalEnv)

        # Definition of processIDAT function starts here
        processIDAT <- function() {
            dataDirectory <- directory
            containsSampleNames <- FALSE
            if (useSampleSheet == TRUE) {
                testDf <- read.csv("Sample_Sheet.csv", header = TRUE)
                if ("Sample_Name" %in% colnames(testDf))
                {
                    containsSampleNames <- TRUE
                }
            }

            if (containsSampleNames == TRUE) {
                pdata <- minfi::read.metharray.sheet(dataDirectory,
                    pattern = "Sample_Sheet.csv"
                )

                pdata$Basename <- pdata$Sample_Name

                pdata$ID <- paste(pdata$Phenotype, pdata$Title, sep = ".")

                pdata$Slide <- gsub("X", "", pdata$Slide)

                write.csv(pdata, file = "Sample_Sheet.csv", row.names = FALSE)


                rgSet <-
                    minfi::read.metharray.exp(targets = pdata, force = TRUE)
                minfi::sampleNames(rgSet) <- pdata$ID
            } else {
                rgSet <- minfi::read.metharray.exp(dataDirectory, force = TRUE)
            }


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

        # Definition of function used to find highly correlated covariates
        corCovariates <- function(x) {
            corDf <- pdataSVs

            corDf <- corDf[seq.default(from = 1, to = (ncol(pdataSVs))), ]

            row.names(corDf) <- colnames(corDf)

            corDf[, ] <- 0

            counter <- 1
            for (i in seq.default(from = 1, to = (ncol(corDf))))
            {
                for (j in seq.default(from = i,
                to = (ncol(corDf)) - 1))
                {
                    corDf[j + 1, i] <- listofCors[counter]
                    counter <- counter + 1
                }
            }


            corDf <- corDf[-nrow(corDf), ]


            counter <- 1
            for (i in seq.default(from = 1, to = (ncol(corDf))))
            {
                for (j in seq.default(from = i, to = (ncol(corDf)) - 1))
                {
                    if (!is.na(corDf[j + 1, i])) {
                        if (corDf[j + 1, i] > 0.6) {
                            covariate1 <- rownames(corDf)[j + 1]
                            covariate2 <- colnames(corDf)[i]
                            x <- paste(
                                x,
                                "\n",
                                covariate1,
                                " and ",
                                covariate2,
                                " are highly correlated: ",
                                corDf[j + 1, i],
                                "\n"
                            )
                            message("\n")
                            message(
                                covariate1,
                                " and ",
                                covariate2,
                                " are highly correlated: ",
                                corDf[j + 1, i]
                            )
                            message("\n")
                            if (!(covariate1 %in% corsToRemove)) {
                                .GlobalEnv$corsToRemove <- append(corsToRemove,
                                                                    covariate1)
                            }
                            if (!(covariate2 %in% corsToRemove)) {
                                .GlobalEnv$corsToRemove <- append(corsToRemove,
                                                                    covariate2)
                            }
                        }
                    }
                }
            }
            return(x)
        }

        # Definition for function used in matrix generation
        twolines <- function(x, y) {
            points(x, y, pch = 20)
            abline(lm(y ~ x), col = "#CC7178")
        }

        # Definition for function used for residual generation
        residGeneration <- function(pdata, output) {
            columns <- colnames(pdata)
            columnsUsed <- columns[seq.default(from = 1, to = length(columns))]
            string <- paste(columnsUsed[-1], collapse = "    +    ")

            rowIndex <- which(columns == "Row")
            columnIndex <- which(columns == "Column")
            slideIndex <- which(columns == "Slide")
            batchIndex <- which(columns == "Batch")

            if (!("Column" %in% columnsUsed)) {
                formula_string <- paste(
                    columns[1],
                    "~",
                    string,
                    "    +    ",
                    "(Row|Slide)",
                    "    +    ",
                    "(1|Batch)"
                )
            } else if (!("Slide" %in% columnsUsed)) {
                formula_string <- paste(
                    columns[1],
                    "~",
                    string,
                    "    +    ",
                    "(Row&Column)",
                    "    +    ",
                    "(1|Batch)"
                )
            } else {
                formula_string <- paste(
                    columns[1],
                    "~",
                    string,
                    "    +    ",
                    "(1    |    Slide)",
                    "    +    ",
                    "(Row    +    Column    |    Slide)",
                    "    +    ",
                    "(1    |    Batch)"
                )
            }

            if (!("Row" %in% columnsUsed) | !("Batch" %in% columnsUsed)) {
                formula_string <- paste(columns[1], "~", string)
            }

            runlme <- function(formula) {
                lme1 <- glmmTMB::glmmTMB(formula,
                    data = pdata,
                    family = "gaussian",
                    control = glmmTMB::glmmTMBControl(optCtrl = list(
                        iter.max    =    10000,
                        eval.max    =    10000
                    ))
                )
                smodel <- lme1
                return(smodel)
            }

            lme_formula <- formula_string
            lme_formula <- as.formula(lme_formula)
            message(lme_formula)

            lme_summary <- try(runlme(lme_formula), silent = FALSE)

            resids <- residuals(lme_summary)

            for (i in seq.default(from = 1, to = length(resids)))
            {
                output <- paste(output,
                                rownames(pdata)[i],
                                "\t",
                                resids[i],
                                "\n")
            }

            return(output)
        }

        # Definition for function used for user specified covariate removal
        removeCovariates <- function(df) {
            message("The following covariates
                    were found to be highly correlated: \n")
            for (i in corsToRemove)
            {
                message(i)
            }
            message("\nTo remove",
            "one of the covariates or several,",
            "enter 1 for",
            "each index  want to  remove,",
            "0   to  keep")
            userInput <- scan(file = "", n = length(corsToRemove))
            userInput <- as.numeric(userInput)
            for (i in seq.default(from = 1, to = length(userInput))) {
                if (userInput[i] == 1) {
                    column <- corsToRemove[i]
                    message(column)
                    df <- df[, !grepl(column, names(df))]
                }
                message(i)
            }
            return(df)
        }

        # Definition for function used for removing outlier samples
        removeOutliers <- function(df, isSampleSheet = FALSE) {
            if (isSampleSheet == FALSE) {
                df <- df[!(rownames(df) %in% outlier), ]
                return(df)
            } else {
                df <- df[!(df$ID %in% outlier), ]
                return(df)
            }
        }

        # Definition for function used for processing epigenetic age measures
        processAgeType <- function(data, age_type, output) {
            .GlobalEnv$pdataSVs$Clock <- as.numeric(data[[age_type]])
            for (i in seq.default(from = 1, to = ncol(pdataSVs)))
            {
                if(length(unique(pdataSVs[,i])) == 1)
                {
                    .GlobalEnv$pdataSVs[,i] <-
                        jitter(pdataSVs[,i], factor = 0.00001, amount = 0.1)
                }
            }
            diag.labels[1] <- age_type
            grDevices::cairo_pdf(
                paste("matrixplot",
                    age_type,
                    ".pdf",
                    sep = ""
                ),
                width = 14,
                height = 14,
                fallback_resolution = 1000
            )
            pairs(plot.formula,
                data = pdataSVs,
                upper.panel = twolines,
                labels = diag.labels,
                diag.panel = mydiag.panel,
                lower.panel = panel.cor,
                label.pos = 0.5,
                main = ""
            )
            grDevices::dev.off()

            finalOutput <- paste(output, "\n", age_type, "Covariates\n")
            finalOutput <- corCovariates(finalOutput)
            covariate_data <- pdataSVs


            if (generateResiduals == TRUE) {
                covariate_data <- removeCovariates(covariate_data)
                finalOutput <- paste(
                    finalOutput,
                    "\n",
                    age_type,
                    "Residuals Based on Epigenetic Age",
                    "\n"
                )
                finalOutput <- residGeneration(
                    pdata = covariate_data,
                    output = finalOutput
                )
                if ("Age" %in% colnames(pdataSVs)) {
                    finalOutput <- paste(
                        finalOutput,
                        "\n",
                        age_type,
                        "Residuals Based on Epigenetic Age Acceleration",
                        "\n"
                    )
                    covariate_data$Clock <-
                        covariate_data$Clock - covariate_data$Age
                    finalOutput <- residGeneration(
                        pdata = covariate_data,
                        output = finalOutput
                    )
                }
            }
            .GlobalEnv$listofCors <- c()
            .GlobalEnv$corsToRemove <- c()
            return(finalOutput)
        }

        #Definition for function used for generating the grouped bar chart
        createGroupedBarChart <- function(data, x, y, fill, title) {
            melted_df <- reshape2::melt(data, id.vars = x, variable.name = fill)
            custom_palette <- c(
                "age" = "red",
                "horvath" = "#66c2a5",
                "skinhorvath" = "#3288bd",
                "hannum" = "#5e4fa2",
                "levine" = "#3288dd"
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
                ggplot2::theme(plot.background =
                    ggplot2::element_rect(fill = "white"))

            ggplot2::ggsave("SampleIDandAge.png",
                plot = plot,
                width = 10000,
                height = 1000,
                units = "px",
                dpi = 300
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


        if (useBeta == TRUE) {
            if (typeof(bVals) != "list") {
                bVals <- read.csv("betaValues.csv")
            }
        } else {
            processIDAT()
        }


        #    PCA    Generation    ####

        bVals <- na.omit(bVals)

        bValst <- NULL
        if (useBeta == TRUE) {
            bValst <- t(bVals[, -1])
        } else {
            bValst <- t(bVals)
        }

        bpca <- prcomp(bValst, center = TRUE, scale = FALSE)

        pca_scores <- as.data.frame(bpca$x)

        constant <- 3

        sample_outliers <- c()
        alloutliers <- c()
        if (ncol(pca_scores) < 5) {
            loopNum <- ncol(pca_scores)
        }
        else {
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
        outlier <- unique(alloutliers)
        outlierIndexes <- c()
        for (i in seq.default(from = 1, to = ncol(bVals)))
        {
            if (colnames(bVals)[i] %in% outlier) {
                outlierIndexes <- append(outlierIndexes, i)
            }
        }
        bVals <- bVals[, !(colnames(bVals) %in% outlier)]
        if (useBeta == TRUE) {
            bValst <- t(bVals[, -1])
        } else {
            bValst <- t(bVals)
        }
        bpca <- prcomp(bValst, center = TRUE, scale = FALSE)
        pca_scores <- as.data.frame(bpca$x)


        #    determining    cell    composition
        if (!is.numeric(rgSet) & arrayType != "27K") {
            FlowSorted.CordBlood.450k::FlowSorted.CordBlood.450k

            CC <- minfi::estimateCellCounts(rgSet,
                compositeCellType = "CordBlood",
                processMethod = "auto", probeSelect = "auto",
                cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC"),
                referencePlatform = c("IlluminaHumanMethylation450k"),
                returnAll = FALSE, meanPlot = FALSE, verbose = TRUE
            )

            CC <- removeOutliers(CC)
        }

        write.csv(bVals, file = "extractedBetaValues.csv")

        #    Correlation    Matrix    Construction    ####

        if (useSampleSheet == TRUE) {
            sampleData <- read.csv("Sample_Sheet.csv", header = TRUE)
            if ("ID" %in% colnames(sampleData)) {
                sampleData <- removeOutliers(sampleData, TRUE)
            } else {
                if (colnames(bVals)[1] == "X") {
                    if (!length(outlierIndexes) == 0) {
                        sampleData <-
                            sampleData[-(outlierIndexes - 1), ]
                    }
                } else {
                    if (!length(outlierIndexes) == 0) {
                        sampleData <-
                            sampleData[-outlierIndexes, ]
                    }
                }
            }
            if (colnames(bVals)[1] == "X") {
                .GlobalEnv$pdataSVs <- data.frame(Clock = seq.default(
                    from = 1,
                    to = ncol(bVals) - 1
                ))
                rownames(.GlobalEnv$pdataSVs) <- colnames(bVals)[-1]
            } else {
                .GlobalEnv$pdataSVs <- data.frame(Clock = seq.default(from = 1,
                to = ncol(bVals)))
                rownames(.GlobalEnv$pdataSVs) <- colnames(bVals)
            }

            sampleData <- as.data.frame(sampleData)
            columnData <- read.csv("Sample_Sheet.csv", header = TRUE)
            colnames(sampleData) <- colnames(columnData)

            for (i in colnames(sampleData))
            {
                switch(i,
                    Age = .GlobalEnv$pdataSVs$Age <- as.numeric(sampleData$Age),
                    Sex = .GlobalEnv$pdataSVs$Sex <- as.factor(sampleData$Sex),
                    Smoking_Status =
                        .GlobalEnv$pdataSVs$Smoking_Status <-
                        as.factor(sampleData$Smoking_Status),
                    Batch =
                    .GlobalEnv$pdataSVs$Batch <- as.factor(sampleData$Batch),
                    Slide =
                    .GlobalEnv$pdataSVs$Slide <- as.factor(sampleData$Slide),
                    Bcell =
                    .GlobalEnv$pdataSVs$Bcell <- as.numeric(sampleData$Bcell),
                    CD4T =
                    .GlobalEnv$pdataSVs$CD4T <- as.numeric(sampleData$CD4T),
                    CD8T =
                    .GlobalEnv$pdataSVs$CD8T <- as.numeric(sampleData$CD8T),
                    Gran =
                    .GlobalEnv$pdataSVs$Gran <- as.numeric(sampleData$Gran),
                    Mono =
                    .GlobalEnv$pdataSVs$Mono <- as.numeric(sampleData$Mono),
                    nRBC =
                    .GlobalEnv$pdataSVs$nRBC <- as.numeric(sampleData$nRBC),
                    Array = {
                        row <- as.factor(gsub("R(\\d+).*",
                                                "\\1",
                                                sampleData$Array))
                        column <- as.factor(gsub(".*C(\\d+)",
                                                "\\1",
                                                sampleData$Array))
                        .GlobalEnv$pdataSVs$Row <- row
                        .GlobalEnv$pdataSVs$Column <- column
                    },
                    {
                        message("Looks like you have a custom covariate")
                        message(i)
                        message("Enter 0 if this is a numerical variable,",
                            "or 1 if this is a factor,",
                            "or 2 to ignore")
                        userInput <- scan(file = "", nmax = 1)
                        message(userInput)
                        message(class(userInput))
                        if (userInput == 0) {
                            .GlobalEnv$pdataSVs[[i]] <-
                                as.numeric(sampleData[[i]])
                        } else if (userInput == 1) {
                            .GlobalEnv$pdataSVs[[i]] <-
                                as.factor(sampleData[[i]])
                        }
                    }
                )
            }

            if (!is.numeric(rgSet) & arrayType != "27K") {
                .GlobalEnv$pdataSVs$Bcell <- as.numeric(CC[, "Bcell"])
                .GlobalEnv$pdataSVs$CD4T <- as.numeric(CC[, "CD4T"])
                .GlobalEnv$pdataSVs$CD8T <- as.numeric(CC[, "CD8T"])
                .GlobalEnv$pdataSVs$Gran <- as.numeric(CC[, "Gran"])
                .GlobalEnv$pdataSVs$Mono <- as.numeric(CC[, "Mono"])
                .GlobalEnv$pdataSVs$nRBC <- as.numeric(CC[, "nRBC"])
            }
            if ("PC1" %in% names(pca_scores)) {
                .GlobalEnv$pdataSVs$P1 <- as.numeric(pca_scores$PC1)
            }
            if ("PC2" %in% names(pca_scores)) {
                .GlobalEnv$pdataSVs$P2 <- as.numeric(pca_scores$PC2)
            }
            if ("PC3" %in% names(pca_scores)) {
                .GlobalEnv$pdataSVs$P3 <- as.numeric(pca_scores$PC3)
            }
            if ("PC4" %in% names(pca_scores)) {
                .GlobalEnv$pdataSVs$P4 <- as.numeric(pca_scores$PC4)
            }
            if ("PC5" %in% names(pca_scores)) {
                .GlobalEnv$pdataSVs$P5 <- as.numeric(pca_scores$PC5)
            }


            .GlobalEnv$pdataSVs <- pdataSVs[
                ,
                vapply(
                    pdataSVs,
                    function(col) length(unique(col)) >= 2,
                    logical(1)
                )
            ]
            diag.labels <- colnames(pdataSVs)
            pdataColumns <- names(pdataSVs)
            plot.formula <- as.formula(paste("~",
                                            paste(pdataColumns,
                                            collapse = "+")))
        } else {
            if (colnames(bVals)[1] == "X") {
                .GlobalEnv$pdataSVs <- data.frame(Clock = seq.default(
                    from = 1,
                    to = ncol(bVals) - 1
                ))
                rownames(.GlobalEnv$pdataSVs) <- colnames(bVals)[-1]
            } else {
                .GlobalEnv$pdataSVs <- data.frame(Clock = seq.default(from = 1,
                to = ncol(bVals)))
                rownames(.GlobalEnv$pdataSVs) <- colnames(bVals)
            }
            if (!is.numeric(rgSet) & arrayType != "27K") {
                .GlobalEnv$pdataSVs$Bcell <- as.numeric(CC[, "Bcell"])
                .GlobalEnv$pdataSVs$CD4T <- as.numeric(CC[, "CD4T"])
                .GlobalEnv$pdataSVs$CD8T <- as.numeric(CC[, "CD8T"])
                .GlobalEnv$pdataSVs$Gran <- as.numeric(CC[, "Gran"])
                .GlobalEnv$pdataSVs$Mono <- as.numeric(CC[, "Mono"])
                .GlobalEnv$pdataSVs$nRBC <- as.numeric(CC[, "nRBC"])
            }
            if ("PC1" %in% names(pca_scores))
            {
                .GlobalEnv$pdataSVs$P1 <- as.numeric(pca_scores$PC1)
            }

            if ("PC2" %in% names(pca_scores))
            {
                .GlobalEnv$pdataSVs$P2 <- as.numeric(pca_scores$PC2)
            }

            if ("PC3" %in% names(pca_scores))
            {
                .GlobalEnv$pdataSVs$P3 <- as.numeric(pca_scores$PC3)
            }
            if ("PC4" %in% names(pca_scores))
            {
                .GlobalEnv$pdataSVs$P4 <- as.numeric(pca_scores$PC4)
            }

            if ("PC5" %in% names(pca_scores))
            {
                .GlobalEnv$pdataSVs$P5 <- as.numeric(pca_scores$PC5)
            }

            diag.labels <- colnames(pdataSVs)
            pdataColumns <- names(pdataSVs)
            plot.formula <- as.formula(paste("~",
                                                paste(pdataColumns,
                                                collapse = "+")))
        }



        #    Clock    Output####

        results <- NULL
        if ("Age" %in% colnames(pdataSVs)) {
            if (shouldNormalize == TRUE) {
                results <- methylclock::DNAmAge(bVals,
                    normalize = TRUE,
                    age = pdataSVs$Age
                )
            } else {
                results <- methylclock::DNAmAge(bVals,
                    normalize = FALSE,
                    age = pdataSVs$Age
                )
            }
        } else {
            if (shouldNormalize == TRUE) {
                results <- methylclock::DNAmAge(bVals, normalize = TRUE)
            } else {
                results <- methylclock::DNAmAge(bVals, normalize = FALSE)
            }
        }

        finalOutput <- "Raw    Clock    Results\n"

        headers <- c("SampleID", "Horvath", "SkinHorvath", "Hannum", "PhenoAge")
        finalOutput <- sprintf("%-30s\t%-30s\t%-30s\t%-30s\t%-30s\n",
                                headers[1],
                                headers[2],
                                headers[3],
                                headers[4],
                                headers[5])


        for (i in seq_len(nrow(results))) {
            row <- sprintf("%-30s\t%-30s\t%-30s\t%-30s\t%-30s\n",
                            results[i, 1],
                            results$Horvath[i],
                            results$skinHorvath[i],
                            results$Hannum[i],
                            results$Levine[i])
            finalOutput <- paste0(finalOutput, row)
        }


        finalOutput <- processAgeType(results, "Horvath", finalOutput)
        finalOutput <- processAgeType(results, "skinHorvath", finalOutput)
        finalOutput <- processAgeType(results, "Hannum", finalOutput)
        finalOutput <- processAgeType(results, "Levine", finalOutput)

        if ("Age" %in% colnames(pdataSVs)) {
            plotDf <- data.frame(
                sample = sampleData$ID,
                horvath = results$Horvath,
                skinhorvath = results$skinHorvath,
                hannum = results$Hannum,
                levine = results$Levine,
                age = sampleData$Age
            )

            createGroupedBarChart(
                plotDf,
                "sample",
                "value",
                "Age_Measure",
                "Sample    ID    and    Type    of    Age    Measure"
            )
        }

        if ("Age" %in% colnames(pdataSVs)) {
            .GlobalEnv$exportDf <- results[, c(
                "id",
                "Horvath",
                "Hannum",
                "Levine",
                "skinHorvath",
                "age"
            )]
        } else {
            .GlobalEnv$exportDf <- results[,
                                c("id",
                                "Horvath",
                                "Hannum",
                                "Levine",
                                "skinHorvath")]
        }

        finalOutput <- paste(finalOutput, "\n\nDunedinPACE\n")

        if (!is.numeric(bVals[5, 1])) {
            rownames(bVals) <- bVals[, 1]
            bVals <- bVals[, -1]
        }
        clockname <- "DunedinPACE"
        results <- methyAge(
            betas = bVals,
            clock = clockname,
            do_plot = FALSE
        )
        results <- as.data.frame(results)
        for (i in seq.default(from = 1, to = nrow(results)))
        {
            finalOutput <- paste(
                finalOutput,
                results[i, 1],
                "\t",
                results[i, 2],
                "\n"
            )
        }

        results$id <- rownames(results)

        .GlobalEnv$exportDf$DunedinPACE <- results[,2]
        if ("Age" %in% colnames(pdataSVs)) {
            finalOutput <- paste(finalOutput, "\n\nGrimAGE\n")
            grimDf <- data.frame(
                Sample = sampleData$ID,
                Age = sampleData$Age,
                Sex = sampleData$Sex
            )

            grimDf$Sex <- ifelse(grimDf$Sex == "M", "Male", "Female")
            clockname <- "PCGrimAge"
            grimage <- methyAge(
                betas = bVals,
                clock = clockname,
                age_info = as.data.frame(grimDf),
                do_plot = FALSE
            )
            grimage$id <- grimage$Sample
            .GlobalEnv$exportDf$Grimage <- grimage$Age_Acceleration
            for (i in seq.default(from = 1, to = nrow(results)))
            {
                finalOutput <- paste(
                    finalOutput,
                    grimage[i, 1],
                    "\t",
                    grimage$Age_Acceleration[i],
                    "\n"
                )
            }
        }

        .GlobalEnv$exportDf <- as.data.frame(exportDf)
        write.table(as.data.frame(exportDf), file = "epigeneticAge.txt")

        outputString <- paste(finalOutput, collapse = "\n")

        readr::write_file(outputString, file = "output.txt")
    }
