# Function for getting cell counts
estimateCellCountsFlexible <- function(rgSet,
                                       arrayType = c("450K","EPIC","EPICv2","MSA", "27K"),
                                       tissue = c("bloodAdult","bloodCord","saliva",
                                                  "placenta","brainDLPFC", "buccal"),
                                       method = c("CP","RPC"),
                                       placentaTrimester = c("third","first")) {
    arrayType <- match.arg(arrayType)
    tissue <- match.arg(tissue)
    method <- match.arg(method)
    placentaTrimester <- match.arg(placentaTrimester)

    betas <- getBetasForDeconv(rgSet, arrayType)

    # choose reference/compTable
    if (tissue == "bloodAdult") {
        data("FlowSorted.Blood.450k.compTable", package = "FlowSorted.Blood.450k")
        compTable <- myEnv$FlowSorted.Blood.450k.compTable[,c("CD8T","CD4T","NK","Bcell","Mono","Gran")]
    } else if (tissue == "bloodCord") {
        data("FlowSorted.CordBloodCombined.450k.compTable", envir = myEnv,
             package = "FlowSorted.CordBloodCombined.450k")
        compTable <- myEnv$FlowSorted.CordBloodCombined.450k.compTable
    } else if (tissue == "saliva") {
        compTable <- getSalivaCompTable()
    } else if (tissue == "placenta") {
        if (!requireNamespace("planet", quietly = TRUE)) stop("Install 'planet' for placenta.")
        if (placentaTrimester == "third") {
            data("plCellCpGsThird", envir = myEnv, package = "planet")
            compTable <- myEnv$plCellCpGsThird
        } else {
            data("plCellCpGsFirst", envir = myEnv, package = "planet")
            compTable <- myEnv$plCellCpGsFirst
        }
    } else if (tissue == "brainDLPFC") {
        if (!requireNamespace("FlowSorted.DLPFC.450k", quietly = TRUE))
            stop("Install 'FlowSorted.DLPFC.450k' for brain deconvolution.")
        ref <- FlowSorted.DLPFC.450k::FlowSorted.DLPFC.450k
        refBetas <- minfi::getBeta(minfi::preprocessNoob(ref))
        cellType <- as.factor(ref$CellType)
        compTable <- sapply(levels(cellType), function(ct) {
            rowMeans(refBetas[, cellType == ct, drop = FALSE], na.rm = TRUE)
        })
        storage.mode(compTable) <- "double"
    } else if (tissue == "buccal") {
        if (!requireNamespace("EpiDISH", quietly = TRUE))
            stop("Install 'EpiDISH' for buccal deconvolution.")
        # Primary solid-tissue reference: Epi/Fib/Immune
        data("centEpiFibIC.m", envir = myEnv, package = "EpiDISH")
        epiFibIC <- myEnv$centEpiFibIC.m

        # Secondary immune reference (choose per array type)
        if (arrayType %in% c("450K", "27K")) {
            data("cent12CT450k.m", envir = myEnv, package = "EpiDISH")
            immuneRef <- myEnv$cent12CT450k.m
        } else {
            data("cent12CT.m", envir = myEnv, package = "EpiDISH")
            immuneRef <- myEnv$cent12CT.m
        }

        # Coverage checks (both refs)
        pct1 <- round(100 * length(intersect(rownames(betas), rownames(epiFibIC))) / nrow(epiFibIC), 1)
        pct2 <- round(100 * length(intersect(rownames(betas), rownames(immuneRef))) / nrow(immuneRef), 1)
        if (pct1 < 70 || pct2 < 70)
            warning(sprintf("Low buccal reference coverage: EpiFibIC %s%%, ImmuneRef %s%% of CpGs present.",
                            pct1, pct2))

        if (method == "CP") {
            est <- EpiDISH::hepidish(beta.m = betas,
                                     ref1.m = epiFibIC,
                                     ref2.m = immuneRef,
                                     h.CT.idx = 3,
                                     method = "CP")
        } else {
            est <- EpiDISH::hepidish(beta.m = betas,
                                     ref1.m = epiFibIC,
                                     ref2.m = immuneRef,
                                     h.CT.idx = 3,
                                     method = "RPC")
        }
        est <- as.data.frame(est, check.names = FALSE)
        if (is.null(rownames(est))) rownames(est) <- colnames(betas)

        attr(est, "referenceRows") <- c(EpiFibIC = nrow(epiFibIC), ImmuneRef = nrow(immuneRef))
        attr(est, "matchedRows") <- c(EpiFibIC = length(intersect(rownames(betas), rownames(epiFibIC))),
                                      ImmuneRef = length(intersect(rownames(betas), rownames(immuneRef))))
        return(est)
        # -------------------------------------------------------------------------------
    }

    # ensure cg IDs (EPICv2/MSA collapse already done in getBetasForDeconv)
    compTable <- compTable[!duplicated(rownames(compTable)), , drop = FALSE]

    est <- projectWithCompTable(betas, compTable, method = method)

    # coverage + provenance
    attr(est, "referenceRows") <- nrow(compTable)
    attr(est, "matchedRows") <- length(intersect(rownames(betas), rownames(compTable)))
    est
}

getSalivaCompTable <- function() {
    if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
        stop("Install 'ExperimentHub' to retrieve BeadSorted.Saliva.EPIC resources.")
    }
    eh <- ExperimentHub::ExperimentHub()
    resList <- ExperimentHub::loadResources(eh, "BeadSorted.Saliva.EPIC")
    dfIdx <- which(vapply(resList, is.data.frame, logical(1)))
    if (length(dfIdx) == 0) stop("No data.frame resource found in BeadSorted.Saliva.EPIC.")
    df <- resList[[dfIdx[1]]]

    probeCol <- if ("probeName" %in% names(df)) "probeName" else
        if ("CpG" %in% names(df)) "CpG" else
            stop("Neither 'probeName' nor 'CpG' column found in saliva resource.")
    epiCol <- grep("Epithelial", names(df), ignore.case = TRUE, value = TRUE)
    immCol <- grep("\\bImmune\\b", names(df), ignore.case = TRUE, value = TRUE)
    if (length(epiCol) == 0 || length(immCol) == 0) {
        epiCol <- grep("average.*Epithelial", names(df), ignore.case = TRUE, value = TRUE)
        immCol <- grep("average.*Immune", names(df), ignore.case = TRUE, value = TRUE)
    }
    if (length(epiCol) == 0 || length(immCol) == 0) {
        stop("Could not locate epithelial/immune mean columns in saliva resource.")
    }

    mat <- cbind(
        Epithelial = as.numeric(df[[epiCol[1]]]),
        Immune     = as.numeric(df[[immCol[1]]])
    )
    rownames(mat) <- as.character(df[[probeCol]])
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
    storage.mode(mat) <- "double"
    mat
}

getBetasForDeconv <- function(rgSet, arrayType) {
    if (!is.numeric(rgSet)) {
        betas <- minfi::getBeta(minfi::preprocessNoob(rgSet))
        if (arrayType %in% c("EPICv2","MSA")) {
            betas <- sesame::betasCollapseToPfx(betas)
        }
        return(betas)
    }
    return(myEnv$bVals)
}

projectWithCompTable <- function(betas, compTable, method = c("CP","RPC")) {
    method <- match.arg(method)
    common <- intersect(rownames(betas), rownames(compTable))
    pct <- round(100 * length(common) / nrow(compTable), 1)
    if (pct < 70) warning(sprintf("Low reference coverage: %s%% of CpGs present.", pct))

    X <- betas[common, , drop = FALSE]
    W <- compTable[common, , drop = FALSE]

    if (method == "CP") {
        est <- FlowSorted.Blood.EPIC::projectCellType_CP(
            X, W, contrastWBC = NULL, nonnegative = TRUE, lessThanOne = FALSE
        )
        # CP returns samples × cellTypes
        est <- as.data.frame(est, check.names = FALSE)
    } else {
        est <- EpiDISH::epidish(X, ref = W, method = "RPC")$estF
        est <- as.data.frame(est, check.names = FALSE)
    }

    # ensure rownames are sample IDs
    if (is.null(rownames(est))) rownames(est) <- colnames(X)
    est
}

# Adding cell counts to pdataSVs
addCellCountsToPdataSvs <- function(ccMatrix) {
    for (nm in colnames(ccMatrix)) {
        myEnv$pdataSVs[[nm]] <- as.numeric(ccMatrix[, nm])
    }
}
