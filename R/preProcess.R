filterByBeadCount <- function(rgSet,
                              minBeads = 3,
                              sampleMaxFailFrac = 0.05,
                              probeMaxFailFrac  = 0.05,
                              outDir = myEnv$baseDirectory) {
    nb <- tryCatch(minfi::getNBeads(rgSet), error = function(e) NULL)
    if (is.null(nb)) {
        message("Bead counts not available in this object; skipping bead-count filters.")
        return(list(rgSet = rgSet, keepSamples = rep(TRUE, ncol(rgSet)), keepProbes = NULL, nBeads = NULL))
    }

    # observation-level failures
    failObs <- nb < minBeads

    # sample-level: drop if too many low-bead observations
    sampleFailFrac <- colMeans(failObs, na.rm = TRUE)
    keepSamples <- sampleFailFrac <= sampleMaxFailFrac
    if (any(!keepSamples)) {
        failed <- data.frame(
            sample = colnames(nb)[!keepSamples],
            fracLowBead = sampleFailFrac[!keepSamples],
            row.names = NULL
        )
    }
    rgSet2 <- rgSet[, keepSamples, drop = FALSE]
    failObs2 <- failObs[, keepSamples, drop = FALSE]

    # probe-level: drop if too many low-bead observations across remaining samples
    probeFailFrac <- rowMeans(failObs2, na.rm = TRUE)
    keepProbes <- probeFailFrac <= probeMaxFailFrac

    nDroppedProbes <- sum(!keepProbes, na.rm = TRUE)
    message("-----    ","Dropping ", nDroppedProbes,
            " probe(s) failing bead-count threshold (minBeads = ",
            minBeads, ", probeMaxFailFrac = ", probeMaxFailFrac, ").")

    list(rgSet = rgSet2, keepSamples = keepSamples, keepProbes = keepProbes, nBeads = nb)
}

dropWithLog <- function(obj, keep, reason, env = myEnv) {
    if (length(keep) != ncol(obj)) stop("keep length != ncol(obj)")
    failed_now <- colnames(obj)[!keep]
    if (length(failed_now)) {
        message("-----    Dropping ", length(failed_now),
                " sample(s) [", reason, "]:")
        for (nm in failed_now) message("        ", nm)
        env$failedSamples <- unique(c(env$failedSamples, failed_now))
    }
    obj[, keep, drop = FALSE]
}



processIDAT <- function(directory,
                        arrayType,
                        useSampleSheet,
                        sampleSheetFile = "Sample_Sheet.csv",
                        # bead-count thresholds
                        minBeads = 3,
                        beadSampleMaxFailFrac = 0.05,
                        beadProbeMaxFailFrac  = 0.05) {

    dataDirectory <- directory
    myEnv$rgSet <- minfi::read.metharray.exp(dataDirectory, force = TRUE, extended = TRUE)
    myEnv$origSampleNames <- colnames(myEnv$rgSet)
    myEnv$failedSamples <- character(0)

    if (isTRUE(useSampleSheet)) {
        message("Checking dimensionality against ", sampleSheetFile)
        if (!file.exists(file.path(dataDirectory, sampleSheetFile))) {
            stop("Sample sheet file not found: ", file.path(dataDirectory, sampleSheetFile))
        }
        sampleData <- utils::read.csv(file.path(dataDirectory, sampleSheetFile), header = TRUE, check.names = FALSE)
        message("Number of samples in sample sheet: ", nrow(sampleData))
        idatCount <- length(list.files(path = dataDirectory, pattern = "\\.idat$", ignore.case = TRUE)) / 2
        message("Number of samples from IDAT Files: ", idatCount)
        if (idatCount != nrow(sampleData)) {
            message("The number of samples in the sample sheet is not equal to the number of IDAT pairs.")
            stop("Sample sheet / IDAT count mismatch.")
        }
    }

    # Keep annotation in sync with selected array
    if (arrayType == "MSA" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylationMSA" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilm10a1.hg38")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylationMSA", annotation = "ilm10a1.hg38")
    }
    if (arrayType == "EPICv2" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylationEPICv2" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "20a1.hg38")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
    }
    if (arrayType == "EPIC" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylationEPIC" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilm10b4.hg19")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
    }
    if (arrayType == "450K" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylation450k" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilmn12.hg19")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
    }
    if (arrayType == "27K" &&
        (minfi::annotation(myEnv$rgSet)["array"] != "IlluminaHumanMethylation27k" ||
         minfi::annotation(myEnv$rgSet)["annotation"] != "ilmn12.hg19")) {
        minfi::annotation(myEnv$rgSet) <- c(array = "IlluminaHumanMethylation27k", annotation = "ilmn12.hg19")
    }

    ## --- Illumina control metrics ---
    Mset_raw <- minfi::preprocessRaw(myEnv$rgSet)
    mMed <- apply(log2(pmax(minfi::getMeth(Mset_raw), 1)),   2, median, na.rm = TRUE)
    uMed <- apply(log2(pmax(minfi::getUnmeth(Mset_raw), 1)), 2, median, na.rm = TRUE)
    avgMed <- (mMed + uMed)/2
    bad_intensity <- avgMed < 10.5

    if (any(bad_intensity)) {
        bad <- colnames(Mset_raw)[bad_intensity]
        message("Low-intensity samples (minfi getQC cutoff 10.5):")
        for (nm in bad) message(nm)
        keep <- !(colnames(myEnv$rgSet) %in% bad)
        myEnv$rgSet <- dropWithLog(myEnv$rgSet, keep, reason = "low_intensity_minfi")
    }

    # Calculate    the    detection    p-values
    detP <- detectionP(myEnv$rgSet)
    detP <- detP[complete.cases(detP), ]
    samples_before <- ncol(myEnv$rgSet)
    keep <- colMeans(detP) < 0.05
    myEnv$rgSet <- dropWithLog(myEnv$rgSet, keep, reason = "detP_fail")
    samples_removed <- samples_before - ncol(myEnv$rgSet)
    if (samples_removed != 0) {
        message("The following samples were removed due to detP (see log above).")
    }
    message("-----    ", samples_removed, " sample(s) removed due to low detection p")

    # bead-count filtering (sample & probe)
    fb <- filterByBeadCount(
        myEnv$rgSet,
        minBeads = minBeads,
        sampleMaxFailFrac = beadSampleMaxFailFrac,
        probeMaxFailFrac  = beadProbeMaxFailFrac,
        outDir = myEnv$baseDirectory
    )
    myEnv$rgSet <- dropWithLog(myEnv$rgSet, fb$keepSamples, reason = "low_bead_count")
    if (!is.null(fb$keepProbes)) {
        myEnv$rgSet <- myEnv$rgSet[fb$keepProbes, , drop = FALSE]
    }

    # Raw → Methylation preprocessing
    mSetSq <- myEnv$rgSet
    mSetSq <- minfi::preprocessRaw(mSetSq)
    mSetSq <- minfi::mapToGenome(mSetSq)
    mSetSq <- minfi::ratioConvert(mSetSq)

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

    # Remove probes with SNPs at CpG site
    probes_before <- nrow(mSetSqFlt)
    mSetSqFlt <- minfi::dropLociWithSnps(mSetSqFlt)
    probes_removed <- probes_before - nrow(mSetSqFlt)
    message("-----    ", probes_removed, " probe(s) removed for SNPs at CpG.")

    # Cross-reactive probe removal
    probes_before <- nrow(mSetSqFlt)
    if (arrayType == "450K") {
        data("ChenEtAlList", envir = myEnv, package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$ChenEtAlList
    } else if (arrayType == "27K") {
        data("non_specific_probes_Illumina27k", envir = myEnv, package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$non_specific_probes_Illumina27k
    } else if (arrayType == "EPIC") {
        data("PidsleyCrossReactiveProbesEPIC", envir = myEnv, package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$PidsleyCrossReactiveProbesEPIC
    } else if (arrayType == "EPICv2" || arrayType == "MSA") {
        data("epicV2CR", envir = myEnv, package = "EpigeneticAgePipeline")
        xReactiveProbes <- myEnv$epicV2CR
    } else {
        xReactiveProbes <- data.frame(TargetID = character(0))
    }
    keep <- !(sub("_.*", "", minfi::featureNames(mSetSqFlt)) %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep, ]
    probes_removed <- probes_before - nrow(mSetSqFlt)
    message("-----    ", probes_removed, " probe(s) removed for being cross reactive.")

    # Remove sex-chromosome probes
    probes_before <- nrow(mSetSqFlt)
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
    keep <- !(minfi::featureNames(mSetSqFlt) %in% rownames(sexProbes))
    mSetSqFlt <- mSetSqFlt[keep, ]
    probes_removed <- probes_before - nrow(mSetSqFlt)
    message("-----    ", probes_removed, " probe(s) removed for being on sex chromosomes.")

    # Final count
    message("-----    ", nrow(mSetSqFlt), " probe(s) remaining for analysis.")

    # Compute betas
    myEnv$bVals <- minfi::getBeta(mSetSqFlt)
    myEnv$bVals <- myEnv$bVals[, order(colnames(myEnv$bVals)), drop = FALSE]

    # collapsing replicate probes/converting to legacy format for epicv2/msa
    if (arrayType %in% c("EPICv2","MSA")) {
        myEnv$bVals <- sesame::betasCollapseToPfx(myEnv$bVals)
    }
}
