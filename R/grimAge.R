requireCols <- function(df, cols, msg) {
    if (!all(cols %in% names(df))) stop(msg, call. = FALSE)
}

asFemale01 <- function(sexChr) {
    sexChr <- as.character(sexChr)
    sexChr <- ifelse(sexChr %in% c("1","Male","male","M"), "Male", sexChr)
    sexChr <- ifelse(sexChr %in% c("2","Female","female","F"), "Female", sexChr)
    as.numeric(grepl("^F", sexChr, ignore.case = TRUE))
}

coverageReport <- function(coefsDf, betas) {
    coefCpgs <- setdiff(coefsDf$Probe, c("Intercept","Age","Female"))
    present  <- sum(coefCpgs %in% rownames(betas))
    total    <- length(coefCpgs)
    pct      <- if (total == 0) 100 else round(100 * present / total, 2)
    list(present = present, total = total, pct = pct)
}

parseGrimAgeTable <- function(df) {
    # Components (all DNAm blocks)
    colnames(df) <- c("Model", "Probe", "Coefficient")
    dnamRows  <- grepl("^DNAm", df$Model)
    dnamDf    <- df[dnamRows, ]
    dnamSplit <- split(dnamDf[, c("Probe", "Coefficient")], dnamDf$Model)
    components <- lapply(dnamSplit, function(x) {
        x$Probe <- as.character(x$Probe)
        x$Coefficient <- as.numeric(x$Coefficient)
        x
    })

    # COX
    coxDf  <- df[df$Model == "COX", c("Probe", "Coefficient")]
    coxVec <- setNames(as.numeric(coxDf$Coefficient), as.character(coxDf$Probe))

    # Transform
    transformDf   <- df[df$Model == "transform", c("Probe", "Coefficient")]
    transformVals <- setNames(as.numeric(transformDf$Coefficient), as.character(transformDf$Probe))
    cal <- list(
        mAge  = unname(transformVals[["m_age"]]),
        sdAge = unname(transformVals[["sd_age"]]),
        mCox  = unname(transformVals[["m_cox"]]),
        sdCox = unname(transformVals[["sd_cox"]])
    )

    list(components = components, hazard = coxVec, cal = cal)
}

coxToYears <- function(lp, cal) {
    ((lp - cal$mCox) / cal$sdCox) * cal$sdAge + cal$mAge
}

predictSurrogate <- function(coefsDf, betas, ageInfo) {
    stopifnot(all(c("Probe","Coefficient") %in% names(coefsDf)))

    X <- betas

    # Add design rows ONLY if those rows are actually in the coefs
    if ("Intercept" %in% coefsDf$Probe) {
        X <- rbind(X, Intercept = rep(1, ncol(X)))
    }
    if ("Age" %in% coefsDf$Probe) {
        requireCols(ageInfo, c("Sample","Age"),
                    "To compute this surrogate, ageInfo must have Sample and Age.")
        ageMap <- setNames(ageInfo$Age, ageInfo$Sample)
        X <- rbind(X, Age = as.numeric(ageMap[colnames(X)]))
    }
    if ("Female" %in% coefsDf$Probe) {
        requireCols(ageInfo, c("Sample","Sex"),
                    "To compute this surrogate, ageInfo must have Sample and Sex.")
        sexMap <- setNames(asFemale01(ageInfo$Sex), ageInfo$Sample)
        X <- rbind(X, Female = as.numeric(sexMap[colnames(X)]))
    }

    # Build coefficient vector
    coefVec <- setNames(as.numeric(coefsDf$Coefficient), as.character(coefsDf$Probe))

    # SAFE alignment: keep only common rows and order both sides identically
    common <- intersect(rownames(X), names(coefVec))
    if (length(common) == 0L) {
        # no overlap → return all-NA vector (length = samples)
        return(rep(NA_real_, ncol(X)))
    }
    # order by coef order for deterministic math
    common <- intersect(names(coefVec), rownames(X))  # keeps coef ordering
    X <- X[common, , drop = FALSE]
    coefVec <- coefVec[common]

    # guard NAs
    X[is.na(X)] <- 0

    # samples-length numeric vector
    drop(t(X) %*% coefVec)
}

computeGrimAgeV1 <- function(betas, ageInfo,
                             comps, hazard, cal,
                             coverageWarnPct = 85) {
    requireCols(ageInfo, c("Sample","Age","Sex"),
                "GrimAge v1 requires ageInfo with Sample, Age, Sex.")

    covList <- lapply(comps, coverageReport, betas = betas)
    covDf <- data.frame(
        Surrogate = names(covList),
        present   = sapply(covList, `[[`, "present"),
        total     = sapply(covList, `[[`, "total"),
        pct       = sapply(covList, `[[`, "pct"),
        row.names = NULL
    )
    bad <- covDf$pct < coverageWarnPct
    if (any(bad)) {
        warning(sprintf("Low CpG coverage for: %s",
                        paste(covDf$Surrogate[bad], collapse = ", ")))
    }

    needed <- c("DNAmPACKYRS","DNAmADM","DNAmB2M","DNAmCystatinC",
                "DNAmGDF15","DNAmLeptin","DNAmPAI1","DNAmTIMP1")
    missingSurrogates <- setdiff(needed, names(comps))
    if (length(missingSurrogates) > 0) {
        stop("Missing surrogate coefficient blocks: ",
             paste(missingSurrogates, collapse = ", "))
    }

    S <- sapply(needed, function(nm) {
        predictSurrogate(
            coefsDf = comps[[nm]],
            betas   = betas,
            ageInfo = ageInfo
        )
    })
    S <- as.data.frame(S)
    S$Sample <- colnames(betas)

    S <- merge(S, ageInfo[, c("Sample","Age","Sex")], by = "Sample", sort = FALSE)
    S$Female <- asFemale01(S$Sex)

    S$Intercept <- 1
    X  <- as.matrix(S[, names(hazard), drop = FALSE])
    lp <- drop(X %*% hazard)

    years <- coxToYears(lp, cal)

    out <- data.frame(Sample = S$Sample, mAge = years, row.names = NULL)
    attr(out, "coverage") <- covDf
    return(out)
}

computeGrimAgeV2 <- function(betas, ageInfo,
                             comps, hazard, cal,
                             coverageWarnPct = 85) {
    requireCols(ageInfo, c("Sample","Age","Sex"),
                "GrimAge v2 requires ageInfo with Sample, Age, Sex.")

    needed <- c("DNAmPACKYRS","DNAmADM","DNAmB2M","DNAmCystatinC",
                "DNAmGDF15","DNAmLeptin","DNAmPAI1","DNAmTIMP1",
                "DNAmlogCRP","DNAmlogA1C")
    missingSurrogates <- setdiff(needed, names(comps))
    if (length(missingSurrogates) > 0) {
        stop("Missing v2 surrogate blocks: ", paste(missingSurrogates, collapse = ", "))
    }

    covList <- lapply(comps[needed], coverageReport, betas = betas)
    covDf <- data.frame(
        Surrogate = names(covList),
        present   = sapply(covList, `[[`, "present"),
        total     = sapply(covList, `[[`, "total"),
        pct       = sapply(covList, `[[`, "pct"),
        row.names = NULL
    )
    if (any(covDf$pct < coverageWarnPct)) {
        warning(sprintf("Low CpG coverage for: %s",
                        paste(covDf$Surrogate[covDf$pct < coverageWarnPct], collapse = ", ")))
    }

    S <- sapply(needed, function(nm) {
        predictSurrogate(
            coefsDf = comps[[nm]],
            betas   = betas,
            ageInfo = ageInfo
        )
    })
    S <- as.data.frame(S)
    S$Sample <- colnames(betas)

    S <- merge(S, ageInfo[, c("Sample","Age","Sex")], by = "Sample", sort = FALSE)
    S$Female <- asFemale01(S$Sex)

    S$Intercept <- 1
    X  <- as.matrix(S[, names(hazard), drop = FALSE])
    lp <- drop(X %*% hazard)

    years <- coxToYears(lp, cal)

    out <- data.frame(Sample = S$Sample, mAge = years, row.names = NULL)
    attr(out, "coverage") <- covDf
    return(out)
}
