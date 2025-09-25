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
        if(clock %in% c('PCHorvathS2013', 'PCHorvathS2018')){
            HorvathS2013_transform <- function(x){
                if (x > 0){
                    x <- x *(20 + 1) + 20
                }else{
                    x <- exp(x + log(20 + 1)) - 1
                }
            }
            m_age[,1] <- sapply(m_age[,1], HorvathS2013_transform)
        }
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
                    m_age$PCGrimAge_Acceleration <- getAccel(m_age$Age, m_age$mAge)
                }else{
                    stop(message("\nTo calculate 'PCGrimage', 'age_info'
                        should include a 'Sex' column that contains
                        binary sex annotation, i.e. either Female or Male."))
                }
            }

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
    warning(paste0(back_p, "% of the required probes have been found for DunedinPACE"))
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
