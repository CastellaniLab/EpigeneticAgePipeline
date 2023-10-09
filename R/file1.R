main <- function(directory = getwd(), normalize = TRUE, useBeta = FALSE, arrayType = "450K")
{
  setwd(directory)
  if (!file.exists("output.txt")) {
    file.create("output.txt")
  }

  if (normalize == TRUE) shouldNormalize <- FALSE else shouldNormalize <- TRUE
  bVals <- 0
  rgSet <- 0
  listofCors <- c()
  corsToRemove <- c()
  load(paste0(directory, "/PC-clocks.rda"), envir = .GlobalEnv)
  load(paste0(directory, "/golden_ref.rda"), envir = .GlobalEnv)
  #function for installing packages
  installPackages <- function()
  {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    if (!requireNamespace("minfi", quietly = TRUE)) {
      remotes::install_bioc("minfi")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
    }

    if (!requireNamespace("readr", quietly = TRUE)) {
      install.packages("readr")
    }

    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      install.packages("ggpubr")
    }

    if (!requireNamespace("methylclock", quietly = TRUE)) {
      BiocManager::install("methylclock")
    }

    if (!requireNamespace("methylclockData", quietly = TRUE)) {
      BiocManager::install("methylclockData")
      methylclockData::get_MethylationDataExample()
    }

    if (!requireNamespace("DunedinPACE", quietly = TRUE)) {
      remotes::install_github("danbelsky/DunedinPACE")  # Use 'remotes' for devtools::install_github()
    }

    if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE)) {
      BiocManager::install("IlluminaHumanMethylationEPICmanifest")
      IlluminaHumanMethylationEPICmanifest::IlluminaHumanMethylationEPICmanifest
    }

    if (!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE)) {
      BiocManager::install("IlluminaHumanMethylation450kmanifest")
      IlluminaHumanMethylation450kmanifest::IlluminaHumanMethylation450kmanifest
    }

    if (!requireNamespace("IlluminaHumanMethylation27kmanifest", quietly = TRUE)) {
      BiocManager::install("IlluminaHumanMethylation27kmanifest")
      IlluminaHumanMethylation27kmanifest::IlluminaHumanMethylation27kmanifest
    }

    if (!requireNamespace("FlowSorted.CordBlood.450k", quietly = TRUE)) {
      BiocManager::install("FlowSorted.CordBlood.450k")
    }

    if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
      BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    }

    if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
      BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    }

    if (!requireNamespace("IlluminaHumanMethylation27kanno.ilmn12.hg19", quietly = TRUE)) {
      BiocManager::install("IlluminaHumanMethylation27kanno.ilmn12.hg19")
    }

    if (!requireNamespace("dnaMethyAge", quietly = TRUE)) {
      remotes::install_github("yiluyucheng/dnaMethyAge")
    }
  }

  installPackages()

  #function for processing IDAT files, returns df of beta values
  processIDAT <- function(){
    dataDirectory <-  directory
    pdata <- minfi::read.metharray.sheet(dataDirectory, pattern="Sample_Sheet.csv")
    pdata$Basename <- pdata$Sample_Name

    #Give the samples descriptive names
    pdata$ID <- paste(pdata$Phenotype, pdata$Title,sep=".")

    #Clean up slide names
    pdata$Slide <- gsub("X", "", pdata$Slide)

    write.csv(pdata, file = "Sample_Sheet.csv", row.names = FALSE)

    #Read in the raw data from the IDAT files
    rgSet <- minfi::read.metharray.exp(targets=pdata)
    minfi::sampleNames(rgSet) <- pdata$ID
    #Calculate the detection p-values
    detP <- suppressMessages(minfi::detectionP(rgSet))

    #Track how many samples were removed due to poor quality
    samples_before <- dim(rgSet)[2]

    #Remove poor quality samples
    keep <- colMeans(detP) < 0.05
    rgSet <- rgSet[,keep]

    #Print out number of samples removed due to poor quality
    samples_removed <- samples_before - dim(detP)[2]
    message("----- ", samples_removed, " sample(s) removed due to poor quality") #----- 0 sample(s) removed due to poor quality

    mSetSq <- rgSet

    mSetSq <- minfi::preprocessRaw(mSetSq) #does not perform normalization

    mSetSq <- minfi::mapToGenome(mSetSq)

    mSetSq <- minfi::ratioConvert(mSetSq)


    #Ensure probes are in the same order in the mSetSq and detP objects
    detP <- detP[match(minfi::featureNames(mSetSq),rownames(detP)),]

    #Remove any probes that have failed in one or more samples
    probes_before <- dim(mSetSq)[1]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    mSetSqFlt <- mSetSq[keep,]

    #Print out number of probes removed for failing in one or more samples
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for failing in one or more samples") #----- 5482 probe(s) removed for failing in one or more samples
    probes_before <- dim(mSetSqFlt)[1]

    #Remove probes with SNPs at CpG site
    mSetSqFlt <- minfi::dropLociWithSnps(mSetSqFlt)

    #Print out number of probes removed for having SNPs at CpG site
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for having SNPs at CpG site") #----- 29130 probe(s) removed for having SNPs at CpG site
    probes_before <- dim(mSetSqFlt)[1]

    #Exclude cross reactive probes
    if (arrayType == "450K")
    {
      xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                             "ChenEtAlList.csv",
                                             sep="/"), stringsAsFactors=FALSE)
    }
    else if (arrayType == "27K")
    {
      xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                             "non-specific-probes-Illumina27k.csv",
                                             sep="/"), stringsAsFactors=FALSE)
    }
    else
    {
      xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                             "PidsleyCrossReactiveProbesEPIC.csv",
                                             sep="/"), stringsAsFactors=FALSE)
    }

    keep <- !(minfi::featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep,]

    #Print out number of probes removed for being cross reactive
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for being cross reactive") #----- 25462 probe(s) removed for being cross reactive
    probes_before <- dim(mSetSqFlt)[1]

    #Remove Sex Probes
    if (arrayType == "EPIC")
    {
      ann <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest
    }
    else if(arrayType == "450K")
    {
      ann <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Manifest
    }
    else
    {
      ann <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Manifest
    }

    sexProbes <- ann[which(ann$chr %in% c("chrX", "chrY")),]
    keep <- !(minfi::featureNames(mSetSqFlt) %in% sexProbes$Name)
    mSetSqFlt <- mSetSqFlt[keep,]

    #Print out number of probes removed for being cross reactive
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for being on sex chromosomes") #----- 17936 probe(s) removed for being on sex chromosomes

    #Print out the number of probes remaining
    message("----- ", dim(mSetSqFlt)[1], " probe(s) remaining for analysis") #----- 787849 probe(s) remaining for analysis

    #Calculate methylation beta values
    rgSet <<- rgSet
    bVals <<- minfi::getBeta(mSetSqFlt)
  }

  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor = 2, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="complete.obs"))
    listofCors <<- append(listofCors, r)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor)
  }

  mydiag.panel <- function( x,  labels, ...){
    ll <- par("usr")
    rect(ll[1], ll[3], ll[2], ll[4], col="#CC7178")
  }

  #function for finding the correlation between covariates
  corCovariates = function(x){
    corDf <- pdataSVs
    corDf <- corDf[1:ncol(pdataSVs),]
    row.names(corDf) <- colnames(corDf)
    corDf[,] <- 0

    counter <- 1
    for (i in 1:ncol(corDf))
    {
      for (j in i:(ncol(corDf) - 1)) #rows
      {
        corDf[j + 1, i] <- listofCors[counter]
        counter <- counter + 1
      }
    }

    corDf <- corDf[-nrow(corDf),]
    print(corDf)

    counter <- 1
    for (i in 1:ncol(corDf))
    {
      for (j in i:(ncol(corDf) - 1))
      {
        if (!is.na(corDf[j + 1, i]))
        {
          if (corDf[j + 1, i] > 0.6)
          {
            covariate1 <- rownames(corDf)[j+1]
            covariate2 <- colnames(corDf)[i]
            x <- paste(x, "\n", covariate1, " and ", covariate2, " are highly correlated: ",corDf[j + 1, i],"\n" )
            print(paste("\n", covariate1, " and ", covariate2, " are highly correlated: ",corDf[j + 1, i],"\n"))
            if (!(covariate1 %in% corsToRemove))
            {
              corsToRemove <<- append(corsToRemove, covariate1)
            }
            if (!(covariate2 %in% corsToRemove))
            {
              corsToRemove <<- append(corsToRemove, covariate2)
            }
          }
        }

      }
    }
    return(x)
  }

  #function for plotting correlation line between matrix
  twolines = function(x,y) {
    points(x,y,pch=20)
    abline(lm(y~x),col="#CC7178")
  }

  residGeneration <- function(pdata, output) {
    columns <- colnames(pdata)
    columnsUsed <- columns[1:length(columns)]
    string <- paste(columnsUsed[-1], collapse = " + ")

    rowIndex <- which(columns == "Row")
    columnIndex <- which(columns == "Column")
    slideIndex <- which(columns == "Slide")
    batchIndex <- which(columns == "Batch")

    if (!("Column" %in% columnsUsed))
    {
      formula_string <- paste(columns[1], "~", string, " + ", "(Row|Slide)", " + ", "(1|Batch)")
    } else if (!("Slide" %in% columnsUsed))
    {
      formula_string <- paste(columns[1], "~",  string, " + ", "(Row&Column)", " + ", "(1|Batch)")
    } else
    {
      formula_string <- paste(columns[1], "~", string, " + ", "(1 | Slide)",  " + ",  "(Row + Column | Slide)", " + ",  "(1 | Batch)")
    }

    runlme <- function(formula) {
      lme1 <- glmmTMB::glmmTMB(formula, data = pdata, family = "gaussian", control = glmmTMB::glmmTMBControl(optCtrl=list(iter.max=10000,eval.max=10000)))
      smodel <- lme1
      return(smodel)
    }

    lme_formula <- formula_string
    lme_formula <- as.formula(lme_formula)
    print(lme_formula)

    lme_summary <- try(runlme(lme_formula), silent = FALSE)

    resids <- residuals(lme_summary)
    print(resids)

    for (i in 1:length(resids))
    {
      output <- paste(output, resids[i], "\n")
    }

    return(output)
  }

  removeCovariates <- function(df) {
    print("The following covariates were found to be highly correlated: \n")
    print(corsToRemove)
    print("\nTo remove one of the covariates or several, enter 1 for each index want to remove, 0 to keep")
    userInput <- scan(file = "", n = length(corsToRemove))
    userInput <- as.numeric(userInput)
    print(userInput)
    for (i in 1:length(userInput)){
      if (userInput[i] == 1) {
        column <- corsToRemove[i]
        print(column)
        df <- df[, !grepl(column, names(df))]
      }
      print(i)
    }
    return(df)
  }

  removeOutliers <- function(df,isSampleSheet = FALSE) {
    if (isSampleSheet == FALSE)
    {
      df <- df[!(rownames(df) %in% outlier), ]
      return (df)
    } else {
      df <- df[!(df$ID %in% outlier), ]
      return (df)
    }

  }

  processAgeType <- function(data, age_type, output) {
    pdataSVs$Clock <- as.numeric(data[[age_type]])
    diag.labels[1] <- age_type
    grDevices::cairo_pdf(paste("matrixplot", age_type, ".pdf", sep = ""), width = 14, height = 14, fallback_resolution = 1000)
    print(pairs(plot.formula, data = pdataSVs, upper.panel = twolines, labels = diag.labels, diag.panel = mydiag.panel, lower.panel = panel.cor, label.pos = 0.5, main = ""))
    grDevices::dev.off()

    finalOutput <- paste(output, "\n", age_type, "Covariates\n")
    finalOutput <- corCovariates(finalOutput)
    covariate_data <- pdataSVs

    covariate_data <- removeCovariates(covariate_data)

    finalOutput <- paste(finalOutput, "\n", age_type, "Residuals Based on Epigenetic Age", "\n")
    finalOutput <- residGeneration(pdata = covariate_data, output = finalOutput)
    finalOutput <- paste(finalOutput, "\n", age_type, "Residuals Based on Epigenetic Age Acceleration", "\n")
    covariate_data$Clock <- covariate_data$Clock - covariate_data$Age
    finalOutput <- residGeneration(pdata = covariate_data, output = finalOutput)

    listofCors <<- c()
    corsToRemove <<- c()

    return(finalOutput)
  }

  createGroupedBarChart <- function(data, x, y, fill, title) {
    melted_df <- reshape2::melt(data, id.vars = x, variable.name = fill)
    custom_palette <- c("age" = "red", "horvath" = "#66c2a5", "skinhorvath" = "#3288bd", "hannum" = "#5e4fa2", "levine"="#3288dd")

    plot <- ggplot2::ggplot(data = melted_df, ggplot2::aes_string(x = x, y = y, fill = fill)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::labs(x = x, y = "Age", title = title) +
      ggplot2::scale_fill_manual(values = custom_palette) +  # Apply the custom palette
      ggplot2::theme_minimal()

    plot <- plot + ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"))

    ggplot2::ggsave("SampleIDandAge.png", plot = plot, width = 10000, height = 1000, units = "px", dpi = 300)

    for (i in 2:(ncol(data) - 1))
    {
      plot <- ggpubr::ggscatter(data, x = "age", y = colnames(data)[i],
                                add = "reg.line",  #
                                add.params = list(color = "blue", fill = "lightgray"))
      plot <- plot + ggpubr::stat_cor(method = "pearson")
      ggplot2::ggsave(filename = paste0("plot_", colnames(data)[i], ".png"),
                      plot = plot,
                      width = 1000, height = 1000,
                      units = "px")
    }

  }

  if (file.exists("betaValues.csv")){
    if (useBeta == TRUE)
    {
      bVals <- read.csv("betaValues.csv")
    }
    else
    {
      processIDAT()
    }
  } else {
    processIDAT()
  }


  #PCA Generation ####

  bVals <- na.omit(bVals)

  bValst <- NULL
  if (useBeta == TRUE)
  {
    bValst <- t(bVals[,-1])
  } else
  {
    bValst <- t(bVals)
  }

  bpca <- prcomp(bValst, center=TRUE, scale=FALSE)

  pca_scores <- as.data.frame(bpca$x)

  constant = 3

  sample_outliers=c()
  alloutliers=c()
  for(i in 1:5){
    a<-subset(rownames(bpca$x), bpca$x[,i] > (mean(bpca$x[,i])+constant*sd(bpca$x[,i])))
    b<-subset(rownames(bpca$x), bpca$x[,i] < (mean(bpca$x[,i])-constant*sd(bpca$x[,i])))
    out<-c(a,b)
    sample_outliers <- c(sample_outliers,out)
    print(paste("outliers in PCA",i,":",sep=""))
    print(sample_outliers)
    alloutliers=c(alloutliers,sample_outliers)
    sample_outliers=c()
  }
  outlier<-unique(alloutliers)
  bVals <- bVals[,!( colnames(bVals) %in% outlier)]
  if (useBeta == TRUE)
  {
    bValst <- t(bVals[,-1])
  } else
  {
    bValst <- t(bVals)
  }
  bpca <- prcomp(bValst, center=TRUE, scale=FALSE)
  pca_scores <- as.data.frame(bpca$x)


  #determining cell composition
  if (class(rgSet) != "numeric" & arrayType != "27K")
  {
    FlowSorted.CordBlood.450k::FlowSorted.CordBlood.450k

    CC <- minfi::estimateCellCounts(rgSet, compositeCellType = "CordBlood",
                                    processMethod = "auto", probeSelect = "auto",
                                    cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC"),
                                    referencePlatform = c("IlluminaHumanMethylation450k"),
                                    returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)

    CC <- removeOutliers(CC)
  }

  #Correlation Matrix Construction ####

  sampleData <- read.csv("Sample_Sheet.csv")
  sampleData <- removeOutliers(sampleData, TRUE)
  pdataSVs <- data.frame(Clock = 1:nrow(sampleData), Age = 1:nrow(sampleData), Sex = 1:nrow(sampleData), Smoking_Status = 1:nrow(sampleData), Batch = 1:nrow(sampleData), Slide = 1:nrow(sampleData), Row = 1:nrow(sampleData), Column = 1:nrow(sampleData),
                         Bcell = 1:nrow(sampleData), CD4T = 1:nrow(sampleData), CD8T = 1:nrow(sampleData), Gran = 1:nrow(sampleData), Mono = 1:nrow(sampleData), nRBC = 1:nrow(sampleData),
                         P1 = 1:nrow(sampleData), P2 = 1:nrow(sampleData), P3 = 1:nrow(sampleData), P4 = 1:nrow(sampleData), P5 = 1:nrow(sampleData))
  pdataSVs[,-1] <- NA
  if ("Array" %in% names(sampleData))
  {
    row <- gsub("R(\\d+).*", "\\1", sampleData$Array)
    column <- gsub(".*C(\\d+)", "\\1", sampleData$Array)
    pdataSVs$Row <- as.factor(row)
    pdataSVs$Column <- as.factor(column)
  }

  if ("Age" %in% names(sampleData)) pdataSVs$Age <- as.numeric(sampleData$Age)
  if ("Sex" %in% names(sampleData)) pdataSVs$Sex <- as.factor(sampleData$Sex)
  if ("Smoking_Status" %in% names(sampleData)) pdataSVs$Smoking_Status <- as.factor(sampleData$Smoking_Status)
  if ("Batch" %in% names(sampleData)) pdataSVs$Batch <- as.factor(sampleData$Batch)
  if ("Slide" %in% names(sampleData)) pdataSVs$Slide <- as.factor(sampleData$Slide)
  if (class(rgSet) != "numeric" & arrayType != "27K")
  {
    pdataSVs$Bcell <- as.numeric(CC[, "Bcell"])
    pdataSVs$CD4T <- as.numeric(CC[, "CD4T"])
    pdataSVs$CD8T <- as.numeric(CC[, "CD8T"])
    pdataSVs$Gran <- as.numeric(CC[, "Gran"])
    pdataSVs$Mono <- as.numeric(CC[, "Mono"])
    pdataSVs$nRBC <- as.numeric(CC[, "nRBC"])
  } else
  {
    if ("Bcell" %in% names(sampleData)) pdataSVs$Bcell <- as.numeric(sampleData$Bcell)
    if ("CD4T" %in% names(sampleData)) pdataSVs$CD4T <- as.numeric(sampleData$CD4T)
    if ("CD8T" %in% names(sampleData)) pdataSVs$CD8T <- as.numeric(sampleData$CD8T)
    if ("Gran" %in% names(sampleData)) pdataSVs$Gran <- as.numeric(sampleData$Gran)
    if ("Mono" %in% names(sampleData)) pdataSVs$Mono <- as.numeric(sampleData$Mono)
    if ("nRBC" %in% names(sampleData)) pdataSVs$nRBC <- as.numeric(sampleData$nRBC)
  }
  if ("PC1" %in% names(pca_scores)) pdataSVs$P1 <- as.numeric(pca_scores$PC1)
  if ("PC2" %in% names(pca_scores)) pdataSVs$P2 <- as.numeric(pca_scores$PC2)
  if ("PC3" %in% names(pca_scores)) pdataSVs$P3 <- as.numeric(pca_scores$PC3)
  if ("PC4" %in% names(pca_scores)) pdataSVs$P4 <- as.numeric(pca_scores$PC4)
  if ("PC5" %in% names(pca_scores)) pdataSVs$P5 <- as.numeric(pca_scores$PC5)


  pdataSVs <- pdataSVs[, sapply(pdataSVs, function(col) length(unique(col)) >= 2)]
  print(pdataSVs)
  print(sampleData)

  diag.labels=colnames(pdataSVs)
  pdataColumns <- names(pdataSVs)
  plot.formula=as.formula(paste("~", paste(pdataColumns, collapse = "+")))


  #Clock Output####

  results <- NULL

  if (shouldNormalize == TRUE) results <- methylclock::DNAmAge(bVals, normalize = TRUE, age = sampleData$Age) else results <- methylclock::DNAmAge(bVals, normalize = FALSE, age = sampleData$Age)
  print(results)
  print(results$Hannum)

  finalOutput <- "Raw Clock Results\n"
  finalOutput <- paste(finalOutput, "SampleID", "\t", "Horvath", "\t", "        SkinHorvath", "\t", "        Hannum", "\t", "        PhenoAge","\n" )
  for (i in 1:nrow(results)) {
    finalOutput <- paste(finalOutput, results[i, 1], "\t", results$Horvath[i], "\t", results$skinHorvath[i], "\t", results$Hannum[i], "\t", results$Levine[i], "\n")
  }

  finalOutput <- processAgeType(results, "Horvath", finalOutput)
  finalOutput <- processAgeType(results, "skinHorvath", finalOutput)
  finalOutput <- processAgeType(results, "Hannum", finalOutput)
  finalOutput <- processAgeType(results, "Levine", finalOutput)

  plotDf <- data.frame(
    sample = sampleData$ID,
    horvath = results$Horvath,
    skinhorvath = results$skinHorvath,
    hannum = results$Hannum,
    levine = results$Levine,
    age = sampleData$Age
  )


  exportDf <- results[,c("id", "Horvath", "Hannum", "Levine", "skinHorvath","age")]

  createGroupedBarChart(plotDf, "sample", "value", "Age_Measure", "Sample ID and Type of Age Measure")

  finalOutput <- paste(finalOutput, "\n\nDunedinPACE\n")

  if (!is.numeric(bVals[5,1])) {
    rownames(bVals) <- bVals[,1]
    bVals <- bVals[,-1]
  }

  results <- DunedinPACE::PACEProjector(as.matrix(bVals), proportionOfProbesRequired = 0.8)
  results <- as.data.frame(results)
  for (i in 1:nrow(results))
  {
    finalOutput <- paste(finalOutput, rownames(results)[i], "\t", results$DunedinPACE[i], "\n")
  }

  results$id <- rownames(results)

  exportDf <- merge(exportDf, results, by="id")

  finalOutput <- paste(finalOutput, "\n\nGrimAGE\n")
  grimDf <- data.frame(
    Sample = sampleData$ID,
    Age = sampleData$Age,
    Sex = sampleData$Sex
  )

  grimDf$Sex <- ifelse(grimDf$Sex == "M", "Male", "Female")
  clockname <- "PCGrimAge"
  grimage <- dnaMethyAge::methyAge(beta = bVals, clock = clockname, age_info = as.data.frame(grimDf), do_plot = FALSE)
  grimage$id <- grimage$Sample
  exportDf <- merge(exportDf, grimage, by="id")
  for (i in 1:nrow(results))
  {
    finalOutput <- paste(finalOutput, grimage[i,1], "\t", grimage$Age_Acceleration[i], "\n")
  }

  write.table(as.data.frame(exportDf), file = "epigeneticAge.txt")

  outputString <- paste(finalOutput, collapse = "\n")

  readr::write_file(outputString, file = "output.txt")

}


main(directory = "C:/Users/stanl/Desktop/Development/R/DNAm-age-pipeline/data", useBeta = TRUE, arrayType = "450K")






