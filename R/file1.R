main <- function(normalize = TRUE, useBeta = FALSE, directory = getwd(), arrayType = "450K")
{
  setwd(directory)
  print(directory)
  start <- Sys.time()
  print(getwd())
  if (!file.exists("output.txt")) {
    file.create("output.txt")
  }
  
  if (normalize == TRUE) shouldNormalize <- FALSE else shouldNormalize <- TRUE
  
  bVals <- NULL
  
  #installing packages 
  
  packages_to_install <- c(
    "BiocManager",
    "minfi",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "tidyverse",
    "ggplot2",
    "ggpubr",
    "umap",
    "IlluminaHumanMethylationEPICmanifest",
    "methylclock",
    "DunedinPACE",
    "FlowSorted.CordBlood.450k",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "lme4",
    "dnaMethyAge"
  )
  
  install <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      message(paste("Installing", package_name))
      if (package_name == "BiocManager") {
        install.packages(package_name)
        BiocManager::install()
      } else {
        install.packages(package_name)
      }
    } else {
      message(paste(package_name, "is already installed"))
    }
  }
  
  for (package in packages_to_install) {
    install(package)
  }
  
  #if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
  #BiocManager::install("minfi")
  #BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  #install.packages(c("tidyverse", "ggplot2", "ggpubr", "umap"))
  #BiocManager::install("IlluminaHumanMethylationEPICmanifest")
  #BiocManager::install("methylclock")
  #devtools::install_github("danbelsky/DunedinPACE")
  #BiocManager::install("FlowSorted.CordBlood.450k")
  #BiocManager::install("IlluminaHumanMethylation450kmanifest")
  #BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  #install.packages("lme4")
  #BiocManager::install("yiluyucheng/dnaMethyAge")
  #BiocManager::install("FlowSorted.CordBlood.450k")
  
  #Load packages
  suppressMessages(library(minfi))
  suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
  suppressMessages(library(IlluminaHumanMethylation450kmanifest))
  suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  suppressMessages(library(methylclock))
  suppressMessages(library("DunedinPACE"))
  suppressMessages(library(FlowSorted.CordBlood.450k))
  suppressMessages(library(tidyverse))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggpubr))
  suppressMessages(library(umap))
  suppressMessages(library(lme4))
  suppressMessages(library('dnaMethyAge'))
  suppressMessages(library(glmmTMB))
  
  
  processIDAT <- function(){
    dataDirectory <-  directory
    pdata <- read.metharray.sheet(dataDirectory, pattern="Sample_Sheet.csv")
    pdata$Basename <- pdata$Sample_Name
    
    #Give the samples descriptive names
    pdata$ID <- paste(pdata$Phenotype, pdata$Title,sep=".")
    
    #Clean up slide names
    pdata$Slide <- gsub("X", "", pdata$Slide)
    
    write.csv(pdata, file = "Sample_Sheet.csv", row.names = FALSE)
    
    #Read in the raw data from the IDAT files
    rgSet <- read.metharray.exp(targets=pdata)
    sampleNames(rgSet) <- pdata$ID
    
    #Calculate the detection p-values
    detP <- suppressMessages(detectionP(rgSet))
    
    #Track how many samples were removed due to poor quality
    samples_before <- dim(rgSet)[2]
    
    #Remove poor quality samples
    keep <- colMeans(detP) < 0.05
    rgSet <<- rgSet[,keep]
    
    #Print out number of samples removed due to poor quality
    samples_removed <- samples_before - dim(detP)[2]
    message("----- ", samples_removed, " sample(s) removed due to poor quality") #----- 0 sample(s) removed due to poor quality
    
    #Normalize the data *** REMOVING THIS STEP *** 
    #mSetSq <- suppressMessages(preprocessFunnorm(rgSet))
    
    mSetSq <- rgSet
    #class(mSetSq)
    mSetSq <- preprocessRaw(mSetSq) #does not perform normalization 
    #class(mSetSq)
    mSetSq <- mapToGenome(mSetSq)
    #class(mSetSq)
    mSetSq <- ratioConvert(mSetSq)
    #class(mSetSq)
    
    #Ensure probes are in the same order in the mSetSq and detP objects
    detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
    
    #Remove any probes that have failed in one or more samples
    probes_before <- dim(mSetSq)[1]
    keep <- rowSums(detP < 0.01) == ncol(mSetSq)
    mSetSqFlt <- mSetSq[keep,]
    
    #Print out number of probes removed for failing in one or more samples
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for failing in one or more samples") #----- 5482 probe(s) removed for failing in one or more samples
    probes_before <- dim(mSetSqFlt)[1]
    
    #Remove probes with SNPs at CpG site
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
    
    #Print out number of probes removed for having SNPs at CpG site
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for having SNPs at CpG site") #----- 29130 probe(s) removed for having SNPs at CpG site
    probes_before <- dim(mSetSqFlt)[1]
    
    #Exclude cross reactive probes
    xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                           "ChenEtAlList.csv",
                                           sep="/"), stringsAsFactors=FALSE)
    keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep,]
    
    #Print out number of probes removed for being cross reactive
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for being cross reactive") #----- 25462 probe(s) removed for being cross reactive
    probes_before <- dim(mSetSqFlt)[1]
    
    #Remove Sex Probes
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    sexProbes <- ann[which(ann$chr %in% c("chrX", "chrY")),]
    keep <- !(featureNames(mSetSqFlt) %in% sexProbes$Name)
    mSetSqFlt <- mSetSqFlt[keep,]
    
    #Print out number of probes removed for being cross reactive
    probes_removed <- probes_before - dim(mSetSqFlt)[1]
    message("----- ", probes_removed, " probe(s) removed for being on sex chromosomes") #----- 17936 probe(s) removed for being on sex chromosomes
    
    #Print out the number of probes remaining
    message("----- ", dim(mSetSqFlt)[1], " probe(s) remaining for analysis") #----- 787849 probe(s) remaining for analysis
    
    #Calculate methylation beta values
    bVals <<- getBeta(mSetSqFlt)
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
      lme1 <- glmmTMB(formula, data = pdata, family = "gaussian")
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
    cairo_pdf(paste("matrixplot", age_type, ".pdf", sep = ""), width = 14, height = 14, fallback_resolution = 1000)
    print(pairs(plot.formula, data = pdataSVs, upper.panel = twolines, labels = diag.labels, diag.panel = mydiag.panel, lower.panel = panel.cor, label.pos = 0.5, main = ""))
    dev.off()
    
    finalOutput <- paste(output, "\n", age_type, "Covariates\n")
    finalOutput <- corCovariates(finalOutput)
    covariate_data <- pdataSVs
    
    covariate_data <- removeCovariates(covariate_data)
    
    finalOutput <- paste(finalOutput, "\n", age_type, "Residuals Based on Epigenetic Age", "\n")
    finalOutput <- residGeneration(pdata = covariate_data, output = finalOutput)
    finalOutput <- paste(finalOutput, "\n", age_type, "Residuals Based on Epigenetic Age Acceleration", "\n")
    covariate_data$Clock <- covariate_data$Clock - covariate_data$Age
    finalOutput <- residGeneration(pdata = covariate_data, output = finalOutput)
    
    assign("listofCors", c(), envir = .GlobalEnv)
    assign("corsToRemove", c(), envir = .GlobalEnv)
    
    return(finalOutput)
  }
  
  createGroupedBarChart <- function(data, x, y, fill, title) {
    melted_df <- reshape2::melt(data, id.vars = x, variable.name = fill)
    custom_palette <- c("age" = "red", "horvath" = "#66c2a5", "skinhorvath" = "#3288bd", "hannum" = "#5e4fa2", "levine"="#3288dd")
    
    plot <- ggplot(data = melted_df, aes_string(x = x, y = y, fill = fill)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = x, y = "Age", title = title) +
      scale_fill_manual(values = custom_palette) +  # Apply the custom palette
      theme_minimal()
    
    plot <- plot + theme(plot.background = element_rect(fill = "white"))
    
    ggsave("SampleIDandAge.png", plot = plot, width = 10000, height = 1000, units = "px", dpi = 300)
    
    for (i in 2:(ncol(data) - 1))
    {
      plot <- ggscatter(data, x = "age", y = colnames(data)[i],
                        add = "reg.line",  # 
                        add.params = list(color = "blue", fill = "lightgray")) 
      plot <- plot + stat_cor(method = "pearson")
      ggsave(filename = paste0("plot_", colnames(data)[i], ".png"),
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
  
  print(bVals)
  #PCA Generation####################################
  #Remove NAs
  #which(is.na(bVals)) # no NA values in data frame
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
  
  
  #Cell Type Composition
  
  if (exists("rgSet"))
  {
    library(FlowSorted.CordBlood.450k)
    CC <- estimateCellCounts(rgSet, compositeCellType = "CordBlood",
                             processMethod = "auto", probeSelect = "auto",
                             cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "nRBC"),
                             referencePlatform = "IlluminaHumanMethylation450k",
                             returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)
    
    CC <- removeOutliers(CC)
  }
  
  #Correlation Matrix##########################
  
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
  if (exists("rgSet"))
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
  
  pdataSVs <- pdataSVs %>%
    select_if(~n_distinct(.) >= 2)
  numRemoved <- ncol(pdataSVs) - ncol(sampleData)
  warning(paste(numRemoved," variables were removed because you either did not follow proper naming conventions or your variable had only 1 level "))
  
  assign("listofCors", c(), envir = .GlobalEnv)
  assign("corsToRemove", c(), envir = .GlobalEnv)
  
  diag.labels=colnames(pdataSVs)
  pdataColumns <- names(pdataSVs)
  plot.formula=as.formula(paste("~", paste(pdataColumns, collapse = "+")))
  
  
  #Clock Output#######################################
  
  results <- NULL
  
  if (shouldNormalize == TRUE) results <- DNAmAge(bVals, normalize = TRUE, age = sampleData$Age) else results <- DNAmAge(bVals, normalize = FALSE, age = sampleData$Age)
  
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
  
  results <- PACEProjector(as.matrix(bVals), proportionOfProbesRequired = 0.8)
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
  grimage <- methyAge(bVals, clock = clockname, age_info = grimDf)
  grimage$id <- grimage$Sample
  exportDf <- merge(exportDf, grimage, by="id")
  for (i in 1:nrow(results))
  {
    finalOutput <- paste(finalOutput, grimage[i,1], "\t", grimage$Age_Acceleration[i], "\n")
  }
  
  write.table(as.data.frame(exportDf), file = "epigeneticAge.txt")
  
  outputString <- paste(finalOutput, collapse = "\n")
  
  write_file(outputString, file = "output.txt")
  
  print(Sys.time() - start)
  
}

# 1 = 27K, 2 = 450K, 3 = EPIC

#main(directory = "C:/Users/stanl/Desktop/Development/R/DNAm-age-pipeline/data", useBeta = FALSE)




