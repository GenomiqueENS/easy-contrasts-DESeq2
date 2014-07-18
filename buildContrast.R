#!/usr/bin/env Rscript
###############################################################################
## R script for generation of contrast matrix using by  
## DESeq2 for complex comparison 
# 
## Author : Xavier Bauquet
###############################################################################

###############################################################################
# 
## Main
# 
###############################################################################
args <- commandArgs(TRUE) 
# Parameters
designPath <- args[1]
deseqModel <- args[2]
comparisonPath <- args[3]
contrastFile <- args[4]
verbose <- toupper(args[5])

if(verbose=="TRUE"){            
    # Print date, hours, packages version and params for log file 
    cat("\n\n########################################################\n")
    cat("Start of the DESeq2 contrast matrix generator script version 1.0\n")
    cat("########################################################\n\n")
    cat(format(Sys.time(), "%Y/%m/%d %H:%M:%S\n"))
    cat(R.version.string)
    cat(paste("\nDESeq2 version", packageVersion("DESeq2")))
    
    # Print params
    cat("\n\n########################\n")
    cat("Params\n")
    cat("########################\n")
    cat(paste("Path of the design file  =", args[1]))
    cat(paste("\nDESeq2 modele            =", args[2]))
    cat(paste("\nComparison file          =", args[3]))
    cat("\n\n########################\n\n")
}            
# loading of DESeq2 library
library(DESeq2)                
# -----------------------------------------------------------------------------                
    if(verbose=="TRUE")cat("\n\n########################\n")
    if(verbose=="TRUE")cat("1 - Read design file\n") 
# loading of the design file
target <- read.table(designPath, sep="\t", header=T, dec=".", stringsAsFactors=F)

# -----------------------------------------------------------------------------
    if(verbose=="TRUE")cat("2 - Artificial count matrix building\n")
# creation of the random count matrix for running the DESeq2 minimal script
forCounts <- rep(1, (nrow(target)*3))
for (i in 1:(nrow(target)*3)){
    forCounts[i] <- round(runif(1, min=1, max=1000000000), digits=0)
}
counts <- matrix(data=forCounts, nrow=3)

# -----------------------------------------------------------------------------
    if(verbose=="TRUE")cat("3 - DESeq2 function runing for B factors names\n") 
# DESeq2 minimal script                
dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, design=as.formula(deseqModel))
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, modelMatrixType="expanded")
# B factors 
Bfactors <- resultsNames(dds)

if(verbose=="TRUE"){
    cat("\n###############################################################################\n")
    cat("Beta factors:\n")
    print(Bfactors)
    cat("\n###############################################################################\n\n")
}
# Zero matrix model creation 
B <- rep(0, (length(Bfactors)-1))
# -----------------------------------------------------------------------------
    if(verbose=="TRUE")cat("3 - DESeq2 function runing for B factors names\n") 
# creation of the model matrix for each condition
factorData <- data.frame(condition= Bfactors[2:length(Bfactors)])
for(i in 1:length(Bfactors)){
    factorData <- data.frame(factorData, B)
    if(i > 1){
        factorData[i-1,i+1] <- 1
        factorData[i-1,2] <- 1
    }
}
# Creation of the final matrix to be export for DESeq2 analysis
contrastData <- factorData[1,]
names(contrastData)[1] <- "comparison"

# -----------------------------------------------------------------------------                
    if(verbose=="TRUE")cat("4 - Read comparison file\n") 
# loading of the comparison file
form <- read.table(comparisonPath, sep="\t", header=F, dec=".", stringsAsFactors=F)

# -----------------------------------------------------------------------------
    if(verbose=="TRUE")cat("5 - Add of model matrix for each condition of the interaction\n")
# For each comparison in the comparison file
for(x in 1:nrow(form)){
    # vs séparations
    results <- unlist(strsplit(as.vector(form[x,]), "_vs_")) 
    # Creation of a temporary matrix for saving
    tmpData <- factorData[1,2:ncol(factorData)]      
    for(y in 1:length(results)){
        # interaction separations
        toAdd <- unlist(strsplit(results[[y]], "%")) 
        # vector for the add of the contrast vectors in interaction
        add <- B
        # Test the presence of interaction factor (ex: statut:disease and add of contrast vector) 
        # for interaction factor to add vector
        for(i in 1:length(toAdd)-1){
            for(j in 2:length(toAdd)){
                # Test the presence of condition1:condition2
                test <- paste(toAdd[i], toAdd[j], sep=".")
                if(test %in% factorData$condition){
                    add <- add + subset(factorData, condition == test)[2:length(names(factorData))]
                }
                # Test the presence of condition2:condition1
                test <- paste(toAdd[j], toAdd[i], sep=".")
                if(test %in% factorData$condition){
                    add <- add + subset(factorData, condition == test)[2:length(names(factorData))]
                }
            }
        }
        # Add vectors of normal factors to the add vector
        for(matrix in toAdd){
            add <- add + subset(factorData, condition == matrix)[2:length(names(factorData))]
        }
        # saving of the add vector into a data frame
        tmpData <- rbind(tmpData, add)
    }
    
# -----------------------------------------------------------------------------                
    if(verbose=="TRUE")cat("6 - Substraction of matrix for _vs_ comparison\n")
    # suppression of the first row of the tmpData data frame kip during the creation of the data frame
    tmpData <- tmpData[2:nrow(tmpData),]
    # substraction for a comparison of 2 group (ex: group1_vs_group2)
    if(nrow(tmpData) == 2){
        sum <- tmpData[1,] - tmpData[2,]
    # substraction for a comparison of 4 group (ex: group1_vs_group2_vs_group3_vs_group4)
    }else if(nrow(tmpData) == 4){
        sum <- (tmpData[1,] - tmpData[2,])-(tmpData[3,]-tmpData[4,])
    }
    # saving into a data.frame
    tmpFrame <- data.frame(form[x,], sum)
    names(tmpFrame)[1] <- "comparison"
    contrastData <- rbind(contrastData, tmpFrame)
}

# -----------------------------------------------------------------------------
    if(verbose=="TRUE")cat("7 - Final comparison matrix saving\n")   
# creation of the final data frame for saving the contrast matrix
finalData <- data.frame(comparisons="", matrix="")
for(i in 1:nrow(contrastData)){
    # formatting of the contrast vector like the model (0,1,-1,0)
    newTmpFrame <- data.frame(comparisons=contrastData[i,1], matrix=paste("(", paste(contrastData[i,2:length(names(contrastData))], collapse=","), ")", sep=""))
    finalData <- rbind(finalData, newTmpFrame)
}
# suppression of the 2 first row of the data frame due to the construction
finalData <- finalData[3:nrow(finalData),]
# saving of the contrast matrix
write.table(finalData, contrastFile,sep="\t",row.names=F, quote=F)

                
    # Print date, hours
if(verbose=="TRUE"){
    cat("\n############\n")
    cat(format(Sys.time(), "%Y/%m/%d %H:%M:%S\n"))
    cat("Successful end of the analysis\n")
}
