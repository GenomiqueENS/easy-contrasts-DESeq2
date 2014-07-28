#!/usr/bin/env Rscript
###############################################################################
## R script for generation of contrast matrix using by  
## DESeq2 for complex comparison 
# 
##  Version 1.3 (07/27/2014)
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

         
    # Print date, hours, packages version and parameters for log file 
    cat("\n\n########################################################\n")
    cat("Start of the DESeq2 contrast matrix generator script version 1.2\n")
    cat("########################################################\n\n")
    cat(format(Sys.time(), "%Y/%m/%d %H:%M:%S\n"))
    cat(R.version.string)
    cat(paste("\nDESeq2 version", packageVersion("DESeq2")))
    
    # Print parameters
    cat("\n\n########################\n")
    cat("Params\n")
    cat("########################\n")
    cat(paste("Path of the design file  =", args[1]))
    cat(paste("\nDESeq2 modele            =", args[2]))
    cat(paste("\nComparison file          =", args[3]))
    cat("\n\n########################\n\n")
         
# Loading of DESeq2 library
library(DESeq2)                
# -----------------------------------------------------------------------------                
    cat("\n\n########################\n")
    cat("1 - Read design file\n") 
# Loading of the design file
target <- read.table(designPath, sep="\t", header=T, dec=".", stringsAsFactors=F)

# -----------------------------------------------------------------------------
    cat("2 - Artificial count matrix building\n")
# Creation of the random count matrix for running the DESeq2 minimal script to pick up the B factors
forCounts <- rep(1, (nrow(target)*3))
for (i in 1:(nrow(target)*3)){
    forCounts[i] <- round(runif(1, min=1, max=1000000000), digits=0)
}
counts <- matrix(data=forCounts, nrow=3)

# -----------------------------------------------------------------------------
    cat("3 - DESeq2 function runing for B factors names\n") 
# DESeq2 minimal script for pick up the B factors                
dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, design=as.formula(deseqModel))
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, modelMatrixType="expanded")
# B factors 
Bfactors <- resultsNames(dds)


    cat("\n###############################################################################\n")
    cat("Beta factors:\n")
    print(Bfactors)
    cat("\n###############################################################################\n\n")

# Zero vector of the length of Bfactors-1 (- the Intercept factor) that will be used in column for the construction of the factorData
B <- rep(0, (length(Bfactors)-1))
# -----------------------------------------------------------------------------
    cat("3 - DESeq2 function runing for B factors names\n") 
# Creation of the model matrix for each condition
factorData <- data.frame(condition= Bfactors[2:length(Bfactors)])
for(i in 1:length(Bfactors)){
    factorData <- data.frame(factorData, B)
    if(i > 1){
        factorData[i-1,i+1] <- 1
    }
}
# Creation of the final data frame to be export for DESeq2 analysis. This data frame is construct from the factorData data frame
# to have same number and names of column. The new data frame contrastData include a first row that will be remove after filling
contrastData <- factorData[1,]
names(contrastData)[1] <- "comparison"

# -----------------------------------------------------------------------------                
    cat("4 - Read comparison file\n") 
# Loading of the comparison file
comparisonRow <- read.table(comparisonPath, sep="\t", header=F, dec=".", stringsAsFactors=F)

# -----------------------------------------------------------------------------
    cat("5 - Addition of model matrix for each condition of the interaction\n")
# For each row of the comparison file (for each comparison)
for(x in 1:nrow(comparisonRow)){
    # Separation of the comparison formula by the "_vs_" to obtain a vector including all %-interaction to compare
    comparisonFormulas <- unlist(strsplit(as.vector(comparisonRow[x,2]), "_vs_")) 
    # Creation of a temporary data frame for saving the contrast vectors from each %-interaction. This data frame is 
    # build from the factorData data frame kipping the first row without the "condition" column to kip only column
    # corresponding to the contrast vector. The first line kip during the construction of the data frame will be 
    # remove after filling
    tmpData <- factorData[1,2:ncol(factorData)]
    # For each %-interaction of the "_vs_" comparison
    for(y in 1:length(comparisonFormulas)){
        # %-interaction separations to obtain a vector 
        toAdd <- unlist(strsplit(comparisonFormulas[y], "%")) 
        # Creation of the add vector for the addition of the contrast vectors in %-interaction
        add <- 0
        # Test the presence of :-interaction factor (ex: status:disease) and add of contrast vector 
        # for :-interaction factor to the add vector
        for(i in 1:length(toAdd)-1){
            for(j in 2:length(toAdd)){
                # Test the presence of condition1:condition2
                test <- paste(toAdd[i], toAdd[j], sep=".")
                if(test %in% factorData$condition){
                    # Addition of the contrast vector of :-interaction factor to the add vector
                    add <- add + subset(factorData, condition == test)[2:length(names(factorData))]
                }
                # Test the presence of condition2:condition1
                test <- paste(toAdd[j], toAdd[i], sep=".")
                if(test %in% factorData$condition){
                    # Addition of the contrast vector of :-interaction factor to the add vector
                    add <- add + subset(factorData, condition == test)[2:length(names(factorData))]
                }
            }
        }
        # Pick up column names of interest in the deseq model formula. To select only the names of the interest columns of
        # the design file, the deseqModel is separated by the "+" and the first character "~" is remove. In the columnNames 
        # vector obtained names including the ":" character are remove to kip only column names of interest. 
        columnNames <- unlist(strsplit(as.vector(substr(deseqModel, 2, nchar(deseqModel))), "+",fixed = TRUE))
        columnNames <- columnNames[setdiff(seq(columnNames), grep(":", columnNames))]

        # For each names on columnNames
        for( w in 1:length(columnNames)){
            # Selection into the toDiv vector of all %-interaction condition that begin be the names of the column. This 
            # step aims to test the number of condition from the same column into the %-interaction. If there is more than 
            # one condition from the column each B factor have to be divide by the number of condition from the same column.
            toDiv <- grep(paste("^", columnNames[w], sep=""), toAdd, value = TRUE)
            # If at least one condition used from the column
            if(length(toDiv)>0){
                # For each condition of the column
                for( p in 1:length(toDiv)){
                    # Pick up of the contrast vector of the condition 
                    vec <- subset(factorData, condition == toDiv[p])[2:length(names(factorData))]
                    # The vector is divide by the number of conditions from the same column
                    vec <- vec/length(toDiv)
                    # Addition of the vector to the add vector
                    add <- add + vec
                }
            }
        }
        # Save the add vector into the temporary data frame
        tmpData <- rbind(tmpData, add)
    }
    
# -----------------------------------------------------------------------------                
    cat("6 - Substraction of matrix for _vs_ comparison\n")
    # Remove of the first row of the tmpData data frame kip during its construction
    tmpData <- tmpData[2:nrow(tmpData),]
    
    # Subtraction for a comparison of 2 group (ex: group1_vs_group2)
    if(nrow(tmpData) == 2){
        sum <- tmpData[1,] - tmpData[2,]
    # Subtraction for a comparison of 4 group (ex: group1_vs_group2_vs_group3_vs_group4)
    }else if(nrow(tmpData) == 4){
        sum <- (tmpData[1,] - tmpData[2,])-(tmpData[3,]-tmpData[4,])
    }
    # Saving into the contrastData  data.frame
    tmpFrame <- data.frame(comparisonRow[x,2], sum)
    names(tmpFrame)[1] <- "comparison"
    contrastData <- rbind(contrastData, tmpFrame)
}
# Remove of the first row of the contrastData that was kip during its construction
contrastData <- contrastData[2:nrow(contrastData),]

# -----------------------------------------------------------------------------
    cat("7 - Final comparison matrix saving\n")   
# Creation of the final data frame for saving the contrast matrix. This step create an empty row that will be remove after filling
finalData <- data.frame(name="", comparisons="", matrix="")

for(i in 1:nrow(contrastData)){
    # Formatting of the contrast vectors like the model (0,1,-1,0)
    newTmpFrame <- data.frame(name=as.character(comparisonRow[i,1]) ,comparisons=contrastData[i,1], matrix=paste("(", paste(contrastData[i,2:length(names(contrastData))], collapse=","), ")", sep=""), stringsAsFactors=F)
    # Saving into the final data frame
    finalData <- rbind(finalData, newTmpFrame)
}
# Remove of the first row of the finalData that was kip during its construction
finalData <- finalData[2:nrow(finalData),]
# Saving of the contrast matrix
write.table(finalData, contrastFile,sep="\t",row.names=F, quote=F)

                
    # Print date, hours

    cat("\n############\n")
    cat(format(Sys.time(), "%Y/%m/%d %H:%M:%S\n"))
    cat("Successful end of the contrast vectors building\n")

