#!/usr/bin/env Rscript
###############################################################################
## Normalisation and Differential analysis using DESeq2 R package 
## for RNA sequencing analysis
# 
## Author : Xavier Bauquet
###############################################################################
# -----------------------------------------------------------------------------
# buildCountMatrix
# Create a matrix of reads count
#
# Input:
#   files : a vector of files names
#   sampleLabel : a vector of sample names
#   projectPath: path to the project directory
#
# Ouput:
#   countMatrix : a reads count matrix
#
# Original author : Vivien DESHAIES
# -----------------------------------------------------------------------------
buildCountMatrix <- function(files, sampleLabel,expHeader){
    
    # read first file and create countMatrix
    countMatrix <- read.table(files[1],header=as.logical(expHeader),stringsAsFactors=F,quote="")[,c("Id","Count")]
    # lowercase countMatrix columns names
    colnames(countMatrix) <- tolower(colnames(countMatrix))
    
    # read and merge all remaining files with the first
    for(i in 2:length(files)){
        # read files
        exp <- read.table(files[i],header=as.logical(expHeader),stringsAsFactors=F,quote="")[,c("Id","Count")]
        # lowercase exp columns names
        colnames(exp) <- c("id", paste("count", i, sep=""))
        # merge file data to count matrix by id
        countMatrix <- merge(countMatrix, exp, by="id", suffixes="_") 
    }
    # name rows
    rownames(countMatrix) <- countMatrix[,1]
    # delete first row containing row names
    countMatrix <- countMatrix[,-1]
    
    # name columns
    colnames(countMatrix) <- sampleLabel
    return(countMatrix)
}
###############################################################################
# ----------------------------------------------------------------------------
# deleteUnexpressedGene
# 
# Delete genes row whithout read count in a count matrix
#
# Original author : Vivien DESHAIES
# -----------------------------------------------------------------------------
deleteUnexpressedGene <- function(countMatrix){
    # Delete empty rows
    countMatrix <- countMatrix[rowSums(countMatrix)>0,]
    return(countMatrix)
}
###############################################################################
# -----------------------------------------------------------------------------
# saveCountMatrix
# 
# Input: 
#   countMatrix : a read count matrix
#   outpath : the path where to save file
#   fileName : the file name
#
# Original author : Vivien DESHAIES
# -----------------------------------------------------------------------------

saveCountMatrix <- function(countMatrix,fileName){
    # Add ID column to make the file openable by calc or excel
    i <- length(countMatrix[1,])
    countMatrix <- cbind(countMatrix, rep(0, length(countMatrix[,1])))
    while(i >= 1){
        countMatrix[,i+1] <- countMatrix[,i]
        colnames(countMatrix)[i+1] <- colnames(countMatrix)[i]
        i <- i-1
    }
    countMatrix[,1] <- rownames(countMatrix)
    colnames(countMatrix)[1] <- "Id"
    # Write the count matrix in a file
    write.table(countMatrix, paste(fileName, sep=""),sep="\t",row.names=F, quote=F)
}
###############################################################################
# -----------------------------------------------------------------------------
# printInformationStart
#
#   This function print several informations about the script at the begining: 
#   the time, the R version used, the versions of packages used in the script 
#   and the parameters
#
# -----------------------------------------------------------------------------
printInformationStart <- function(args){
    cat("\n\n########################################################\n")
    cat("Start of the DESeq2 script version 1.0\n")
    cat("########################################################\n\n")
    cat(format(Sys.time(), "%Y/%m/%d %H:%M:%S\n"))
    cat(R.version.string)
    cat(paste("\nDESeq2 version", packageVersion("DESeq2")))
    cat(paste("\nRColorBrewer version", packageVersion("RColorBrewer")))
    cat(paste("\nFactoMineR version", packageVersion("FactoMineR")))
    cat("\n\n########################\n")
    cat("Params\n")
    cat("########################\n")
    cat(paste("Figures                                        =", toupper(args[1])))
    cat(paste("\nDifferential analysis                          =", toupper(args[2])))
    cat(paste("\nFigures of the differential analysis           =", toupper(args[3])))
    cat(paste("\nContrast matrix for differential analysis      =", args[4]))
    cat(paste("\nName of the design file                        =", args[5]))
    cat(paste("\nDESeq2 modele                                  =", args[6]))
    cat(paste("\nProject name                                   =",args[7]))
    cat(paste("\nHeader on expression files                     =",toupper(args[8])))
    cat("\n\n########################\n\n")
}
###############################################################################
# -----------------------------------------------------------------------------
# buildColorVector
#   
#   This function aims to prepare a vector of colors, one for each unique 
#   condition of the design
#
#   input: design -> data frame (the design file)
#   output: coLors -> vector (the vector of colors)
#
# -----------------------------------------------------------------------------
buildColorVector <- function(design){
    # for a 2 conditions analysis
    if(length(unique(design$Condition))<= 2){
        coLors <- c("#A6CEE3","#1F78B4")
    # for a 3-12 conditions analysis, using of a paired set of colors
    }else if(2 < length(unique(design$Condition)) &&  length(unique(design$Condition))<= 12){
        # download of the "Paired" set of colors from the RColorBrewer library
        uniqueColors <- brewer.pal(length(unique(design$Condition)), "Paired")
        # selection of the good number of colors for the analysis 
        test <- lapply(design$Condition ,function(x){x == unique(design$Condition)})
        coLors <- c()
        for (result in test){
            coLors <- c(coLors, uniqueColors[result])
        }
    # for an analysis with more than 12 conditions
    }else{
        uniqueColors <- rainbow(length(unique(design$Condition)))
        # selection of the good number of colors for the analysis
        test <- lapply(design$Condition ,function(x){x == unique(design$Condition)})
        coLors <- c()        
        for (result in test){
            coLors <- c(coLors, uniqueColors[result])
        }
    }
    return(coLors)
}
###############################################################################
# -----------------------------------------------------------------------------
# firstPlots
#
#   This function create 3 plots for raw count matrix: unpooled clustering plot, 
#   unpooled PCA plot and unpooled null counts barplot
#
#   input: verbose -> booleen (for verbose mode)
#          projectName -> character (name of the project)
#          count_mat -> data frame (the count matrix)
#   output: 3 plots -> png
#
# -----------------------------------------------------------------------------
firstPlots <- function(verbose, projectName, count_mat){
    if(verbose=="TRUE")cat("      Fig 1 - Unpooled clustering\n") 
    png(paste("unpooled_clustering-", projectName,".png", sep=""),width=1000, height=600)
        # calculation of the dispersion
        dist.mat <- dist(1 - cor(count_mat)/2)
        plot(hclust(dist.mat), main=paste("non_pooled cluster dendrogram - ", projectName, sep=""), xlab="")
    dev.off()
    
    if(verbose=="TRUE")cat("      Fig 2 - Unpooled PCA\n") 
    pcaCount <- PCA(t(count_mat), graph=FALSE)
    png(paste("PCA - ",projectName,".png",sep=""),width=1000, height=600)
        par(mar=c(5,5,5,20))
        plot.PCA(pcaCount, choix="ind", col.ind=as.character(design$coLors), title = paste("PCA ", projectName, sep=""))
        cor<-par('usr')
        par(xpd=NA)
        # add of legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(design$Condition), col=unique(as.character(design$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
 
    if(verbose=="TRUE")cat("      Fig 3 - Null counts barplot\n") 
    png(paste("Null_counts-", projectName,".png"),width=1000, height=600)
        par(mar=c(15,8,5,20))
        barplot(100*colMeans(count_mat==0), cex.lab=2, las=3, col=as.character(design$coLors) ,main=paste("Proportion of null counts per sample -", projectName, sep=" "), ylab="Proportion of null counts (%)")
        cor<-par('usr')
        par(xpd=NA)
        # add of legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(design$Condition), col=unique(as.character(design$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
}
###############################################################################
# -----------------------------------------------------------------------------
# secondPlots
#
#   This function create 2 plots for raw count marix after deletion of unexpressed
#   genes and convertion of the count matrix and the design file in a DESeq object:
#   unpooled counts barplot, unpooled counts boxplot
#
#   input: verbose -> booleen (for verbose mode)
#          projectName -> character (name of the project)
#          dds -> DESeq object
#   output: 2 plots -> png
#
# -----------------------------------------------------------------------------
secondPlots <- function(verbose, projectName, dds){
    if(verbose=="TRUE")cat("      Fig 4 - Unpooled counts barplot\n") 
    png(paste("barplot_counts-", projectName,".png", sep=""),width=1000, height=600)
        par(mar=c(15,8,5,20))
        barplot(colSums(counts(dds)), main=paste("Read counts - ", projectName, sep=""), col=as.character(colData(dds)$coLors), names.arg =colData(dds)$Name,cex.lab=2, las=3, ylab="Total expression counts")
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()

    if(verbose=="TRUE")cat("      Fig 5 - Unpooled counts boxplot\n") 
    png(paste("boxplot_count-",projectName,".png", sep=""),width=1000, height=600)
        par(mar=c(15,8,5,20))
        boxplot(log2(counts(dds)+1),main=paste("Count distribution - ", projectName, sep=""), col=as.character(colData(dds)$coLors), names =colData(dds)$Name,cex.lab=2, las=3, ylab="log2 (counts+1)")
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
}
###############################################################################
# -----------------------------------------------------------------------------
# pooledPlots
#
#   This function create 2 plots for count matrix after collapsing of technical 
#   replicates: pooled counts barplot and pooled counts boxplot
#
#   input: verbose -> booleen (for verbose mode)
#          projectName -> character (name of the project)
#          dds -> DESeq object
#   output: 2 plots -> png
# -----------------------------------------------------------------------------
pooledPlots <- function(verbose, projectName, dds){

    if(verbose=="TRUE")cat("      Fig 6 - Pooled counts barplot\n") 
    png(paste("barplot_counts_pooled-", projectName,".png", sep=""),width=1000, height=600)
        par(mar=c(15,8,5,20))
        barplot(colSums(counts(dds)), main=paste("Pooled expression counts - ", projectName, sep=""), col=as.character(colData(dds)$coLors), names.arg =colData(dds)$Name,cex.lab=2, las=3, ylab="Total expression counts")
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()

    if(verbose=="TRUE")cat("      Fig 7 - Pooled counts boxplot\n") 
    png(paste("boxplot_count_pooled", projectName,".png", sep=""),width=1000, height=600)
        par(mar=c(15,8,5,20))
        boxplot(log2(counts(dds)+1),main=paste("Pooled count distribution - ", projectName, sep=""), col=as.character(colData(dds)$coLors), names =colData(dds)$Name,cex.lab=2, las=3, ylab="log2 (counts+1)")
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
}
###############################################################################
# -----------------------------------------------------------------------------
# normPlots
#
#   This function create 4 plots for count matrix after collapsing of technical
#   replicates and normalisation: pooled and normalised clustering, pooled and
#   normalised PCA, pooled and normalised boxplot and most expressed sequence plot
#
#   input: verbose -> booleen (for verbose mode)
#          projectName -> character (name of the project)
#          dds -> DESeq object
#   output: 2 plots -> png  
#
# -----------------------------------------------------------------------------
normPlots <- function(verbose, projectName, dds){

    if(verbose=="TRUE")cat("      Fig 8 - Pooled and Normalised clustering\n") 
    png(paste("pooled_clustering-",projectName,".png", sep=""),width=1000, height=600)
        # calculation of the dispertion 
        ddsStabilized <- assay(varianceStabilizingTransformation(dds))
        dist.mat <- dist(1 - cor(ddsStabilized)/2)
        plot(hclust(dist.mat), main=paste("Pooled and Normalised cluster dendrogram - ", projectName, sep=""), xlab="")
    dev.off()

    if(verbose=="TRUE")cat("      Fig 9 - Pooled and Normalised PCA\n") 
    pcaCount <- PCA(t(counts(dds, normalized=TRUE)), graph=FALSE)
    png(paste("normalised_PCA - ",projectName,".png",sep=""),width=1000, height=600)
        par(mar=c(5,5,5,20))
        plot.PCA(pcaCount, choix="ind", col.ind=as.character(colData(dds)$coLors), title = paste("Pooled and normalised PCA ", projectName, sep=""))
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
        
    if(verbose=="TRUE")cat("      Fig 10 - Pooled and Normalised boxplot\n") 
    png(paste("normalised_boxplot_count-",projectName,".png",sep=""),width=1000, height=600)
        par(mar=c(15,8,5,20))
        boxplot(log2(counts(dds,normalized=TRUE)+1),main=paste("Pooled and Normalised count distribution - ", projectName, sep=""), col=as.character(colData(dds)$coLors), names =colData(dds)$Name,cex.lab=2, las=3, ylab="log2 (counts+1)")
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
        
    if(verbose=="TRUE")cat("      Fig 11 - Most expressed sequence plot\n") 
    # preparation of 2 data frame with the same number of column than the dds count matrix 
    maxCounts <- counts(dds)[1,]
    transcriptNames <- counts(dds)[1,]
    # for each sample (column)
    for(i in 1:ncol(counts(dds))){
        # selection of the maximum number of count
        maxCounts[i] <- (max(counts(dds, normalized=TRUE)[,i])/sum(counts(dds, normalized=TRUE)[,i]))*100
        # selection of the name of the gene/transcript this the maximum of count
        transcriptNames[i] <- row.names(subset(counts(dds, normalized=TRUE), counts(dds, normalized=TRUE)[,i]==max(counts(dds, normalized=TRUE)[,i])))
    }
    png(paste("most_expressed_sequence-",projectName,".png",sep=""),width=1000, height=600)
        par(mar=c(15,8,5,20))
        x <- barplot(maxCounts, main=paste("Most expressed sequences - ", projectName, sep=""), col=as.character(colData(dds)$coLors), names.arg =colData(dds)$Name,cex.lab=2, las=3, ylab="Proportion of reads (%)")
        # add names of the gene/transcript on the plot bars
        text(x, 0, labels= transcriptNames, srt=90, adj=0)
        cor<-par('usr')
        par(xpd=NA)
        # add legends
        legend(cor[2]*1.01,cor[4], title="Legend", legend=unique(colData(dds)$Condition), col=unique(as.character(colData(dds)$coLors)), pch=15, pt.cex=3, cex=1.2)
    dev.off()
}
###############################################################################
# -----------------------------------------------------------------------------
# anadiff
#
#   This function perform the classical differential analysis without contrast 
#   vector: comparison of 2 conditions from the Condition colomne of the design
#   file
#
#   input: dds -> DESeq object
#          condition1, condition2 -> vectors (names of the conditions to compare)
#          param -> booleen (FALSE to escape plots)
#          projectName -> character (name of the project)
#   output: 4 plots -> png (using the anadiffPlots function)
#           diffana matrix of the comparison -> tsv (file containing results of 
#               the differential analysis between both conditions)
#   
# -----------------------------------------------------------------------------
anadiff <- function(dds, condition1, condition2, param, projectName){
    # selection of results of the comparison
    res <- results(dds, contrast=c("Condition", condition1, condition2))
    # function for plots
    if(param=="TRUE")anadiffPlots(paste(condition1,"_vs_", condition2,sep=""), projectName,res)
    saveCountMatrix(as.data.frame(res),paste("diffana_",condition1,"_vs_",condition2,"-",projectName,".tsv", sep=""))
}
###############################################################################
# -----------------------------------------------------------------------------
# contrastAnadiff
#
#   This function perform the differential analysis using contrast vector
#
#   input: dds -> DESeq object
#          nameContrastVec -> character (name of the contrast vector)
#          contrastVec -> vector (the contrast vector)
#          param -> booleen (FALSE to escape plots)
#          projectName -> character (name of the project)
#   output: 4 plots -> png (using the anadiffPlots function)
#           diffana matrix of the comparison -> tsv (file containing results of 
#               the differential analysis between both conditions)
# 
# -----------------------------------------------------------------------------
contrastAnadiff <- function(dds, nameContrastVec, contrastVec, param, projectName){
    # separation of the character string in a vector
    contrastVec <- unlist(strsplit(substr(contrastVec, 2, (nchar(contrastVec)-1)), ","))
    # selection of results with the contrast vector
    res <- results(dds, contrast= as.numeric(contrastVec))
    # change the % in the name in -
    nameContrastVec <- gsub("%","-",nameContrastVec)
    # function for plots
    if(param=="TRUE")anadiffPlots(nameContrastVec, projectName,res)   
    saveCountMatrix(as.data.frame(res),paste("diffana_",nameContrastVec,"-",projectName,".tsv", sep=""))
    cat(paste("Comparison: ",nameContrastVec," finish\n", sep=""))
}
###############################################################################
# -----------------------------------------------------------------------------
# anadiffPlots
#
#   This function create 4 plots for the both differential analysis function
#   (anadiff and contrastAnadiff): p-value plot, adjusted p-value plot, MA-plot
#   and differentially expressed gense according adjusted P-value plot
#
#   input: condName -> character (name of the comparison or of the contrast vector)
#          projectName -> character (name of the project)
#          res -> DESeq object (results of the differential analysis)
#   output: 4 plots -> png
#
# -----------------------------------------------------------------------------
anadiffPlots <- function(condName, projectName,res){
    # Raw Pvalue plot
    png(paste("diffana_plot_pvalue_",condName,"-",projectName,".png",sep=""),width=1000, height=600)
        hist(res$pvalue,main=paste("Raw P-value plot - ", projectName,sep=""), col="#1F78B4")
    dev.off()
    # Adjusted Pvalue plot
    png(paste("diffana_plot_padj_",condName,"-",projectName,".png",sep=""),width=1000, height=600)
        hist(res$padj,main=paste("Adjusted P-value plot - ", projectName,sep=""), col="#FF7F00")
    dev.off()
    # MA-plot      
    png(paste("diffana_MA_plot_",condName,"-",projectName,".png",sep=""),width=1000, height=600) 
        plotMA(res, alpha=0.05, ylim=c(-10,10), main= paste("MA-plot of ",condName," - ", projectName,sep=""))
    dev.off()
    # Number of differentially expressed genes according to padj
    # adjusted p-value to plot
    p <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05)
    value <- c(0,0,0,0,0,0)
    # calculation of the number of genes/transcript differentially expressed for each adjusted p-value
    for(i in 1:length(p)){
        value[i] <- nrow(subset(res, res$padj < p[i]))
    }
    png(paste("diffana_plot_differentially_expressed_genes_",condName,"-",projectName,".png",sep=""),width=600, height=600) 
        par(mar=c(8,8,5,5))
        x <- barplot(value, main=paste("Differentially expressed genes according adjusted P-value - ", projectName,sep=""), col="tan", names.arg =as.character(p),cex.lab=2, ylab="Number of genes differentially expressed")
        # add of the number of genes differentially expressed on the plot bars
        text(x, 0, labels= value, srt=90, adj=0)
        cor<-par('usr')
        par(xpd=NA)
    dev.off()    
}
###############################################################################
# -----------------------------------------------------------------------------
# printInformationEnd
#
#   This function print the time at the end of the script 
#
# -----------------------------------------------------------------------------
printInformationEnd <- function(){
    cat("\n############\n")
    cat(format(Sys.time(), "%Y/%m/%d %H:%M:%S"))
    cat("\nSuccessful end of the analysis\n")
    cat("\n############\n")
}

###############################################################################
# 
## Main
# 
###############################################################################
args <- commandArgs(TRUE)
# -----------------------------------------------------------------------------
# Parameters
normFigTest <- toupper(args[1])
diffanaTest <- toupper(args[2])
diffanaFigTest <- toupper(args[3])
contrastTest <- args[4]
designPath <- args[5]
deseqModel <- args[6]
projectName <- args[7]
expHeader <- toupper(args[8])
verbose <- toupper(args[9])
# -----------------------------------------------------------------------------
        if(verbose =="TRUE"){
            printInformationStart(args)
        }
# load libraries
library(DESeq2)
library(RColorBrewer)
library(FactoMineR)
# -----------------------------------------------------------------------------
        if(verbose=="TRUE"){
            cat("\n\n########################\n")
            cat("1 - Read design file\n")
        }
# load design file
design <- read.table(designPath, sep="\t", header=T, dec=".", stringsAsFactors=F)

        if(verbose=="TRUE")cat("2 - Creation of the color vector\n")
coLors <- buildColorVector(design)
# add colors to the design
design <- data.frame(design, coLors)

        if(verbose=="TRUE")cat("3 - Count matrix building\n") 
fileNames <- paste("expression_",design$SampleNumber,".tsv",sep="") 
# computing of expression files in one unique file
# count_mat <- buildCountMatrix(expressionFiles, design$Name)
count_mat <- buildCountMatrix(fileNames, design$Name, expHeader)

### plots: unpooled clustering plot, unpooled PCA plot and unpooled null counts barplot
if(normFigTest=="TRUE")firstPlots(verbose, projectName, count_mat)
###

        if(verbose=="TRUE")cat("4 - Unexpressed gene deletion\n") 
count_mat <- deleteUnexpressedGene(count_mat)

        if(verbose=="TRUE")cat("5 - DESeq2 object building\n")
# creation of the DESeq object including the count matrix and the design file
dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=design, design=as.formula(deseqModel))

### plots: unpooled counts barplot, unpooled counts boxplot
if(normFigTest=="TRUE")secondPlots(verbose, projectName, dds)
###

        if(verbose=="TRUE")cat("6 - Saving of rawCountMatrix\n")
saveCountMatrix(counts(dds),paste("diffana_", projectName,"_rawCountMatrix.tsv", sep=""))

# -----------------------------------------------------------------------------
#
#   Collapsing technical replicates
#
# -----------------------------------------------------------------------------
        if(verbose=="TRUE")cat("7 - Collapsing of technical replicas\n") 
dds <- collapseReplicates(dds, groupby=dds$RepTechGroup)

### plots: pooled counts barplot and pooled counts boxplot
if(normFigTest=="TRUE")pooledPlots(verbose, projectName, dds)
###

        if(verbose=="TRUE")cat("8 - Saving of rawPooledCountMatrix\n")
saveCountMatrix(counts(dds),paste("diffana_", projectName, "_rawPooledCountMatrix.tsv", sep=""))

# -----------------------------------------------------------------------------
#
#   Normalisation
#
# -----------------------------------------------------------------------------
        if(verbose=="TRUE")cat("9 - Normalisation\n")
dds <- estimateSizeFactors(dds)

### plots: pooled and normalised clustering, pooled and normalised PCA, pooled and normalised boxplot and most expressed sequence plot
if(normFigTest=="TRUE")normPlots(verbose, projectName, dds)
###

        if(verbose=="TRUE")cat("10 - Saving of normalisedCountMatrix\n")
saveCountMatrix(counts(dds,normalized=TRUE),paste("diffana_", projectName,"_normalisedCountMatrix.tsv", sep=""))

# -----------------------------------------------------------------------------
#
#   Differential analysis
# 
# -----------------------------------------------------------------------------

if(diffanaTest=="TRUE"){
        if(verbose=="TRUE")cat("11 - Dispersion estimations\n")
    dds <- estimateDispersions(dds)
    if(diffanaFigTest=="TRUE"){
        if(verbose=="TRUE")cat("      Fig 12 - Dispersion plot\n") 
        png("plot_disp.png",width=1000, height=600)
            plotDispEsts(dds)
        dev.off()
    }
    # if differential analysis using contrast matrix
    if(contrastTest != "FALSE"){
        if(verbose=="TRUE")cat("12 - Differential analysis using contrast matrix\n")
        # statistical analysis
        dds <- nbinomWaldTest(dds, modelMatrixType="expanded")
        contrastMatrix <- read.table(contrastTest, sep="\t", header=T, dec=".", stringsAsFactors=F)
        for(i in 1:nrow(contrastMatrix)){
            contrastAnadiff(dds, contrastMatrix[i,1],contrastMatrix[i,2], diffanaFigTest, projectName)
        }
        
    }else{    
        # statistical analysis
        dds <- DESeq(dds)
        # selection of the reference condition
        TrueRef <- subset(design, Reference %in% "true")
        # selection of conditions
        unique_condition <- unique(design$Condition)
        # if no reference condition
        if(nrow(TrueRef)<1){
            if(verbose=="TRUE")cat("12 - Differential analysis without contrast matrix and without reference condition\n")        
            for(i in 1:(length(unique_condition)-1)){
                for(j in (i+1):length(unique_condition)){
                    if(unique_condition[i] != unique_condition[j]){
                        anadiff(dds, unique_condition[i], unique_condition[j], diffanaFigTest, projectName)
                    }
                }
            }
        # if reference condition
        }else{
            if(verbose=="TRUE")cat("12 - Differential analysis without contrast matrix and with reference condition\n") 
            for(i in 1:length(unique_condition)){
                if(unique_condition[i] != TrueRef$Condition[1]){
                    anadiff(dds, TrueRef$Condition[1], unique_condition[i], diffanaFigTest, projectName) 
                }
            }
        }
    }
}

# End printing
if(verbose =="TRUE"){
    printInformationEnd()
}