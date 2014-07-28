#!/usr/bin/python
# -*- coding: utf-8 -*-

##############################################################
# runDESeq2
#
# Script to run and formatte options for both scripts 
# buildContrast.R and normDiffana.R
# 
# Version 1.3 (07/28/2014)
#
# Author: Xavier Bauquet
##############################################################
import os, sys
from optparse import OptionParser, OptionGroup
usage = "usage: %prog [options] arg1 arg2"

parser = OptionParser()

group = OptionGroup(parser, "Compulsory option")
group.add_option("-m", "--model", type = "string", dest = "model", help ="The DESeq2 model to be used for the differential analysis as a formula: ~condition1+condition2+condition1:condition2 for more information about the formula please refer to the DESeq2 documentation http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html")
parser.add_option_group(group)
group = OptionGroup(parser, "Option for analysis with contrast vectors:")
group.add_option("-c", "--contrast", type = "string", dest = "contrast", default = "FALSE", help ="logical TRUE/FALSE to escape or not the use of contrast vectors for differential anlysis. Default = FALSE")
group.add_option("-b", "--buildContrast", type = "string", dest = "buildContrast", default = "FALSE", help ="logical TRUE/FALSE to escape or not the building of contrast vectors. Default = FALSE")
parser.add_option_group(group)
group = OptionGroup(parser, "Files options")
group.add_option("-f", "--designFile", type = "string", dest = "designFile", default = "deseqDesign.txt", help ="Name of the design file\nDefault = deseqDesign.txt")
group.add_option("-C", "--comparisonFile", type = "string", dest = "comparisonFile", default = "comparisonFile.txt", help ="Name of the file including the comparison to be computed in contrast vector. Comparisons must have this form: condition1%condition2_vs_condition3%condition4 where the % mean interaction between conditions and _vs_ mean a comparison between interactions. Default = comparisonFile.txt")

parser.add_option_group(group)
group = OptionGroup(parser, "Other options")
group.add_option("-n", "--normFig", type = "string", dest = "normFig", default = "TRUE", help ="logical TRUE/FALSE to escape or not figures from the normalisation part. Default = TRUE")
group.add_option("-D", "--diffana", type = "string", dest = "diffana", default = "TRUE", help ="logical TRUE/FALSE to escape or not the differential analysis. Default = TRUE")
group.add_option("-d", "--diffanaFig", type = "string", dest = "diffanaFig", default = "TRUE", help ="logical TRUE/FALSE to escape or not figures from the differential analysis. Default = TRUE")
group.add_option("-p", "--projectName", type = "string", dest = "projectName", default = "exp1", help ="Name of the project. Default = exp1")
group.add_option("-H", "--expHeader", type = "string", dest = "expHeader", default = "TRUE", help ="logical TRUE/FALSE, TRUE if expression files have a header. Default = TRUE")
group.add_option("-N", "--normDiffana", type = "string", dest = "normDiffana", default = "TRUE", help ="logical TRUE/FALSE, FALSE to escape the normalisation and differential analysis steps. Default = TRUE")
parser.add_option_group(group)

options, arguments = parser.parse_args()

##############################################################
# Error messages
if not options.model:
    sys.exit("Model option is missing please use -m or --model") 
###    
if options.contrast.upper() != "TRUE" and options.contrast.upper() != "FALSE":
   sys.exit("error: -c --contrast option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")
if options.buildContrast.upper() != "TRUE" and options.buildContrast.upper() != "FALSE":
   sys.exit("error: -b --buildContrast option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section") 
###   
if options.normFig.upper() != "TRUE" and options.normFig.upper() != "FALSE":
    sys.exit("error: -n --normFig option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")
if options.diffana.upper() != "TRUE" and options.diffana.upper() != "FALSE":
    sys.exit("error: -D --diffana option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")
if options.diffanaFig.upper() != "TRUE" and options.diffanaFig.upper() != "FALSE":
    sys.exit("error: -d --diffanaFig option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")
if options.expHeader.upper() != "TRUE" and options.expHeader.upper() != "FALSE":
    sys.exit("error: -H --expHeader option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")
if options.normDiffana.upper() != "TRUE" and options.normDiffana.upper() != "FALSE":
    sys.exit("error: -N --normDiffana option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")
    
##############################################################   
command = ["/bin/bash -c '("]   
    
if options.buildContrast.upper() == "TRUE":
    # build of the command for the buildContrast script running
    commandBuilContrast = []
    commandBuilContrast.append("./buildContrast.R")
    commandBuilContrast.append(options.designFile)
    commandBuilContrast.append("\""+options.model+"\"")
    commandBuilContrast.append(options.comparisonFile)
    commandBuilContrast.append("{0}-contrastFile.tsv".format(options.projectName))
    commandBuilContrast = " ".join(commandBuilContrast)
    command.append(commandBuilContrast)
    command.append(";")
        
if options.normDiffana.upper() == "TRUE":
    # build of the command for the normDiffana script running   
    commandNormDiffana = []
    commandNormDiffana.append("./normDiffana.R")
    commandNormDiffana.append(options.normFig)
    commandNormDiffana.append(options.diffana)
    commandNormDiffana.append(options.diffanaFig)
    if options.contrast.upper() == "TRUE":
        commandNormDiffana.append("{0}-contrastFile.tsv".format(options.projectName))
    else:
        commandNormDiffana.append(options.contrast.upper())
    commandNormDiffana.append(options.designFile)
    commandNormDiffana.append("\""+options.model+"\"")
    commandNormDiffana.append(options.projectName)
    commandNormDiffana.append(options.expHeader)
    commandNormDiffana = " ".join(commandNormDiffana)
    command.append(commandNormDiffana)
    command.append(") &> ")
    command.append("{0}-deseq2.log".format(options.projectName))
    command.append("'")
    command = " ".join(command)

os.system(command)
