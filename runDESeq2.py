#!/usr/bin/python
# -*- coding: utf-8 -*-

##############################################################
# gff to centil_gff
#
# Python script  for gene coverage
#
# Author: Xavier Bauquet
##############################################################
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--normFig", type = "string", dest = "normFig", default = "TRUE", help ="booleen TRUE/FALSE to escape or not figures from the normalisation part\nDefault = TRUE")
parser.add_option("-D", "--diffana", type = "string", dest = "diffana", default = "TRUE", help ="booleen TRUE/FALSE to escape or not the differential analysis\nDefault = TRUE")
parser.add_option("-d", "--diffanaFig", type = "string", dest = "diffanaFig", default = "TRUE", help ="booleen TRUE/FALSE to escape or not figures from the differential analysis\nDefault = TRUE")
parser.add_option("-c", "--contrast", type = "string", dest = "contrast", default = "FALSE", help ="Name of the file including contrats vector for the differential analysis, if FALSE the differential analysis will be performed without contrast vectors, comparing all samples each other or to the reference sample\nDefault = FALSE")
parser.add_option("-f", "--designFile", type = "string", dest = "designFile", default = "design.txt", help ="Name of the design file\nDefault = design.txt")
parser.add_option("-m", "--model", type = "string", dest = "model", help ="The DESeq2 model to be used for the differential analysis\nas a formula: ~condition1+condition2+condition1:condition2 for more information about the formula please refer to the DESeq2 documentation http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html")
parser.add_option("-C", "--comparisonFile", type = "string", dest = "comparisonFile", default = "comparisonFile.txt", help ="Name of the file including the comparison to be compute in contrast vector.\nComparisons must have this form: condition1%condition2_vs_condition3%condition4 where the % mean interaction between conditions and _vs_ mean a comparison between interactions\nDefault = comparisonFile.txt")
parser.add_option("-N", "--projectName", type = "string", dest = "projectName", default = "exp1", help ="Name of the project\nDefault = exp1")
parser.add_option("-v", "--verbose", type = "string", dest = "verbose", default = "TRUE", help ="Verbose mode\nDefault = TRUE")
parser.add_option("-b", "--buildContrast", type = "string", dest = "buildContrast", default = "FALSE", help ="booleen TRUE/FALSE to run the building of contrast vectors\nDefault = FALSE")
parser.add_option("-l", "--log", type = "string", dest = "log", default = "deseq2.log", help ="name of the Log file\nDefault = deseq2.log")
options, arguments = parser.parse_args()

if options.model:
    if options.buildContrast == "TRUE":
        commandBuilContrast = []
        commandBuilContrast.append("./buildContrast.R")
        commandBuilContrast.append(options.designFile)
        commandBuilContrast.append(options.model)
        commandBuilContrast.append(options.comparisonFile)
        commandBuilContrast.append(options.verbose)
        commandBuilContrast.append("&>")
        commandBuilContrast.append(options.log)
        commandBuilContrast = " ".join(commandBuilContrast)
        print("Start building of contrast vectors")
        print str(commandBuilContrast)
        os.system(commandBuilContrast)
        
    commandNormDiffana = []
    commandNormDiffana.append("./normDiffana.R")
    commandNormDiffana.append(options.normFig)
    commandNormDiffana.append(options.diffana)
    commandNormDiffana.append(options.diffanaFig)
    commandNormDiffana.append(options.contrast)
    commandNormDiffana.append(options.designFile)
    commandNormDiffana.append(options.model)
    commandNormDiffana.append(options.projectName)
    commandNormDiffana.append(options.verbose)
    if options.buildContrast == "TRUE":
        commandNormDiffana.append("&>>")
    else:
        commandNormDiffana.append("&>")
    commandNormDiffana.append(options.log)
    commandNormDiffana = " ".join(commandNormDiffana)
    print("Start normalisation and differential analysis")
    print str(commandNormDiffana)
    os.system(commandNormDiffana)