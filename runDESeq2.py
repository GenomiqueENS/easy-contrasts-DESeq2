#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys
from optparse import OptionParser, OptionGroup

##############################################################
# runDESeq2
#
# Script to run and formatte options for both scripts
# buildContrast.R and normDiffana.R
#
# Version 1.6 (06/18/2015)
#
# Author: Xavier Bauquet
##############################################################

def createParser():
	""" Create a parser object. """

	# Create parser object
	parser = OptionParser()

	# Compulsory options group
	group = OptionGroup(parser, "Compulsory option")
	group.add_option("-m", "--model", type = "string", dest = "model", help ="The DESeq2 model to be used for the differential analysis as a formula: ~condition1+condition2+condition1:condition2 for more information about the formula please refer to the DESeq2 documentation http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html")
	parser.add_option_group(group)

	# Options for analysis with contrast vectors
	group = OptionGroup(parser, "Option for analysis with contrast vectors:")
	group.add_option("-c", "--contrast", type = "string", dest = "contrast", default = "FALSE", help ="logical TRUE/FALSE to escape or not the use of contrast vectors for differential anlysis. Default = FALSE")
	group.add_option("-b", "--buildContrast", type = "string", dest = "buildContrast", default = "FALSE", help ="logical TRUE/FALSE to escape or not the building of contrast vectors. Default = FALSE")
	parser.add_option_group(group)

	# Files options
	group = OptionGroup(parser, "Files options")
	group.add_option("-f", "--designFile", type = "string", dest = "designFile", default = "deseqDesign.txt", help ="Name of the design file\nDefault = deseqDesign.txt")
	group.add_option("-C", "--comparisonFile", type = "string", dest = "comparisonFile", default = "comparisonFile.txt", help ="Name of the file including the comparison to be computed in contrast vector. Comparisons must have this form: condition1%condition2_vs_condition3%condition4 where the % mean interaction between conditions and _vs_ mean a comparison between interactions. Default = comparisonFile.txt")
	group.add_option("--contrastFile", type = "string", dest = "contrastFile", default = "", help ="contrast file including the contrast vectors for the differential analysis.  Default = contrastFile.txt")
	parser.add_option_group(group)

	# Other options
	group = OptionGroup(parser, "Other options")
	group.add_option("-n", "--normFig", type = "string", dest = "normFig", default = "TRUE", help ="logical TRUE/FALSE to escape or not figures from the normalisation part. Default = TRUE")
	group.add_option("-D", "--diffana", type = "string", dest = "diffana", default = "TRUE", help ="logical TRUE/FALSE to escape or not the differential analysis. Default = TRUE")
	group.add_option("-d", "--diffanaFig", type = "string", dest = "diffanaFig", default = "TRUE", help ="logical TRUE/FALSE to escape or not figures from the differential analysis. Default = TRUE")
	group.add_option("-p", "--projectName", type = "string", dest = "projectName", default = "exp1", help ="Name of the project. Default = exp1")
	group.add_option("-H", "--expHeader", type = "string", dest = "expHeader", default = "TRUE", help ="logical TRUE/FALSE, TRUE if expression files have a header. Default = TRUE")
	group.add_option("-N", "--normDiffana", type = "string", dest = "normDiffana", default = "TRUE", help ="logical TRUE/FALSE, FALSE to escape the normalisation and differential analysis steps. Default = TRUE")
	group.add_option("--sizeFactorsType", type = "string", dest = "sizeFactorType", default = "ratio", help ="Type of size factor estimation to do. Two values possible: ratio or iterate.  Default = ratio")
	group.add_option("--fitType", type = "string", dest = "fitType", default = "parametric", help ="Fit type to be used in the dispersions estimation. Three values possible: parametric, local or mean.  Default = parametric")
	group.add_option("--statisticTest", type = "string", dest = "statisticTest", default = "Wald", help ="Satistic test to be used for the differential analysis. Two values possible: Wald or LRT.  Default = Wald")
	group.add_option("--prefix", type = "string", dest = "prefix", default = "", help ="Prefix for each file output by the module.  Default = """)

	parser.add_option_group(group)

	# Return parser object
	return parser

# Parse arguments
options, arguments = createParser().parse_args()

##############################################################
# Check command line arguments
if not options.model:
    sys.exit("Model option is missing please use -m or --model")

if options.contrast.upper() != "TRUE" and options.contrast.upper() != "FALSE":
   sys.exit("error: -c --contrast option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")

if options.buildContrast.upper() != "TRUE" and options.buildContrast.upper() != "FALSE":
   sys.exit("error: -b --buildContrast option only accepts 2 values TRUE or FALSE. Please use the -h option and refer to the help section")

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

if options.sizeFactorType != "ratio" and options.sizeFactorType != "iterate":
    sys.exit("error: --sizeFactorType option only accepts 2 values ratio or iterate. Please use the -h option and refer to the help section")

if options.fitType != "parametric" and options.fitType != "local" and options.fitType != "mean":
    sys.exit("error: --fitType option only accepts 3 values parametric, local or mean. Please use the -h option and refer to the help section")

if options.statisticTest != "Wald" and options.statisticTest != "LRT":
    sys.exit("error: --statisticTest option only accepts 2 values Wald or LRT. Please use the -h option and refer to the help section")


##############################################################
# Build command line to execute
command = ["/bin/bash -c '("]

if options.buildContrast.upper() == "TRUE":

    # build of the command for the buildContrast script running
    commandBuilContrast = []
    commandBuilContrast.append("./buildContrast.R")
    commandBuilContrast.append(options.designFile)
    commandBuilContrast.append("\""+options.model.replace(' ','')+"\"")
    commandBuilContrast.append(options.comparisonFile)
    commandBuilContrast.append("{0}-contrastFile.tsv".format(options.projectName))
    commandBuilContrast.append(options.prefix)
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
	commandNormDiffana.append("\""+options.model.replace(' ','')+"\"")
	commandNormDiffana.append(options.projectName)
	commandNormDiffana.append(options.expHeader)
	commandNormDiffana.append(options.sizeFactorType)
	commandNormDiffana.append(options.fitType)
	commandNormDiffana.append(options.statisticTest)
	commandNormDiffana.append(options.contrastFile)
	commandBuilContrast.append(options.prefix)
	commandNormDiffana = " ".join(commandNormDiffana)
	command.append(commandNormDiffana)
	command.append(") &> ")
	command.append("{0}-deseq2.log".format(options.projectName))
	command.append("'")
	command = " ".join(command)

# Execute command line
os.system(command)

