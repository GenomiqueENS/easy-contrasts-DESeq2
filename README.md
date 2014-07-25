easy-contrasts-DESeq2
=====================
**Easy-contrast-DEseq2** is a module for analysis of count data from RNA-seq. It performs both Normalization and Differential analysis using expression count matrix. This module uses the DESeq2 bioconductor R-package  and perform the construction of contrast vectors used by [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).

##Use:
Easy-contrast-DEseq2 can be used in 3 modes: 
To run the Easy-contrast-DEseq2 module you should make sure that the 3 scripts ( buildContrast.R, normDiffana.R,  runDESeq2.py ) are executable . For this open the Terminal and use the ls -l command: 
```
I have no name!@b6b0df808204:/root$ ls -l
-rwxrwxr-x 1 2743 users 10434 Jul 25 13:38 buildContrast.R
-rwxrwxr-x 1 2743 users 29189 Jul 25 13:38 normDiffana.R
-rwxrwxr-x 1 2743 users  6936 Jul 25 13:38 runDESeq2.py
I have no name!@b6b0df808204:/root$
```
You should have a 'x' at the end of the first expression for all the 3 files, as on this example. If you don't have this 'x' you must use the chmod command: 
```
I have no name!@b6b0df808204:/root$ chmod +x buildContrast.R normDiffana.R runDESeq2.py
I have no name!@b6b0df808204:/root$
```
Now you can run the Easy-contrast-DEseq2 with the following command: 
```
./runDESeq2.py -m '~Condition'
```

###Options:
All options available on the Easy-contrast-DEseq2 are presented here:
  * **-m --model** : deseqModel, the only compulsory option. This option should contain the deseq formula as on the previous example.
  * **-c --contrast** : TRUE/FALSE. If this option have the “TRUE” value the differential analysis will be performed using contrast vectors. **Default=** FALSE.
  * **-b --buildContrast** : TRUE/FALSE. If this option have the “TRUE” value the comparisonFile.tsv will be load and the buildContrast.R script will generate the 'projectName'-contrastMatrix.tsv file including the contrast vectors. **Default=** FALSE.
  * **-f --designFile** : the name of the design file. **Default=** deseqDesign.txt.
  * **-C --comparisonFile** : the name of the file including the comparison to be compute in contrast vector (see the Contrast file section). **Default=** comparisonFile.txt.
  * **-n --normFig** : TRUE/FALSE. If this option have the “FLASE” value figures from the normalization will be escaped. **Default=** TRUE.
  * **-N --normDiffana** : TRUE/FALSE. If this option have the “FLASE” value the normalization and differential analysis steps will be escaped. This option can be use to only build contrast vectors. **Default=** TRUE.
  * **-d --diffanaFig** :  TRUE/FALSE. If this option have the “FLASE” value figures from the differential analysis will be escaped. **Default=** TRUE.
  * **-D --diffana** : TRUE/FALSE. If this option have the “FLASE” value the differential analysis step will be escaped.**Default=** TRUE.
  * **-p --projectName** : The name of the project. **Default=**exp1.
  * **-H --expHeader** :  TRUE/FALSE. “TRUE” if the expression files have a header. **Default=** TRUE.
  * **-v --verbose** : TRUE/FALSE. Verbose mode. **Default=** TRUE.

###Classical mode:
The Classical mode performs the differential analysis comparing all biological replicates from the “Condition” design column each other.
To run the Easy-contrast-DEseq2 module on the Classical mode use this command: 
```
./runDESeq2.py -m '~Condition'
```
For this mode you don't need options -c and -b, and you don't need the comparison file.

###Reference mode: 
The Reference mode performs the differential analysis comparing all biological replicates from the “Condition” design column to the reference condition. You should specify the reference condition using the “true” value on the design “Reference” column for the reference condition (see the design file section).
To run the Easy-contrast-DEseq2 module on the Reference mode use this command: 
```
./runDESeq2.py -m '~Condition'
```
For this mode you don't need options -c and -b, and you don't need the comparison file.

###Contrast mode:
The Contrast mode  performs the differential analysis from the comparison file (see comparison file section). 
To run the Easy-contrast-DEseq2 module on the Contrast mode use this command: 
```
./runDESeq2.py -m '~type+day' -c TRUE -b TRUE
```

##Installation: 
This module was code using R version 3.1.0, DESeq2 1.4.5, and two other R packages RcolorBrewer 1.0.5 and FactoMineR 1.26. 
To use Easy-contrast-DEseq2 you should: 
  * Install R in the good version, and all packages
  * Use the docker file available on [genomicpariscentre/deseq2](https://registry.hub.docker.com/u/genomicpariscentre/deseq2/)
  * Use the Dockerfile present in  Easy-contrast-DEseq2 to install the docker directly on your computer

##Input files:
**CAUTION:** All file use by Easy-contrast-DEseq2 should be tabulated files.

###Expression files: 
Expression files should include a first column with names of the features (For example genes names, transcript ensembl id …) and a second column with counts. These expression files can include a header or not, but this information should be specified into the options of the Easy-contrast-DEseq2.

###Design file:
The design file should include at least the following column: Name, Condition, RepTechGroup, Reference and expressionFiles.
  * *Name:* the names of your samples
  * *Condition:* the biological replicates. All biological replicates should have the same condition name
  * *RepTechGroup:* the technical replicates. All technical replicates should have the same RepTechGroup name to by pooled during the normalization step
  * *Reference:* the reference condition use for the differential analysis. All sample should have the “false” value in the Reference column except the reference that should have the “true” value. You can use the “true” value on all sample of the same biological replicates group (Condition column) or only on one. The result will be the same.

###Comparison file:
The comparison file is used to generate the contrast vectors. It should include 2 columns: 
  * the name of the comparison
  * the formula of the comparison

The comparison file must have no header. The formula of the comparison is construct with the name of the column on the design file paste to the name of the condition, each “columncondition” should be separated by the “%” symbol to notify an association between  “columnconditions” and separated by the “_vs_” symbol to notify a comparison. 
Exemple: 
With the model design available on the Easy-contrast-DEseq2 github: 

Name    |Condition |RepTechGroup |Reference |expressionFiles         |type |day
--------|----------|-------------|----------|------------------------|-----|---
sample1 |WT-day1a  |WT-day1      |false     |expression_WT-day1a.tsv |WT     |1
sample2 |WT-day1b  |WT-day1      |false     |expression_WT-day1b.tsv |WT     |1
sample3 |KO-day1a  |KO-day1      |false     |expression_KO-day1a.tsv |KO     |1
sample4 |KO-day1b  |KO-day1      |false     |expression_KO-day1b.tsv |KO     |1
sample5 |WT-day2a  |WT-day2      |false     |expression_WT-day2a.tsv |WT     |2
sample6 |WT-day2b  |WT-day2      |false     |expression_WT-day2b.tsv |WT     |2
sample7 |KO-day2a  |KO-day2      |false     |expression_KO-day2a.tsv |KO     |2
sample8 |KO-day2b  |KO-day2      |false     |expression_KO-day2b.tsv |KO     |2

With the deseq model: 
```
~type+day
```
We want to compare WT at the day 1 to WT at the day 2, the comparison formula will be: 
```
typeWT%day1_vs_typeWT%day2
```
**CAUTION:** You have to respect the letter case as the one of the design file 