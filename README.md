
![License: GNU GPL 2.1](https://img.shields.io/badge/License-GNU_GPL--v2.1-lightgray?logo=GPLv3&logoColor=white)
![Written in R](https://img.shields.io/badge/Written_in-R-dodgerblue?logo=R&logoColor=white)
[![Container: Docker](https://img.shields.io/badge/Docker-V2.0-cyan?logo=docker&logoColor=white)](LINK_TO_ADD)
![Version 2.0](https://img.shields.io/badge/Version-V2.0-green4)

# easy-contrasts-DESeq2

**Easy-contrast-DEseq2** is a module to analyse count data from bulk RNA-Sequencing experiments. It performs quality control, normalisation and differential analysis using the expression count files as input. This module uses the `DESeq2` bioconductor R package. You will find in the [Beginner's guide to using the DESeq2 package](http://www.bioconductor.org/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf) basic informations about `DESeq2` and how produce the expression files it requires.

This module is part of the [Eoulsan](https://github.com/genomicpariscentre/eoulsan) pipeline analysis.

## Overview

The module is made of 4 RMarkdown files:

* 01_normDiffana.Rmd
* 02_child_collapse_replicates.Rmd (02 in the diagram)
* 03_child_differential_expression.Rmd (03 in the diagram)
* 04_child_pairwise_comparison.Rmd (04 in the diagram)

<p align="center">
<img src="./pipeline_overview.jpg?raw=true" alt="pipeline_overview" width="70%"/>
</p>

Notes, for the diagram:

(1) Load more packages in the *Settings > Environment* section. Use the `package::function` nomenclature as much as possible for function traceability.<br>
(2) Write more custom functions in the *Settings > Custom functions* section.<br>
(3) Add more parameters in the RMarkdown file header, and parse them in the *Settings > Parameters* section.<br>
(4) Add more columns to the design file, but explain them in the *Data preparation > Load metadata* section.


The three child files are conditionally executed within the main document based on parameter values:

- The file `02_child_collapse_replicates.Rmd` is executed when there are at least two rows (samples) in the design file with the same value in the `RepTechGroup` column. It is executed only once.
- The file `03_child_differential_expression.Rmd` is executed when the input parameter `diffanaTest` is set to TRUE. It is executed only once.
- The file `04_child_pairwise_comparison.Rmd` is executed when:
    * **Simple mode**: There are at least two distinct values in the `Reference` column of the design file.
    * **Complex mode**: There are comparisons in the comparison file, that match the ones automatically identified based on the provided formula (`deseqModel`) and that lead to unique constrast vectors.<br>
    This file is executed as many times as there are pairwise comparisons.

## Use

Either render the main `01_normDiffana.Rmd` document within RStudio, or use the command line.

```{bash}
docker run \
-ti --rm \
-v /:/ \
-w $(readlink -f .) \
-u $(id -u):$(id -g) \
name_of_a_docker_image \
Rscript -e "rmarkdown::render(
    input = '01_normDiffana.Rmd',
    output_file = 'myProject.html',
    params = list(projectName = 'myProject',
                  moreParameter = parameterValue))"
```

## Input parameters

They are represented in orange (except `complexMode` which is internally determined) on the pipeline overview.

Parameter | Default value | Definition
----------|---------------|--------------
projectName       |  "Experiment1"                             |  name of the project, display in the titles (document, figures)
designPath        |  "deseq2_Experiment1-deseq2Design.txt"     |  path to the design file
comparisonPath    |  "deseq2_Experiment1-comparisonFile.txt"   |  path to the comparison file*
diffanaTest       |  TRUE                          |  Whether to perform the differential expression analysis or not
expHeader         |  TRUE                          |  Whether the design table has a header or not
deseqModel        |  "~Condition"                  |  DESeq2 model
sizeFactorType    |  "ratio"                       |  Size factor type
fitType           |  "parametric"                  |  Fit type
statisticTest     |  "Wald"                        |  Statistic test
weightContrast    |  FALSE                         |  Whether to weight the contrast vector by the number of samples or not*
prefix            |  "deseq2_"                     |  Prefix to save the output TSV files
plotInteractive   |  FALSE                         |  Whether to make interactive volcano plots or not
logoUrl           |  "logo-GenomiqueENS-90pxh.png" | Link to the logo displayed in the top left corner of the output document
authorName        |  "Eoulsan"                     |  Name of the author, for the bottom left corner of the document
authorMail        |  "eoulsan@biologie.ens.fr"     |  Mail of the author, for the bottom left corner of the document
leaveOnError      |  TRUE                          |  Whether to stop the rendering in case of error

*: The `comparisonPath` and `weightContrast` are considered only when `complexMode` is TRUE.

## Input files

**CAUTION**: All the input files should be tabulated files.

### Design and comparison files

The design file should include at least the following columns: *Name*, *Condition*, *RepTechGroup*, *Reference* and *expressionFiles*.

* **Name**: the names of your samples
* **Condition**: the biological replicates. All biological replicates should have the same **Condition** value
* **RepTechGroup**: the technical replicates. All technical replicates should have the same **RepTechGroup** value to be pooled together
* **Reference**: conditions to use as a reference in the differential expression analysis should have a **Reference** value set to 0. Conditions to ignore should have a negative value. Other conditions to be compared to the reference must have a positive value.

Notes:

- More columns can be used for the contrast mode (see the model design with column type and day).
- None of the values must starts by a digit or a symbol.

#### Example: simple mode

The **Reference** column should be filled:

- Reference < 0: the samples are considered for the normalisation step only but not by the differential expression analysis
- Reference = 0: the condition associated with these samples is considered as a reference in the differential expression analysis
- Reference > 0: the condition associated with these samples is compared to each condition having a lower reference value

Example:

Name    |Condition |RepTechGroup |Reference |expressionFiles         |Condition 
--------|----------|-------------|----------|------------------------|----------
sample1 |WT-day1a  |WT-day1      |0         |expression_WT-day1a.tsv |WT        
sample2 |WT-day1b  |WT-day1      |0         |expression_WT-day1b.tsv |WT        
sample3 |KO-day1a  |KO-day1      |1         |expression_KO-day1a.tsv |KO        
sample4 |KO-day1b  |KO-day1      |1         |expression_KO-day1b.tsv |KO        
sample5 |WT-day2a  |WT-day2      |0         |expression_WT-day2a.tsv |WT        
sample6 |WT-day2b  |WT-day2      |0         |expression_WT-day2b.tsv |WT        
sample7 |KO-day2a  |KO-day2      |2         |expression_KO-day2a.tsv |KO        
sample8 |KO-day2b  |KO-day2      |2         |expression_KO-day2b.tsv |KO        

The `deseqModel` parameter could be **"~ Condition"**.

#### Example: complex mode

The **Reference** column should NOT be filled.

Example:

Name    |Condition |RepTechGroup |expressionFiles         |type |day
--------|----------|-------------|------------------------|-----|---
sample1 |WT-day1a  |WT-day1      |expression_WT-day1a.tsv |WT   |d1
sample2 |WT-day1b  |WT-day1      |expression_WT-day1b.tsv |WT   |d1
sample3 |KO-day1a  |KO-day1      |expression_KO-day1a.tsv |KO   |d1
sample4 |KO-day1b  |KO-day1      |expression_KO-day1b.tsv |KO   |d1
sample5 |WT-day2a  |WT-day2      |expression_WT-day2a.tsv |WT   |d2
sample6 |WT-day2b  |WT-day2      |expression_WT-day2b.tsv |WT   |d2
sample7 |KO-day2a  |KO-day2      |expression_KO-day2a.tsv |KO   |d2
sample8 |KO-day2b  |KO-day2      |expression_KO-day2b.tsv |KO   |d2

The `deseqModel` parameter could be **"~type+day+type:day"**.

The comparison file (beind the `comparisonPath` parameter) could be:

```
comparison1    dayd1_vs_dayd2
comparison2    typeKO_vs_typeWT
comparison3    typeWT%dayd1_vs_typeWT%dayd2
comparison4    typeKO%dayd1_vs_typeWT%dayd1
comparison5    typeKO%dayd2_vs_typeWT%dayd2
```

Notes:

- The comparison file doesn't contain a header.
- The comparison file contains two columns: name of the comparison and comparison formula.
- The comparison formula is case sensitive. Respect the column names and values from the design file.

#### Count files

Count files should include a first column with features description (identifiers or symbols) and a second column with counts. These count files can include a header or not (`expHeader` parameter). For each sample, the file path is retrieved through the **expressionFile** column of the design file (`designPath` parameter).

Example:

```
Id      Count
ENST00000000233	569
ENST00000000412	119
ENST00000000442	25
ENST00000001008	42
ENST00000001146	0
ENST00000002125	19
...
```

## Output files

Rendering the main `01_normDiffana.Rmd` document output three types of files:

* the rendered HTML notebook ;
* various figures, saved in the `${prefix}_figures` folder, and highlighted in blue on the pipeline overview ;
* various TSV files, along with the HTML report, and highlighted in purple on the pipeline overview.

## Environment

This module has been developed using R version 4.5.1 and the following main packages:

Package                |  Version 
-----------------------|----------
Biobase                |  2.68.0
BiocGenerics           |  0.54.0
ComplexHeatmap         |  2.24.0
circlize               |  0.4.16
DESeq2                 |  1.48.0
dplyr                  |  1.1.4
FactoMineR             |  2.12
generics               |  0.1.3
ggdendro               |  0.2.0
ggplot2                |  4.0.0
ggrepel                |  0.9.6
GenomeInfoDb           |  1.44.0
GenomicRanges          |  1.60.0
gridtext               |  0.1.5
IRanges                |  2.42.0
kableExtra             |  1.4.0
knitr                  |  1.50
Matrix                 |  1.7-3
MatrixGenerics         |  1.20.0
matrixStats            |  1.5.0
plotly                 |  4.11.0
RColorBrewer           |  1.1-3
rmdformats             |  1.0.4
rmarkdown              |  2.29
reshape2               |  1.4.4
S4Vectors              |  0.46.0
stringr                |  1.5.2
SummarizedExperiment   |  1.38.0


A Docker container containing all the packages of interest along with their dependencies is available at **XXX**.
