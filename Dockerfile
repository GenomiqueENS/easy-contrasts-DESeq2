############################################################
# Dockerfile to build Eoulsan container images with DEseq2
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu:12.04

# File Author / Maintainer
MAINTAINER Bauquet Xavier <xavier.bauquet@gmail.com>

# Add CRAN source to apt
RUN echo "deb http://cran.r-project.org/bin/linux/ubuntu precise/" > /etc/apt/sources.list.d/cran.list

# Add CRAN apt key
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Update the repository sources list
RUN apt-get update

############################################################
# Install DEseq2 and the dependencies 
############################################################
# Install Latex
RUN apt-get install --yes texlive-latex-base texlive-latex-extra 

# Install libxml2 
RUN apt-get install --yes libxml2-dev

# Install R 3.2.0
RUN apt-get install --yes r-base=3.2.0-4precise0 r-base-core=3.2.0-4precise0 r-base-dev=3.2.0-4precise0 r-base-html=3.2.0-4precise0 r-doc-html=3.2.0-4precise0 r-recommended=3.2.0-4precise0 r-cran-boot=1.3-15-1precise0 r-cran-class=7.3-12-1precise0 r-cran-cluster=2.0.1-1precise0 r-cran-codetools=0.2-11-1cran1precise0 r-cran-foreign=0.8.63-1precise0 r-cran-kernsmooth=2.23-14-1precise0 r-cran-lattice=0.20-31-1precise0 r-cran-mass=7.3-40-1precise0 r-cran-matrix=1.2-0-1precise0 r-cran-mgcv=1.8-6-1cran1precise0 r-cran-nlme=3.1.120-1precise0 r-cran-nnet=7.3-9-1precise0 r-cran-rpart=4.1-9-1precise0 r-cran-spatial=7.3-6-1precise0 r-cran-survival=2.38-1-1precise0

# Set CRAN repository to use
RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})' > ~/.Rprofile

# DEseq2 1.8.1 and dependencies 
RUN apt-get install --yes wget



# lme4 dependencies
RUN wget http://cran.r-project.org/src/contrib/Rcpp_0.11.6.tar.gz \
         http://cran.r-project.org/src/contrib/minqa_1.2.4.tar.gz \
         http://cran.r-project.org/src/contrib/nloptr_1.0.4.tar.gz \
		 http://cran.r-project.org/src/contrib/RcppEigen_0.3.2.4.0.tar.gz \
		 http://cran.r-project.org/src/contrib/lme4_1.1-7.tar.gz \
		 http://cran.r-project.org/src/contrib/SparseM_1.6.tar.gz \
		 http://cran.r-project.org/src/contrib/pbkrtest_0.4-2.tar.gz \
		 http://cran.r-project.org/src/contrib/quantreg_5.11.tar.gz \
		 http://cran.r-project.org/src/contrib/car_2.0-25.tar.gz \
		 http://cran.r-project.org/src/contrib/ellipse_0.3-8.tar.gz \
		 http://cran.r-project.org/src/contrib/scatterplot3d_0.3-35.tar.gz \
		 http://cran.r-project.org/src/contrib/leaps_2.9.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/S4Vectors_0.6.0.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/IRanges_2.2.3.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/XVector_0.8.0.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/GenomeInfoDb_1.4.0.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/GenomicRanges_1.20.4.tar.gz \
		 http://cran.r-project.org/src/contrib/RcppArmadillo_0.5.200.1.0.tar.gz \
		 http://bioconductor.org/packages/release/bioc/src/contrib/BiocGenerics_0.14.0.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/Biobase_2.28.0.tar.gz \
		 http://cran.r-project.org/src/contrib/lambda.r_1.1.7.tar.gz \
		 http://cran.r-project.org/src/contrib/futile.options_1.0.0.tar.gz \
		 http://cran.r-project.org/src/contrib/futile.logger_1.4.1.tar.gz \
		 http://cran.r-project.org/src/contrib/snow_0.3-13.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/BiocParallel_1.2.2.tar.gz \
		 http://cran.r-project.org/src/contrib/DBI_0.3.1.tar.gz \
		 http://cran.r-project.org/src/contrib/RSQLite_1.0.0.tar.gz \
		 http://bioconductor.org/packages/release/bioc/src/contrib/AnnotationDbi_1.30.1.tar.gz \
		 http://cran.r-project.org/src/contrib/XML_3.98-1.2.tar.gz \
		 http://cran.r-project.org/src/contrib/xtable_1.7-4.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/annotate_1.46.0.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/genefilter_1.50.0.tar.gz \
		 http://cran.r-project.org/src/contrib/locfit_1.5-9.1.tar.gz \
		 http://cran.r-project.org/src/contrib/RColorBrewer_1.1-2.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/geneplotter_1.46.0.tar.gz \
		 http://cran.r-project.org/src/contrib/plyr_1.8.2.tar.gz \
		 http://cran.r-project.org/src/contrib/digest_0.6.8.tar.gz \
		 http://cran.r-project.org/src/contrib/gtable_0.1.2.tar.gz \
		 http://cran.r-project.org/src/contrib/stringr_1.0.0.tar.gz \
		 http://cran.r-project.org/src/contrib/stringi_0.4-1.tar.gz \
		 http://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz \
		 http://cran.r-project.org/src/contrib/reshape2_1.4.1.tar.gz \
		 http://cran.r-project.org/src/contrib/dichromat_2.0-0.tar.gz \
		 http://cran.r-project.org/src/contrib/colorspace_1.2-6.tar.gz \
		 http://cran.r-project.org/src/contrib/munsell_0.4.2.tar.gz \
		 http://cran.r-project.org/src/contrib/labeling_0.3.tar.gz \
		 http://cran.r-project.org/src/contrib/scales_0.2.4.tar.gz \
		 http://cran.r-project.org/src/contrib/proto_0.3-10.tar.gz \
		 http://cran.r-project.org/src/contrib/ggplot2_1.0.1.tar.gz \
		 http://cran.r-project.org/src/contrib/Formula_1.2-1.tar.gz \
		 http://cran.r-project.org/src/contrib/latticeExtra_0.6-26.tar.gz \
		 http://cran.r-project.org/src/contrib/acepack_1.3-3.3.tar.gz \
		 http://cran.r-project.org/src/contrib/gridExtra_0.9.1.tar.gz \
		 http://cran.r-project.org/src/contrib/Hmisc_3.16-0.tar.gz \
		 http://cran.r-project.org/src/contrib/Archive/FactoMineR/FactoMineR_1.28.tar.gz \
		 http://www.bioconductor.org/packages/release/bioc/src/contrib/DESeq2_1.8.1.tar.gz

RUN R CMD INSTALL Rcpp_*.tar.gz && \
    R CMD INSTALL minqa_*.tar.gz && \
    R CMD INSTALL nloptr_*.tar.gz && \
	R CMD INSTALL RcppEigen_*.tar.gz && \
	R CMD INSTALL lme4_*.tar.gz && \
	R CMD INSTALL SparseM_*.tar.gz && \
	R CMD INSTALL pbkrtest_*.tar.gz && \
	R CMD INSTALL quantreg_*.tar.gz && \
	R CMD INSTALL car_*.tar.gz && \
	R CMD INSTALL ellipse_*.tar.gz && \
	R CMD INSTALL scatterplot3d_*.tar.gz && \
	R CMD INSTALL leaps_*.tar.gz && \
	R CMD INSTALL BiocGenerics_*.tar.gz && \
	R CMD INSTALL S4Vectors_*.tar.gz && \
	R CMD INSTALL IRanges_*.tar.gz && \
	R CMD INSTALL XVector_*.tar.gz && \
	R CMD INSTALL GenomeInfoDb_*.tar.gz && \
	R CMD INSTALL GenomicRanges_*.tar.gz && \
	R CMD INSTALL RcppArmadillo_*.tar.gz && \
	R CMD INSTALL Biobase_*.tar.gz && \
	R CMD INSTALL lambda.r_*.tar.gz && \
	R CMD INSTALL futile.options_*.tar.gz && \
	R CMD INSTALL futile.logger_*.tar.gz && \
	R CMD INSTALL snow_*.tar.gz && \
	R CMD INSTALL BiocParallel_*.tar.gz && \
	R CMD INSTALL DBI_*.tar.gz && \
	R CMD INSTALL RSQLite_*.tar.gz && \
	R CMD INSTALL AnnotationDbi_*.tar.gz && \
	R CMD INSTALL XML_*.tar.gz && \
	R CMD INSTALL xtable_*.tar.gz && \
	R CMD INSTALL annotate_*.tar.gz && \
	R CMD INSTALL genefilter_*.tar.gz && \
	R CMD INSTALL locfit_*.tar.gz && \
	R CMD INSTALL RColorBrewer_*.tar.gz && \
	R CMD INSTALL geneplotter_*.tar.gz && \
	R CMD INSTALL plyr_*.tar.gz && \
	R CMD INSTALL digest_*.tar.gz && \
	R CMD INSTALL gtable_*.tar.gz && \
	R CMD INSTALL stringi_*.tar.gz && \
	R CMD INSTALL magrittr_*.tar.gz && \
	R CMD INSTALL stringr_*.tar.gz && \
	R CMD INSTALL reshape2_*.tar.gz && \
	R CMD INSTALL dichromat_*.tar.gz && \
	R CMD INSTALL colorspace_*.tar.gz && \
	R CMD INSTALL munsell_*.tar.gz && \
	R CMD INSTALL labeling_*.tar.gz && \
	R CMD INSTALL scales_*.tar.gz && \
	R CMD INSTALL proto_*.tar.gz && \
	R CMD INSTALL ggplot2_*.tar.gz && \
	R CMD INSTALL Formula_*.tar.gz && \
	R CMD INSTALL latticeExtra_*.tar.gz && \
	R CMD INSTALL acepack_*.tar.gz && \
	R CMD INSTALL gridExtra_*.tar.gz && \
	R CMD INSTALL Hmisc_*.tar.gz && \
	R CMD INSTALL FactoMineR_*.tar.gz && \
	R CMD INSTALL DESeq2_*.tar.gz

# Remove the tar.gz files
RUN rm *.tar.gz



