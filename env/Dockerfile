FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

# Set environment variables
ENV PATH=/opt/tools/bin:$PATH \
    RETICULATE_MINICONDA_ENABLED=FALSE \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US \
    LC_ALL=en_US.UTF-8 \
    LC_CTYPE=C \
    LC_TIME=en_US.UTF-8 \
    LC_MONETARY=en_US.UTF-8 \
    LC_PAPER=en_US.UTF-8 \
    LC_MEASUREMENT=en_US.UTF-8 \
    R_LIBS_USER='/usr/local/lib/R/library' \
    R_VERSION=4.5.1

# Copy required files into the container
COPY ./info_to_install_2026_01_13.txt /package_version_url.txt
COPY ./install_R_packages.R /install_R_packages.R
COPY ./2026_01_13/* /packages_dir/
COPY ./Dockerfile /Dockerfile

RUN apt update --yes && \
    # ---------------------------
    # Install libraries
    # ---------------------------
    apt install --yes \
    libxml2-dev \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    libcurl4-openssl-dev \
    libssl-dev \
    wget && \
    # Related to R language
    apt install --yes \
    g++ \
    gcc-13 \
    gfortran \
    zlib1g-dev \
    zlib1g \
    libbz2-dev \
    libreadline-dev \
    libbz2-dev \
    liblzma-dev \
    libpango1.0-dev \
    build-essential \
    libcurl4-gnutls-dev \
    libssl-dev \
    libjpeg-dev \
    libjpeg-turbo8-dev \
    xorg-dev && \
    # Related to nloptr package
    apt install --yes cmake && \
    # Related to rmarkdown package
    wget -O pandoc.deb https://github.com/jgm/pandoc/releases/download/3.1.3/pandoc-3.1.3-1-amd64.deb && \
    dpkg -i pandoc.deb && \
    rm /pandoc.deb && \
    # ---------------------------
    # Installation of R language
    # ---------------------------
    cd /usr/local/src && \
    wget -O R-${R_VERSION}.tar.gz https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz && \
    tar xzvf R-${R_VERSION}.tar.gz && \
    rm R-${R_VERSION}.tar.gz && \
    cd R-${R_VERSION} && \
    bash ./configure \
    --prefix=/usr/local \
    --x-includes=/usr/include/X11 \
    --x-libraries=/usr/lib/X11 \
    --enable-R-shlib \
    --enable-prebuilt-html && \
    make && make install && \
    # ---------------------------
    # Installation of R packages
    # ---------------------------
    Rscript /install_R_packages.R /package_version_url.txt && \
    # ---------------------------
    # Cleaning
    # ---------------------------
    rm /install_R_packages.R && \
    rm /package_version_url.txt && \
    rm -R /packages_dir && \
    apt autoremove --purge --yes && \
    apt clean
