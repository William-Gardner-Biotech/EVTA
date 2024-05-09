FROM --platform=linux/x86_64 mambaorg/micromamba:0.24.0 as app

# ARG sets environment variables during the build stage
ARG LIFTOFF_VER="1.6.3"
ARG MINIMAP2_VER="2.24"
ARG SAMTOOLS_VER="1.20"
ARG HTSLIBVER="1.12"
ARG IVARVER="1.4.2"

USER root
WORKDIR /

# LABEL instructions tag the image with metadata that might be important to the user
# Optional, but highly recommended
LABEL base.image="mambaorg/micromamba:0.24.0"
LABEL dockerfile.version="1"
LABEL software="liftoff"
LABEL maintainer="Will Gardner"
LABEL maintainer.email="wkgardner@wisc.edu"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    ca-certificates \
    procps \
    autoconf \
    autotools-dev \
    automake \
    python3 \
    bedtools \
    python3-pip \
    bwa \
    minimap2 \
    vim \
    openjdk-11-jdk \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev \
    libssl-dev \
    gcc \
    make \
    ca-certificates \
    perl \
    zlib1g \
    bzip2 \
    gnuplot \
    procps && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda  \
    minimap2=${MINIMAP2_VER} \
    python=3.9.1 \
    Liftoff=${LIFTOFF_VER} && \
    micromamba clean --all --yes

# Install BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz --no-check-certificate \
    && tar -xvzf BBMap_38.90.tar.gz \
    && rm BBMap_38.90.tar.gz \
    && mv bbmap /opt/
# Update PATH for BBTools
ENV PATH="/opt/bbmap:${PATH}"

ENV PATH=/opt/conda/bin:$PATH

# download, compile, and install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure && \
    make && \
    make install && \
    make test

# installing htslib
RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIBVER}/htslib-${HTSLIBVER}.tar.bz2 && \
    tar xvf htslib-${HTSLIBVER}.tar.bz2 && \
    rm htslib-${HTSLIBVER}.tar.bz2 && \
    cd htslib-${HTSLIBVER}/ && \
    ./configure && \
    make && \
    make install

RUN wget -q https://github.com/andersen-lab/ivar/archive/v${IVARVER}.tar.gz && \
    tar -xf v${IVARVER}.tar.gz && \
    rm -rf v${IVARVER}.tar.gz && \
    cd ivar-${IVARVER} && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    ldconfig

RUN wget -qO- https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# make sure shells are bash
CMD ["/bin/bash"]
