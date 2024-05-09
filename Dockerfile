FROM --platform=linux/x86_64 mambaorg/micromamba:0.24.0 as app

# ARG sets environment variables during the build stage
ARG LIFTOFF_VER="1.6.3"
ARG MINIMAP2_VER="2.24"
ARG SAMTOOLS_VER="1.20"

USER root
WORKDIR /

# LABEL instructions tag the image with metadata that might be important to the user
# Optional, but highly recommended
LABEL base.image="mambaorg/micromamba:0.24.0"
LABEL dockerfile.version="1"
LABEL software="liftoff"
LABEL software.version=$LIFTOFF_VER
LABEL description="A tool that accurately maps annotations in GFF or GTF between assemblies of the same, or closely-related species."
LABEL website="https://github.com/agshumate/Liftoff"
LABEL license="https://github.com/agshumate/Liftoff/blob/master/LICENSE.md"
LABEL maintainer="John Arnn"
LABEL maintainer.email="jarnn@utah.gov"

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    perl \
    zlib1g \
    libncurses5 \
    bzip2 \
    liblzma-dev \
    libcurl4-gnutls-dev \
    gnuplot \
    procps && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda  \
    minimap2=${MINIMAP2_VER} \
    python=3.9.1 \
    samtools=${SAMTOOLS_VER} \
    Liftoff=${LIFTOFF_VER} && \
    micromamba clean --all --yes
