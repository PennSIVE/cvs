#!/bin/bash

set -e

R_VERSION_MAJOR=4
R_VERSION_MINOR=1
R_VERSION_PATCH=2
CONFIGURE_OPTIONS="--with-cairo --with-jpeglib --enable-R-shlib --with-blas --with-lapack"

docker pull repronim/neurodocker:master
docker run --rm repronim/neurodocker:master generate docker \
    --pkg-manager apt \
    --base debian:buster \
    --fsl version=6.0.3 \
    --run "apt-get update && apt-get install -y \
            gfortran \
            git \
            g++ \
            libreadline-dev \
            libx11-dev \
            libxt-dev \
            libpng-dev \
            libjpeg-dev \
            libcairo2-dev \   
            libssl-dev \ 
            libxml2-dev \
            libudunits2-dev \
            libgdal-dev \
            libbz2-dev \
            libzstd-dev \
            liblzma-dev \
            libpcre2-dev \
            locales \
            screen \
            texinfo \
            texlive \
            texlive-fonts-extra \
            vim \
            wget \
            xvfb \
            tcl8.6-dev \
            tk8.6-dev \
            cmake \
            curl \
            unzip \
            libcurl4-gnutls-dev \
            libgsl-dev \
            libcgal-dev \
            libglu1-mesa-dev \
            libglu1-mesa-dev \
            libtiff5-dev \
            multiarch-support \
            && wget https://cran.rstudio.com/src/base/R-${R_VERSION_MAJOR}/R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz \
            && tar zxvf R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz \
            && rm R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH}.tar.gz \
            && cd /R-${R_VERSION_MAJOR}.${R_VERSION_MINOR}.${R_VERSION_PATCH} \
            && ./configure ${CONFIGURE_OPTIONS} \ 
            && make \
            && make install" \
            --afni version=latest \
            --ants version=2.3.4 \
    | docker build -t cvs_base -

docker build -t pennsive/cvs:latest .
id=$(docker images --format="{{.Repository}}:{{.Tag}} {{.ID}}" | egrep "^pennsive/cvs:latest " | cut -d' ' -f2)
version=$(cat version)
docker tag $id pennsive/cvs:$version
docker push pennsive/cvs:latest
docker push pennsive/cvs:$version
