################################################################################
# NOTE: We use python3 because astropy is moving away from python2
# Build by running the command:
# docker build -t riptide .
################################################################################
FROM ubuntu:16.04

RUN echo 'deb http://us.archive.ubuntu.com/ubuntu trusty main multiverse' >> /etc/apt/sources.list

# A well-formed RUN instruction to install ubuntu packages, following
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    python3 \
    python3-dev \
    python3-pip \
    python3-tk \
 && rm -rf /var/lib/apt/lists/*

##### Python Packages #####
RUN pip3 install \
    astropy \
    h5py \
    ipython \
    matplotlib \
    numpy \
    pandas \
    pyyaml

##### riptide #####
ARG RIPTIDE_REPOSITORY=https://vmorello@bitbucket.org/vmorello/riptide.git

RUN mkdir software \
 && cd software \
 && git clone ${RIPTIDE_REPOSITORY} \
 && cd /software/riptide/c_src \
 && make all

ENV PYTHONPATH /software/riptide:${PYTHONPATH}
