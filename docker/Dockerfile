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
    vim \
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

# Aliases
RUN echo "alias python='python3'" >> ~/.bashrc && \
    echo "alias pip='pip3'" >> ~/.bashrc

# Smart history search with arrow keys
RUN echo "\"\e[A\":history-search-backward" >> ~/.inputrc && \
    echo "\"\e[B\":history-search-forward" >> ~/.inputrc

# Configure vim to indent with 4 spaces and behave nicely in general
# https://stackoverflow.com/questions/234564/tab-key-4-spaces-and-auto-indent-after-curly-braces-in-vim
RUN echo "filetype plugin indent on" >> ~/.vimrc && \
    echo "set tabstop=4" >> ~/.vimrc && \
    echo "set shiftwidth=4" >> ~/.vimrc && \
    echo "set expandtab" >> ~/.vimrc && \
    echo "set pastetoggle=<F2>" >> ~/.vimrc && \
    echo "set hlsearch" >> ~/.vimrc

##### riptide #####
ARG RIPTIDE_REPOSITORY=https://vmorello@bitbucket.org/vmorello/riptide.git

RUN mkdir software \
 && cd software \
 && git clone ${RIPTIDE_REPOSITORY} \
 && cd /software/riptide/c_src \
 && CC=gcc make all

ENV PYTHONPATH /software/riptide:${PYTHONPATH}
