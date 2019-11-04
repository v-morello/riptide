################################################################################
# NOTE: We use python3 because astropy is moving away from python2
# Build by running the command:
# docker build -t riptide .
################################################################################
FROM python:3.7

# A well-formed RUN instruction to install ubuntu packages, following
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    vim \
 && rm -rf /var/lib/apt/lists/*

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

RUN pip install ipython

##### riptide #####
ARG RIPTIDE_REPOSITORY=https://vmorello@bitbucket.org/vmorello/riptide.git

# TODO FIXME (04/11/19): Don't checkout the pachaging branch once it has been merged
RUN mkdir software \
 && cd software \
 && git clone ${RIPTIDE_REPOSITORY} \
 && cd /software/riptide \
 && git checkout packaging \
 && make install

ENTRYPOINT [ "/bin/bash" ]