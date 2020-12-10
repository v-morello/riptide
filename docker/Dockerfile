FROM python:3.8

# A well-formed RUN instruction to install ubuntu packages, following
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    vim \
 && rm -rf /var/lib/apt/lists/*

# Additional python packages for convenience
RUN pip install ipython

# Install riptide
ENV RIPTIDE_PATH=/software/riptide
RUN mkdir -p ${RIPTIDE_PATH}
WORKDIR ${RIPTIDE_PATH}
COPY . ${RIPTIDE_PATH}
# NOTE: make clean is important, because the stale C++ build files may have been copied from the 
# host. This forces a rebuild of the C++ source.
RUN make clean install

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
    echo "set hlsearch" >> ~/.vimrc && \
    echo "syntax on" >> ~/.vimrc

# Run tests, this will deliberately make the build fail if there are any issues
RUN make tests
ENTRYPOINT [ "/bin/bash" ]