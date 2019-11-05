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

# Additional packages for convenience
RUN pip install ipython

# Install only requirements first, so that this step can be cached
ENV RIPTIDE_PATH=/software/riptide
RUN mkdir -p ${RIPTIDE_PATH}
WORKDIR ${RIPTIDE_PATH}

# Install riptide
COPY . ${RIPTIDE_PATH}
RUN make install
ENTRYPOINT [ "/bin/bash" ]