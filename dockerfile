FROM debian:9.4-slim

ENV DEBIAN_FRONTEND noninteractive

# Setup .bashrc file for convenience during debugging
RUN echo "alias ls='ls -h --color=tty'\n"\
"alias ll='ls -lrt'\n"\
"alias l='less'\n"\
"alias du='du -hP --max-depth=1'\n"\
"alias cwd='readlink -f .'\n"\
"PATH=$PATH\n">> ~/.bashrc

RUN set -ex && \
    # Basix setup \
    apt-get update -y -qq && \
    apt-get install -y -qq git \
    apt-utils \
    wget \
    python3-pip \
    libz-dev \
    && \
    # Python 3 setup \
    pip3 install --upgrade pip setuptools && \
    ln -sf /usr/bin/pip3 /usr/bin/pip && \
    ln -sf /usr/bin/python3 /usr/bin/python && \
    # Install python dependencies \
    pip install cython tabulate numpy biopython && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*

# Install BLAST
RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar -zxvf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    ln -s /ncbi-blast-2.7.1+/bin/blastn /usr/bin/blastn && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz

# Install ResFinder and databases
RUN pip install cgecore && \
    git clone -b 4.0 https://bitbucket.org/genomicepidemiology/resfinder.git && \
    mkdir resfinder/db_resfinder && \
    mkdir resfinder/db_pointfinder && \
    chmod a+x resfinder/run_resfinder.py && \
    ln -s /resfinder/run_resfinder.py /usr/bin/resfinder && \
    # Install KMA \
    apt-get install -y -qq libz-dev && \
    git clone --branch 1.0.1 https://bitbucket.org/genomicepidemiology/kma.git resfinder/cge/kma && \
    cd resfinder/cge/kma && make && cd ../.. && \
    ln -s /resfinder/cge/kma/kma /usr/bin/kma && \
    ln -s /resfinder/cge/kma/kma_index /usr/bin/kma_index && \
    apt-get clean && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*

ENV DEBIAN_FRONTEND Teletype

WORKDIR /workdir

ENTRYPOINT ["/resfinder/run_resfinder.py"]
