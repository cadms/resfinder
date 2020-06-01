### Install BLASTn ###
FROM debian:buster-slim as blast-builder

WORKDIR /app

RUN apt update
RUN apt install -y wget

RUN wget ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/$(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2)
RUN tar -zxf $(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2)
RUN mv $(echo $(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2) | sed 's/-x64-linux.tar.gz//')/bin/blastn /app

### Install KMA ###
FROM debian:buster-slim as kma-builder

WORKDIR /app

COPY vendors/kma /app

RUN apt update
RUN apt install -y \
  gcc \
  libz-dev \
  make

RUN make
RUN chmod +x kma

### Install ResFinder ###
FROM python:3.8-slim

ARG SERVICE_NAME=resfinder
ENV SERVICE_NAME=$SERVICE_NAME
ARG SERVICE_VERSION=v3
ENV SERVICE_VERSION=$SERVICE_VERSION

RUN useradd -u 9999 exec

WORKDIR /app

COPY db db/
COPY src src/
COPY tests tests/
COPY results results/
COPY requirements.txt requirements.txt
COPY --from=kma-builder /app/kma /bin
COPY --from=blast-builder /app/blastn /bin

RUN pip3 install -r requirements.txt
RUN ln -s /app/src/resfinder.py /bin/resfinder
RUN chown -R exec:exec /app

VOLUME /app/results /app/data

USER exec

CMD ["/bin/bash"]
