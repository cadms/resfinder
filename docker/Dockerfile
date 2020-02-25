FROM python:3.8.1-slim

RUN useradd -u 9999 exec

WORKDIR /app

COPY db db/
COPY src src/
COPY tests tests/
COPY results results/
COPY vendors/kma kma/
COPY requirements.txt .

RUN apt update
RUN apt install -y \
  gcc \
  libz-dev \
  make \
  wget \
  vim

# Install KMA
RUN cd kma && make
RUN chmod +x kma
RUN mv kma/kma /bin
RUN rm -rf /app/kma

# Install BLASTn
RUN wget ftp.ncbi.nlm.nih.gov/blast/executables/LATEST
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/$(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2)
RUN tar -zxf $(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2)
RUN mv $(echo $(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2) | sed 's/-x64-linux.tar.gz//')/bin/blastn /bin
RUN rm -rf LATEST $(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2) $(echo $(grep 'linux.tar.gz' LATEST -m 1 | grep -Eo '\".+\"' | cut -d'"' -f2) | sed 's/-x64-linux.tar.gz//')

# Install python modules
RUN pip3 install -r requirements.txt

# Install resfinder
RUN ln -s /app/src/resfinder.py /bin/resfinder

RUN chown -R exec:exec /app

VOLUME /app/results

USER exec

CMD ["/bin/sh"]
