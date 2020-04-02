# Description

ResFinder identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.

# Dependencies
ResFinder uses:  

* [biopython==1.7.6](https://pypi.org/project/biopython/)
* [cgecore==1.5.2](https://pypi.org/project/cgecore/).
* The newest version of [blastn](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* The newest version of [KMA](https://bitbucket.org/genomicepidemiology/kma/src/master/)

# Installation

**Requirements**, below software is required to install resfinder

* [docker](https://docs.docker.com/install/)
* [docker-compose](https://docs.docker.com/compose/install/) - recomended

```bash
# Clone repository
$ git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git --recursive
$ cd resfinder

# Chagne ownershipt on results directory - docker volume trick!
$ sudo chown 9999:9999 results

# Build docker image with docker-compose or manualy
$ docker-compose build
$ docker build -t genomicepidemiology/resfinder:v3 -f docker/Dockerfile .
```

# Program API

```
$ resfinder -h
usage: resfinder.py [-h] -i INPUTFILE [INPUTFILE ...] [-o OUT_PATH]
                    [-tmp TMP_DIR] [-mp METHOD_PATH] [-p DB_PATH]
                    [-d DATABASES] [-l MIN_COV] [-t THRESHOLD]
                    [-ao ACQ_OVERLAP] [-matrix] [-x] [-q]

optional arguments:
  -h, --help            show this help message and exit.
  -i INPUTFILE [INPUTFILE ...], --inputfile INPUTFILE [INPUTFILE ...] FASTA or FASTQ input files.
  -o OUT_PATH, --outputPath OUT_PATH Path to blast output.
  -tmp TMP_DIR, --tmp_dir TMP_DIR Temporary directory for storage of the results from the external software.
  -mp METHOD_PATH, --methodPath METHOD_PATH Path to method to use (kma or blastn).
  -p DB_PATH, --databasePath DB_PATH Path to the databases.
  -d DATABASES, --databases DATABASES Databases chosen to search in - if none are specified all are used.
  -l MIN_COV, --min_cov MIN_COV Minimum coverage.
  -t THRESHOLD, --threshold THRESHOLD Blast threshold for identity.
  -ao ACQ_OVERLAP, --acq_overlap ACQ_OVERLAP Genes are allowed to overlap this number of nucleotides. Default: 30.
  -matrix, --matrix     Gives the counts all all called bases at each position in each mapped template. Columns are: reference base, A count, C count, G count, T count, N count, - count.
  -x, --extented_output Give extented output with allignment files, template and query hits in fasta and a tab seperated file with gene profile results.
  -j, --json JSON is program's default output
  -q, --quiet
```

# Usage

```bash
# One time run
$ docker run --rm genomicepidemiology/resfinder:v3 resfinder -i /app/tests/test.fsa -mp /bin/blastn -t 0.90 -l 0.60 -p /app/db -o /app/results -d aminoglycoside,beta-lactam

# Multiple time run
$ docker-compose up -d
$ docker exec -it resfinder bash
$ resfinder -i /app/tests/test.fsa -p /app/db -mp /bin/blastn -t 0.90 -l 0.60 -o /app/results -x
```

# Web

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here: https://cge.cbs.dtu.dk/services/ResFinder/

# Citation

When using the method please cite:

**Identification of acquired antimicrobial resistance genes.**
*Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV.J*
Antimicrob Chemother. 2012 Jul 10.
PMID: 22782487 doi: 10.1093/jac/dks261
[Epub ahead of print]

# License

Copyright (c) 2014, Ole Lund, Technical University of Denmark All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0) Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
