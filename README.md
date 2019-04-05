ResFinder documentation
=============

The ResFinder service contains one perl script *resfinder.pl* which is the
script of the latest version of the ResFinder service. ResFinder identifies
acquired antimicrobial resistance genes in total or partial sequenced isolates
of bacteria.

This repository also contains a python script *resfinder.py* which is  a new version 
of ResFinder, but not yet running on the CGE server. This program was added because
it uses a newer version of blastn,  which, in contrary from the blastall version 
that the perl script uses, is avail to download.

## Content of the repository
1. resfinder.pl - the program
2. resfinder.py - (same program using an available blastn version - blastn-2.2.26+)
3. test.fsa     - test fasta file

## Installation

Setting up ResFinder script and database
```bash
# Go to wanted location for resfinder
cd /path/to/some/dir

# Clone and enter the resfinder directory
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd resfinder


# Installing up the ResFinder database
# Go to wanted location for resfinder database
cd /path/to/some/dir

# Clone and enter the resfinder directory
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
cd resfinder_db

```
Build Docker container
```bash
# Build container
docker build -t pmlst
```
## Download and install pMLST database
```bash
# Go to the directory where you want to store the ResFinder database
cd /path/to/some/dir
# Clone database from git repository (develop branch)
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
cd resfinder_db
ResFinder_DB=$(pwd)
# Install ResFinder  database with executable kma_index program
python3 INSTALL.py kma_index
```
If kma_index has not bin install please install kma_index from the kma repository: https://bitbucket.org/genomicepidemiology/kma

## Usage 

The BlastAll and FormatDB that the perl script uses are no longer available 
for downloading through ncbi. Therefore we have provided the resfinder.py 
script that uses Blastn instead. Note, this is not not script that is running 
on the CGE server. The CGE server is running the perl script using BlastAll.
The python script can be used through the Dockerfile.

The program can be invoked with the -h option to get help and more information of the service. Run Docker container
   
```bash
# Run example container
docker run --rm -it \
       -v $ResFinder_DB:/database \
       -v $(pwd):/workdir \
       resfinder -i [INPUTFILE] -1 [FASTQ1] -2 [FASTQ2] -o . [-b] [-k] -p [DB_PATH] [-t] [-l]

# Run resfinder container example
docker run --rm -it \
       -v $ResFinder_DB:/database \
       -v $(pwd):/workdir \
       resfinder -i test.fsa -o . -p /path/to/resfinder_db -b blastn -d aminoglycoside -t 0.90 -l 0.60
```
When running the docker file you have to mount 2 directory: 
 1. pmlst_db (pMLST database) downloaded from bitbucket
 2. An output/input folder from where the input file can be reached and an output files can be saved. 
Here we mount the current working directory (using $pwd) and use this as the output directory, 
the input file should be reachable from this directory as well.

Optional arguments:

  `-h, --help     show this help message and exit`
  `-i INPUTFILE, --inputfile inputfile      Input file`
  `-1 FASTQ1, --fastq1 FASTQ1     Raw read data file 1.`
  `-2 FASTQ2, --fastq2 FASTQ2      Raw read data file 2 (only required if data is paired-end).`
  `-o OUT_PATH, --outputPath OUT_PATH      Path to blast output`
  `-b BLAST_PATH, --blastPath BLAST_PATH      Path to blast`
  `-p DB_PATH, --databasePath DB_PATH     Path to the databases`
  `-k KMA_PATH, --kmaPath KMA_PATH      Path to KMA`
  `-q DB_PATH_KMA, --databasePathKMA DB_PATH_KMA      Path to the directories containing the KMA indexed databases. Defaults to the directory 'kma_indexing' inside the databasePath directory.`
  `-d DATABASES, --databases DATABASES      Databases chosen to search in - if none are specified all are used`
  `-l MIN_COV, --min_cov MIN_COV      Minimum coverage default 0.6`
  `-t THRESHOLD, --threshold THRESHOLD      Blast threshold for identity default minimum 0.9 `


### Web-server

A webserver implementing the methods is available at the [CGE 
website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/ResFinder/

### Documentation

The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/resfinder/overview.


Citation
=======

When using the method please cite:

Identification of acquired antimicrobial resistance genes.
Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup 
FM, Larsen MV.
J Antimicrob Chemother. 2012 Jul 10.
PMID: 22782487         doi: 10.1093/jac/dks261
[Epub ahead of print]


License
=======

Copyright (c) 2014, Ole Lund, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
