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
The installation described here will first install the actual ResFinder software, 
then the dependencies, and finally the databases. A more detailed breakdown of the 
installation is provided below:
1. Install ResFinder tool
2. Install python module BioPython
3. Install python module CGECore
4. Install BLAST (optional)
5. install KMA (optional)
6. Download ResFinder database
7. Download PointFinder database
8. Index databases with KMA (if installed)

### ResFinder tool
Setting up ResFinder script and database
```bash
# Go to wanted location for resfinder
cd /path/to/some/dir

# Clone branch 4.0 and enter the resfinder directory
git clone -b 4.0 https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd resfinder


# Installing up the ResFinder database
# Go to wanted location for resfinder database
cd /path/to/some/dir

# Clone and enter the resfinder directory
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
cd resfinder_db

```

### Dependencies:
Depending on how you plan to run ResFinder BLAST and KMA can be optional.
BLAST is used to analyse assemblies (ie. FASTA files).
KMA is used to analyse read data (ie. FASTQ files).

#### BioPython
To install BioPython you can use pip
```bash
pip3 install biopython
```
For more information visit the BioPython website
```url
http://biopython.org
```

#### CGECore
To install CGECore you can use pip
```bash
pip3 install cgecore
```
Source code is available from:
```url
https://bitbucket.org/genomicepidemiology/cge_core_module/src/master/
```

#### BLAST (optional)
If you don't want to specify the path of blastn every time you run
ResFinder, make sure that blastn is in you PATH.

Blastn can be obtained from:
```url
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

#### KMA (optional)
The instructions here will install KMA in the default location ResFinder uses. KMA 
can be installed in another location but the path to KMA will then need to be 
specified every time you run ResFinder.
```bash
git clone https://bitbucket.org/genomicepidemiology/kma.git
```

## Usage 

You can run resfinder command line using python3
   
```

# Example of running resfinder
python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db \
-b /path/to/blastn -d aminoglycoside -t 90.00 -l 0.60

# The program can be invoked with the -h option 
usage: run_resfinder.py [-h] [-ifa INPUTFASTA]
                        [-ifq INPUTFASTQ [INPUTFASTQ ...]] [-scripts SCRIPTS]
                        [-o OUT_PATH] [-b BLAST_PATH] [-k KMA_PATH]
                        [-s SPECIES] [-l MIN_COV] [-t THRESHOLD]
                        [-db_res DB_PATH_RES] [-db_res_kma DB_PATH_RES_KMA]
                        [-d DATABASES] [-acq] [-c] [-db_point DB_PATH_POINT]
                        [-g SPECIFIC_GENE [SPECIFIC_GENE ...]] [-u]

optional arguments:
  -h, --help            show this help message and exit
  -ifa INPUTFASTA, --inputfasta INPUTFASTA
                        Input fasta file.
  -ifq INPUTFASTQ [INPUTFASTQ ...], --inputfastq INPUTFASTQ [INPUTFASTQ ...]
                        Input fastq file(s). Assumed to be single-end fastq if
                        only one file is provided, and assumed to be paired-
                        end data if two files are provided.
  -o OUT_PATH, --outputPath OUT_PATH
                        All output will be stored in this directory.
  -b BLAST_PATH, --blastPath BLAST_PATH
                        Path to blastn
  -k KMA_PATH, --kmaPath KMA_PATH
                        Path to kma
  -s SPECIES, --species SPECIES
                        Species in the sample
  -l MIN_COV, --min_cov MIN_COV
                        Minimum (breadth-of) coverage
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for identity
  -db_res DB_PATH_RES, --db_path_res DB_PATH_RES
                        Path to the databases for ResFinder
  -db_res_kma DB_PATH_RES_KMA, --db_path_res_kma DB_PATH_RES_KMA
                        Path to the ResFinder databases indexed with KMA.
                        Defaults to the 'kma_indexing' directory inside the
                        given database directory.
  -d DATABASES, --databases DATABASES
                        Databases chosen to search in - if none is specified
                        all is used
  -acq, --acquired      Run resfinder for acquired resistance genes
  -c, --point           Run pointfinder for chromosomal mutations
  -db_point DB_PATH_POINT, --db_path_point DB_PATH_POINT
                        Path to the databases for PointFinder
  -g SPECIFIC_GENE [SPECIFIC_GENE ...]
                        Specify genes existing in the database to search for -
                        if none is specified all genes are included in the
                        search.
  -u, --unknown_mut     Show all mutations found even if in unknown to the
                        resistance database
```

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

Not published yet.


License
=======

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
