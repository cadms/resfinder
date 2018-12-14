ResFinder documentation
=============

ResFinder identifies acquired antimicrobial resistance genes in total or partial
sequenced isolates of bacteria.

## Content of the repository
1. run_resfinder.py - Use this script to run ResFinder
2. tests/data       - Contains fasta and fastq data for testing. More information in the "Test data" section
3. scripts/         - All scripts in this directory is unsupported but has been uploaded as they may be useful
4. cge/             - ResFinder code

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
9. Test installation

A small script has been written to automate this process. It is available from the
scripts directory and is named install_resfinder.sh. It is very simple and might
not work in all environments. It is only meant as a supplement and no support will 
be provided for any scripts in this directory. However, specific suggestions (with code) 
for improvement is very welcome.

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
specified every time you run ResFinder unless you add the kma program to your PATH.
```bash
# Go to the directoy in which you installed the ResFinder tool
cd /path/to/some/dir/resfinder
cd cge
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
```

### Databases
This section describes how to install the databases at the ResFinder default locations. 
The database locations can be changed, but must then be specified to ResFinder at run time.

#### ResFinder database
```
# Go to the directoy in which you installed the ResFinder tool
cd /path/to/some/dir/resfinder
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
```

#### PointFinder database
```
# Go to the directoy in which you installed the ResFinder tool
cd /path/to/some/dir/resfinder
git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder
```

#### Indexing databases with KMA
```
# Go to the directoy in which you installed the ResFinder tool
cd /path/to/some/dir/resfinder
mkdir db_resfinder/kma_indexing
cge/kma/kma_index -i db_resfinder/fusidicacid.fsa -o db_resfinder/kma_indexing/fusidicacid
cge/kma/kma_index -i db_resfinder/phenicol.fsa -o db_resfinder/kma_indexing/phenicol
cge/kma/kma_index -i db_resfinder/glycopeptide.fsa -o db_resfinder/kma_indexing/glycopeptide
cge/kma/kma_index -i db_resfinder/trimethoprim.fsa -o db_resfinder/kma_indexing/trimethoprim
cge/kma/kma_index -i db_resfinder/oxazolidinone.fsa -o db_resfinder/kma_indexing/oxazolidinone
cge/kma/kma_index -i db_resfinder/tetracycline.fsa -o db_resfinder/kma_indexing/tetracycline
cge/kma/kma_index -i db_resfinder/quinolone.fsa -o db_resfinder/kma_indexing/quinolone
cge/kma/kma_index -i db_resfinder/nitroimidazole.fsa -o db_resfinder/kma_indexing/nitroimidazole
cge/kma/kma_index -i db_resfinder/fosfomycin.fsa -o db_resfinder/kma_indexing/fosfomycin
cge/kma/kma_index -i db_resfinder/aminoglycoside.fsa -o db_resfinder/kma_indexing/aminoglycoside
cge/kma/kma_index -i db_resfinder/macrolide.fsa -o db_resfinder/kma_indexing/macrolide
cge/kma/kma_index -i db_resfinder/sulphonamide.fsa -o db_resfinder/kma_indexing/sulphonamide
cge/kma/kma_index -i db_resfinder/rifampicin.fsa -o db_resfinder/kma_indexing/rifampicin
cge/kma/kma_index -i db_resfinder/colistin.fsa -o db_resfinder/kma_indexing/colistin
cge/kma/kma_index -i db_resfinder/beta-lactam.fsa -o db_resfinder/kma_indexing/beta-lactam
cge/kma/kma_index -i db_pointfinder/campylobacter/*.fsa -o db_pointfinder/campylobacter/campylobacter
cge/kma/kma_index -i db_pointfinder/escherichia_coli/*.fsa -o db_pointfinder/escherichia_coli/escherichia_coli
cge/kma/kma_index -i db_pointfinder/enterococcus_faecalis/*.fsa -o db_pointfinder/enterococcus_faecalis/enterococcus_faecalis
cge/kma/kma_index -i db_pointfinder/enterococcus_faecium/*.fsa -o db_pointfinder/enterococcus_faecium/enterococcus_faecium
cge/kma/kma_index -i db_pointfinder/neisseria_gonorrhoeae/*.fsa -o db_pointfinder/neisseria_gonorrhoeae/neisseria_gonorrhoeae
```

### Test ResFinder intallation
If you did not install BLAST, test 1 and 3 will fail. If you did not install KMA, test 2
and 4 will fail.
The 4 tests will in total take approximately take 5-60 seconds, depending on your system.
```
# Go to the directoy in which you installed the ResFinder tool
cd /path/to/some/dir/resfinder

# Run tests
python3 tests/functional_test.py

# Output from successful tests
....
----------------------------------------------------------------------
Ran 4 tests in 8.263s

OK
```

### Test data
Test data can be found in the sub-dierectory /tests/data

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
