# DESCRIPTION

ResFinder identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria. Program uses the newest version of **blastn**.

# Installation

```bash
# Clone repository with dependencies and install pip3 modules
$ git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git --recursive
$ pip3 install -r requirements.txt
```

**Install standalone dependencies**

* [BioPython](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

The BlastAll and FormatDB that `resfinder.pl` uses are no longer available for downloading through ncbi. Therefor we have provided the resfinder.py script that uses Blastn instead. Note, this is no script that is running on the CGE server. The CGE server is running perl script using BlastAll.

**Warning:** Due to bugs in the BioPython 1.74, do not use this version if not using Python 3.7.


# Usage
```
# ./resfinder.py -h
$ python3 resfinder.py -h
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
  -q, --quiet
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
