===================
ResFinder
===================

This project documents ResFinder service


Documentation
=============

## What is it?

The ResFinder service contains one perl script *ResFinder-2.1.pl* which is the script of the lates
version of the ResFinder service. ResFinder identifies acquired antimicrobial resistance genes in total
or partial sequenced isolates of bacteria.

## Usage

To use the service some data needs to be pre-installed: *database*, *blast-2.2.26* and *Makefile*

The folder *database* includes the ResFinder datasets and needs to be manually updataed frequently to get the best
results. The datasets can be downloaded from http://cge.cbs.dtu.dk/services/data.php.

The folder *blast-2.2.26* includes blastall and formatdb which are used by the script

The file *Makefile* installs the nesserary perl modules to run the script. It is used by writing:
    make install

The program can be invoked with the -h option to get help and more information of the service.

```bash
Usage: perl ResFinder-2.1.pl [options]

Options:

    -h HELP
                    Prints a message with options and information to the screen
    -d DATABASE
                    The path to where you have located the database folder
    -b BLAST
                    The path to the location of blast-2.2.26
    -i INFILE
                    Your input file which needs to be preassembled partial or complete genomes in fasta format
    -o OUTFOLDER
                    The folder you want to have your output files places
    -a ANTIMICROBIAL
                    Antimicrobial configuration. The options can be found in the file *ResFinder_Antimicrobial*
    -k  THRESHOLD
                    The threshold for % identity for example '95.00' for 95 %
    -l  MIN_LENGHT
                    The minimum length of the overlap ex 0.60 for an overlap of minimum 60 %
```

## Example of use with the *database* and *blast-2.2.26* folder loacted in the current directory
    
    perl ResFinder-2.1.pl -i test.fsa -o OUTFOLDER -a aminoglycoside -k 95.00 -l 0.60

## Example of use with the *database* and *blast-2.2.26* folder loacted in antoher directory

    perl ResFinder-2.1.pl -d path/to/database -blast path/to/blast-2.2.26 -i test.fsa -o OUTFOLDER -a aminoglycoside -k 95.00 -l 0.60
    
## Web-server

A webserver implementing the methods is available at the [CGE website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/ResFinder/


## The Latest Version


The latest version can be found at
https://bitbucket.org/genomicepidemiology/resfinder/overview

## Documentation


The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/resfinder/overview.

Installation
=======

The scripts are self contained. You just have to copy them to where they should
be used. Only the *database* folder needs to be updataed mannually. 

Remember to add the program to your system path if you want to be able to invoke the program without calling the full path.
If you don't do that you have to write the full path to the program when using it.

Citation
=======

When using the method please cite:

Identification of acquired antimicrobial resistance genes.
Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV.
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
