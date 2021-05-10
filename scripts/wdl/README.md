# Quick guide to running ResFinder with Cromwell

### Disclaimer
Support is not offered for running Cromwell and no files in this directory is
guaranteed to work. These files were uploaded as inspiration. Please do not
report issues relating to this directory.

## Prepare input files

Two input files are needed:
1. input_data.tsv
2. input.json

### input_data.tsv
Tab separated file. Should contain three columns in the following order:
1. Absolute path to fastq file 1
2. Absolute path to fastq file 2
3. Species

Each row should contain a single sample. If species cannot be provided put
"other" (cases sensitive).

#### Example
```

 /test/data/test_isolate_01_1.fq	 /test/data/test_isolate_01_2.fq	Escherichia coli
 /test/data/test_isolate_05_1.fq	 /test/data/test_isolate_05_2.fq	Escherichia coli
 /test/data/test_isolate_09a_1.fq	 /test/data/test_isolate_09a_2.fq	Escherichia coli
 /test/data/test_isolate_09b_1.fq	 /test/data/test_isolate_09b_2.fq	Escherichia coli

```

### input.json
JSON formatted file containing input and output information.

The file should consist of a single dict/hash/map with two keys:
* Resistance.inputSamplesFile
* Resistance.outputDir

The values should be the absolute path to the input_data.tsv and the desired
output directory, respectively.

#### Example

```json

{
    "Resistance.inputSamplesFile": "../../tests/data/wdl_input.tsv",
    "Resistance.outputDir": "/home/rkmo/delme/"
}

```

## Set environment variables

```bash

export CGE_PYTHON="python3"
export CGE_KMA="kma"
export CGE_BLASTN="blastn"
export CGE_RESFINDER="/tools/resfinder/run_resfinder.py"
export CGE_RESFINDER_RESGENE_DB="/tools/resfinder/db_resfinder"
export CGE_RESFINDER_RESPOINT_DB="/tools/resfinder/db_pointfinder"
export CGE_RESFINDER_DISINF_DB="/tools/resfinder/db_disinfinder"
export CGE_RESFINDER_GENE_COV="0.6"
export CGE_RESFINDER_GENE_ID="0.8"
export CGE_RESFINDER_POINT_COV="0.6"
export CGE_RESFINDER_POINT_ID="0.8"

```

## Run Cromwell

```bash

java -jar /home/mydir/wdls/resfinder.wdl --inputs input.json

```
