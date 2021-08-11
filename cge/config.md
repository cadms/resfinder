# Test Config class

**Important**: In order for the tests to work, you either need to have the ResFinder and PointFinder database installed at the default locations "db_resfinder" and "db_pointfinder" inside the ResFinder application folder or you need to have set the environment variables "CGE_RESFINDER_RESGENE_DB" and "CGE_RESFINDER_RESPOINT_DB" to valid databases.

## Setup

```python

>>> from cge.config import Config

>>> class DummyArgs():
...     def __init__(self):
...         self.inputfasta = None
...         self.inputfastq = None
...         self.outputPath = "./"
...         self.blastPath = None
...         self.kmaPath = None
...         self.species = None
...         self.ignore_missing_species = None
...         self.db_path_res = None
...         self.db_path_res_kma = None
...         self.databases = None
...         self.acquired = True
...         self.acq_overlap = None
...         self.min_cov = None
...         self.threshold = None
...         self.point = True
...         self.db_path_point = None
...         self.db_path_point_kma = None
...         self.specific_gene = None
...         self.unknown_mut = None
...         self.min_cov_point = None
...         self.threshold_point = None
...         self.pickle = False

>>> args = DummyArgs()
>>> args1 = DummyArgs()
>>> args2 = DummyArgs()
>>> args3 = DummyArgs()
>>> args4 = DummyArgs()
>>> args5 = DummyArgs()
>>> args6 = DummyArgs()
>>> args7 = DummyArgs()
>>> args8 = DummyArgs()
>>> args9 = DummyArgs()

```

## Config(args)

### No input

```python

>>> args.species = "Other"
>>> args.point = False
>>> conf = Config(args)
>>> print(conf.species)
None
>>> conf.amr_abbreviations["Amoxicillin"][0]
'AMO'

```

### FASTA

```python

>>> args3 = DummyArgs()
>>> fasta_filepath = "{}/{}".format(conf.resfinder_root,
...     "tests/data/test_isolate_01.fa")
>>> args3.inputfasta = fasta_filepath
>>> args3.species = "ecoli"
>>> conf_fasta = Config(args3)
>>> conf_fasta.inputfasta
... #doctest: +ELLIPSIS
'/...test_isolate_01.fa'
>>> conf_fasta.sample_name
'test_isolate_01.fa'

```

### FASTQ

```python

>>> fastq_filepath_1 = "{}/{}".format(conf.resfinder_root,
...     "tests/data/test_isolate_01_1.fq")
>>> fastq_filepath_2 = "{}/{}".format(conf.resfinder_root,
...     "tests/data/test_isolate_01_2.fq")
>>> args7.inputfastq = [fastq_filepath_1, fastq_filepath_2]
>>> args7.species = "ecoli"
>>> conf_fastq = Config(args7)
>>> (conf_fastq.inputfastq_1, conf_fastq.inputfastq_2)
... #doctest: +ELLIPSIS
('...test_isolate_01_1.fq', '...test_isolate_01_2.fq')
>>> conf_fastq.sample_name
'test_isolate_01_1.fq'
>>> conf_fastq.kma
... #doctest: +ELLIPSIS
'...kma'

```

## get_abs_path_and_check(path)

```python

>>> Config.get_abs_path_and_check(fasta_filepath)
... #doctest: +ELLIPSIS
'/...tests/data/test_isolate_01.fa'

>>> Config.get_abs_path_and_check("/file/not/found")
Traceback (most recent call last):
SystemExit: ERROR: Path not found: /file/not/found

```

## get_prg_path(args)

This test will fail if 'blastn' isn't in PATH.

```python

>>> blastPath = "blastn"
>>> Config.get_prg_path(blastPath)
'blastn'

```

## get_species(in_species, species_def_filepath)

```python

>>> species_def_filepath = "{}/{}".format(conf.resfinder_root, "README.md")
>>> Config.get_species("Escherichia coli", species_def_filepath)
'escherichia coli'
>>> Config.get_species("ecoli", species_def_filepath)
'escherichia coli'

>>> str(Config.get_species("Other", species_def_filepath))
'None'


```

## set_path_pointdb(args)

Running of test is done via constructor (Config(args)) in 'Fasta' paragraph above. Result of run is asserted here (below).

```python

>>> assert(conf_fasta.db_path_point_root is not None)
>>> assert(conf_fasta.db_path_point is not None)

>>> args9.species = "Too many words"
>>> Config(args9)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: ERROR: Species name must contain 1 or 2 names...


>>> args9 = DummyArgs()
>>> args9.species = "campylobacter"
>>> conf9 = Config(args9)
>>> conf9.db_path_point[-13:]
'campylobacter'

>>> args9 = DummyArgs()
>>> args9.ignore_missing_species = True
>>> args9.species = "missing species"
>>> conf9 = Config(args9)
>>> str(conf9.db_path_point)
'None'
>>> str(conf9.db_path_point_root)
'None'

```

## set_path_resfinderdb(args)

Test is also done through the FASTA and FASTQ sections above.

```python


```

## _parse_species_dir(path_pointdb, species_dir, ignore_missing_species)

Successful test depends on either environment variabel 'CGE_RESFINDER_RESPOINT_DB' is set or database is located at default location.

```python

>>> conf_fasta._parse_species_dir(args7.db_path_point, "escherichia_coli", False)
'escherichia_coli'

>>> args8.species = "campylobacter coli"
>>> conf8 = Config(args8)
>>> conf8._parse_species_dir(args8.db_path_point, "campylobacter_coli", False)
'campylobacter'

>>> conf8.species = "fail_species"
>>> conf8._parse_species_dir(args8.db_path_point, "fail_species", False)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: ERROR: species...

>>> str(conf8._parse_species_dir(args8.db_path_point, "fail_species", True))
'None'

```

## _set_default_values(args)

```python

>>> Config._set_default_values(args1)
>>> args1.kmaPath
'kma'

>>> args1.kmaPath = "some/path/to/kma"
>>> Config._set_default_values(args1)
>>> args1.kmaPath
'some/path/to/kma'

```

## set_default_and_env_vals(args, env_def_filepath)

```python

>>> env_def_filepath = "{}/{}".format(conf.resfinder_root, "README.md")

>>> print(args2.min_cov)
None
>>> import os
>>> os.environ["CGE_RESFINDER_GENE_COV"] = "0.1"
>>> args2.threshold = 0.9

>>> Config.set_default_and_env_vals(args2, env_def_filepath)
>>> args2.min_cov
'0.1'
>>> args2.threshold
0.9
>>> args2.unknown_mut
False

>>> env_def_fail_filepath = "{}/{}".format(conf.resfinder_root,
...     "tests/data/env_var_test_fail.md")
>>> Config.set_default_and_env_vals(args2, env_def_fail_filepath)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: ERROR: A flag set...kmaPathWrong.

```

## set_resfinder_opts(args)

Test errors. Success is tested via constructor (Config(args)) above.

```python

>>> args5.db_path_res = os.path.expanduser("~")
>>> Config(args5)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: Input Error: The database config or notes.txt file could not be...

>>> args5.db_path_res = None
>>> args5.min_cov = 1.2
>>> Config(args5)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: ERROR: Minimum coverage above 1 or below 0...Given value: 1.2.

>>> args5.min_cov = 0.65
>>> args5.threshold = 80
>>> Config(args5)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: ERROR: Threshold for identity of ResFinder...Given value: 80.0.

```

## set_pointfinder_opts(args)

```python

>>> args6.species = "Escherichia coli"
>>> args6.point = True
>>> conf_args6 = Config(args6)
>>> conf_args6.species
'escherichia coli'
>>> conf_args6.point
True

>>> args6.species = None
>>> Config(args6)
... #doctest: +ELLIPSIS
Traceback (most recent call last):
SystemExit: ERROR: Chromosomal point mutations cannot be located if no specie...

>>> args6.point
True
>>> args6.ignore_missing_species = True
>>> conf_args6_no_species = Config(args6)
>>> conf_args6_no_species.point
False
>>> str(conf_args6_no_species.species)
'None'

```
