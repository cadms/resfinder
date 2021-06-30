# GeneResult class test

## initialize

First part just creates some dummy objects needed for testing the class. A Result object containing a database is created. And a dictionary containing
results from a ResFinder hit.

```python

>>> from cgelib.output.result import Result
>>> res = Result.init_software_result(name="ResFinder", gitdir=".")
>>> res.init_database("ResFinder", ".")

>>> rf_dat_blast = {}
>>> rf_dat_blast["sbjct_header"] = "blaOXA-384_1_KF986263"
>>> rf_dat_blast["perc_ident"] = 97
>>> rf_dat_blast["HSP_length"] = 100
>>> rf_dat_blast["sbjct_length"] = 90
>>> rf_dat_blast["sbjct_start"] = 1
>>> rf_dat_blast["sbjct_end"] = 90
>>> rf_dat_blast["contig_name"] = "Contig01"
>>> rf_dat_blast["query_start"] = 701
>>> rf_dat_blast["query_end"] = 801
>>> rf_dat_blast["coverage"] = 1

>>> rf_dat_kma = {}
>>> rf_dat_kma["sbjct_header"] = "blaOXA-384_1_KF986263"
>>> rf_dat_kma["perc_ident"] = 97
>>> rf_dat_kma["HSP_length"] = 100
>>> rf_dat_kma["sbjct_length"] = 90
>>> rf_dat_kma["sbjct_start"] = 1
>>> rf_dat_kma["sbjct_end"] = 90
>>> rf_dat_kma["contig_name"] = "NA"
>>> rf_dat_kma["query_start"] = "NA"
>>> rf_dat_kma["query_end"] = "NA"
>>> rf_dat_kma["perc_coverage"] = 100
>>> rf_dat_kma["depth"] = 21

```

## Initialize GeneResult

```python

>>> from cge.output.gene_result import GeneResult

>>> gene_result_blast = GeneResult(res, rf_dat_blast, "ResFinder")
>>> assert(gene_result_blast["type"] == "seq_region")
>>> assert(gene_result_blast["gene"] is True)
>>> assert(gene_result_blast["ref_id"] == rf_dat_blast["sbjct_header"])
>>> assert(gene_result_blast["name"] == "blaOXA-384")
>>> assert(gene_result_blast["ref_acc"] == "KF986263")
>>> assert(gene_result_blast["identity"] == rf_dat_blast["perc_ident"])
>>> assert(gene_result_blast["alignment_length"] == rf_dat_blast["HSP_length"])
>>> assert(gene_result_blast["ref_seq_lenght"] == rf_dat_blast["sbjct_length"])
>>> assert(gene_result_blast["ref_start_pos"] == rf_dat_blast["sbjct_start"])
>>> assert(gene_result_blast["ref_end_pos"] == rf_dat_blast["sbjct_end"])
>>> assert(gene_result_blast["query_id"] == rf_dat_blast["contig_name"])
>>> assert(gene_result_blast["query_start_pos"] == rf_dat_blast["query_start"])
>>> assert(gene_result_blast["query_end_pos"] == rf_dat_blast["query_end"])
>>> assert(gene_result_blast["key"] == "blaOXA-384;;1;;KF986263")
>>> gene_result_blast["ref_database"]
... #doctest: +ELLIPSIS
'ResFinder-...'
>>> assert(gene_result_blast["coverage"] == 100)

>>> gene_result_kma = GeneResult(res, rf_dat_kma, "ResFinder")

```

## Methods

### remove_NAs()

```python

res["test_NA"] = "NA"
res["test_None"] = None
res.remove_NAs()
assert("test_NA" not in res)
assert("test_None" not in res)

```

### random_string(str_len)

```python

>>> rnd_str = GeneResult.random_string(str_len=7)
>>> import re
>>> assert(re.search("[a-z]{7}", rnd_str))

```

### get_rnd_unique_gene_key(gene_key, res_collection, minimum_gene_key, delimiter)

```python

GeneResult.get_rnd_unique_gene_key()

```

## Private methods

### _split_sbjct_header(header)

```python

>>> t, v, a = GeneResult._split_sbjct_header("blaOXA-384_1_KF986263")
>>> assert(t == "blaOXA-384")
>>> assert(v == "1")
>>> assert(a == "KF986263")
>>> t, v, a = GeneResult._split_sbjct_header("blaOXA-384")
>>> assert(t == "blaOXA-384")
>>> assert(v is None)
>>> assert(a is None)
>>> t, v, a = GeneResult._split_sbjct_header("blaCMY-150_2_NG_060513")
>>> assert(t == "blaCMY-150")
>>> assert(v == "2")
>>> assert(a == "NG_060513")

```

### _get_unique_gene_key(res_collection, delimiter)

```python

>>> import copy
>>> res_unique_gene_test_blast = copy.deepcopy(res)

>>> gene_result_blast._get_unique_gene_key(res_unique_gene_test_blast)
'blaOXA-384;;1;;KF986263'
>>> res_unique_gene_test_blast.add_class(cl="seq_regions", **gene_result_blast)
>>> gene_result_blast["query_id"] = "Contig02"
>>> gene_result_blast._get_unique_gene_key(res_unique_gene_test_blast)
... #doctest: +ELLIPSIS
'blaOXA-384;;1;;KF986263;;...'

```

```python

>>> res_unique_gene_test_kma = copy.deepcopy(res)

>>> gene_result_kma._get_unique_gene_key(res_unique_gene_test_kma)
'blaOXA-384;;1;;KF986263'
>>> res_unique_gene_test_kma.add_class(cl="seq_regions", **gene_result_kma)
>>> gene_result_kma._get_unique_gene_key(res_unique_gene_test_kma)
... #doctest: +ELLIPSIS
'blaOXA-384;;1;;KF986263;;...'

```
