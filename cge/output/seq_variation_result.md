# SeqVariationResult class test

## initialize

First part just creates some dummy objects needed for testing the class. A
Result object containing the databases are created, a list of GeneResult objects
(actually only one object), and a list containing results from a PointFinder
hit.

```python

>>> from cgelib.output.result import Result
>>> from cge.output.gene_result import GeneResult

>>> res = Result.init_software_result(name="ResFinder", gitdir=".")
>>> res.init_database("ResFinder", ".")
>>> res.init_database("PointFinder", ".")

>>> pf_dat_gyra_kma = {}
>>> pf_dat_gyra_kma["sbjct_header"] = "gyrA"
>>> pf_dat_gyra_kma["perc_ident"] = 99.92
>>> pf_dat_gyra_kma["HSP_length"] = 2628
>>> pf_dat_gyra_kma["sbjct_length"] = 2628
>>> pf_dat_gyra_kma["sbjct_start"] = 1
>>> pf_dat_gyra_kma["sbjct_end"] = 2628
>>> pf_dat_gyra_kma["contig_name"] = "gyrA_G81D_GAT_D82G_GGC"
>>> pf_dat_gyra_kma["query_start"] = 1
>>> pf_dat_gyra_kma["query_end"] = 2628
>>> pf_dat_gyra_kma["perc_coverage"] = 100.0
>>> pf_dat_gyra_kma["depth"] = 21

>>> gene_result_kma = GeneResult(res, pf_dat_gyra_kma, "PointFinder")
>>> genes = [ gene_result_kma ]

>>> pf_dat = []
>>> pf_dat.append("sub")
>>> pf_dat.append(81)
>>> pf_dat.append(81)
>>> pf_dat.append("D")  # not used
>>> pf_dat.append("p.G81D")
>>> pf_dat.append("GGT")
>>> pf_dat.append("GAT")
>>> pf_dat.append("G")
>>> pf_dat.append("D")

```


## Initialize SeqVariationResult

### Test with PointFinder hit list > 7

```python

>>> from cge.output.seq_variation_result import SeqVariationResult

>>> seqvar_result = SeqVariationResult(res, pf_dat, genes, "PointFinder")
>>> print(seqvar_result)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
{'type': 'seq_variation',
  'seq_var': 'p.G81D',
  'ref_codon': 'ggt',
  'var_codon': 'gat',
  'codon_change': 'ggt>gat',
  'ref_aa': 'g',
  'var_aa': 'd',
  'ref_start_pos': 81,
  'ref_end_pos': 81,
  'substitution': True,
  'deletion': False,
  'insertion': False,
  'ref_id': 'gyrA;;81;;d',
  'key': 'gyrA;;81;;d',
  'ref_database': 'PointFinder-...',
  'seq_regions': ['gyrA']}

```

### Test with PointFinder hit list <= 7

```python

>>> pf_dat2 = []
>>> pf_dat2.append("sub")
>>> pf_dat2.append(81)
>>> pf_dat2.append(81)
>>> pf_dat2.append("D")  # not used
>>> pf_dat2.append("p.G81D")
>>> pf_dat2.append("GGT")
>>> pf_dat2.append("GAT")

>>> seqvar_result2 = SeqVariationResult(res, pf_dat2, genes, "PointFinder")
>>> print(seqvar_result2)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
{'type': 'seq_variation',
  'seq_var': 'p.G81D',
  'ref_codon': 'ggt',
  'var_codon': 'gat',
  'codon_change': 'ggt>gat',
  'ref_start_pos': 81,
  'ref_end_pos': 81,
  'substitution': True,
  'deletion': False,
  'insertion': False,
  'ref_id': 'gyrA;;81;;gat',
  'key': 'gyrA;;81;;gat',
  'ref_database': 'PointFinder-...',
  'seq_regions': ['gyrA']}

```

## Methods

### load_var_type(type)

```python

>>> seqvar_result.load_var_type("ins")
>>> assert(seqvar_result["substitution"] is False)
>>> assert(seqvar_result["insertion"] is True)
>>> assert(seqvar_result["deletion"] is False)
>>> seqvar_result.load_var_type("del")
>>> assert(seqvar_result["substitution"] is False)
>>> assert(seqvar_result["insertion"] is False)
>>> assert(seqvar_result["deletion"] is True)
>>> seqvar_result.load_var_type("sub")
>>> assert(seqvar_result["substitution"] is True)
>>> assert(seqvar_result["insertion"] is False)
>>> assert(seqvar_result["deletion"] is False)

```

## Private methods

### _get_unique_seqvar_key(res_collection, minimum_key, delimiter)

```python

>>> res.add_class(cl="seq_variations", **seqvar_result)
>>> for k in res["seq_variations"].keys():
...     print(k)
gyrA;;81;;d

>>> SeqVariationResult._get_unique_seqvar_key(res, "minkey", "||")
'minkey'
>>> SeqVariationResult._get_unique_seqvar_key(res, "gyrA;;81;;d")
... #doctest: +ELLIPSIS
'gyrA;;81;;d;;...'

```
