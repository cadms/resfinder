# std_results tests

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

```

ResFinder stores results in a dict of dicts. The internal dicts represents hits
to the sub-database the dict represent. The sub-database is for legacy reasons
named after amr classes and to some extend represents hits to genes in that amr
class. This is not a naming convention that will be guaranteed in future
releases as the structure is not optimal.

Here a dummy hit (obtained using KMA) is created, named rf\_dat\_kma. It is stored in the beta-lactam sub-database/dict in the ResFinder result dict named
rf\_custom\_kma.

```python

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

>>> rf_custom_kma = {}
>>> rf_custom_kma["excluded"] = {}
>>> rf_custom_kma["aminoglycoside"] = "No hit found"
>>> rf_custom_kma["beta-lactam"] = {"unique_hit_key": rf_dat_kma}

```

The PointFinder part of ResFinder unfortunately stores its results differently than ResFinder. The results are stored in a dict (here stored in the variable 'pf\_custom\_blast'), which contains a dict entry for each gene searched for mutations in the given species (corresponding dummy below is named 'gyrA'), but only if the gene is found. If a gene is not found its not a dict but the string "No hit found". A gene dict then has an entry for each hit to the reference gene in the database. The hit entry itself is also a dict (dummy hit below is named gyrA_hit).

The actual point mutations are stored in the gene dict (dummy: gyrA) in a list
of lists, where each list in the list describes a mutation.

```python

>>> gyrA_hit = {}
>>> gyrA_hit["evalue"] = 0.0
>>> gyrA_hit["sbjct_header"] = "gyrA"
>>> gyrA_hit["bit"] = 4843.04
>>> gyrA_hit["perc_ident"] = 99.92
>>> gyrA_hit["sbjct_length"] = 2628
>>> gyrA_hit["sbjct_start"] = 1
>>> gyrA_hit["sbjct_end"] = 2628
>>> gyrA_hit["gaps"] = 0
>>> gyrA_hit["contig_name"] = "gyrA_G81D_GAT_D82G_GGC"
>>> gyrA_hit["query_start"] = 1
>>> gyrA_hit["query_end"] = 2628
>>> gyrA_hit["HSP_length"] = 2628
>>> gyrA_hit["coverage"] = 1.0
>>> gyrA_hit["cal_score"] = 99.92
>>> gyrA_hit["hit_id"] = "gyrA_G81D_GAT_D82G_GGC:1..2628:gyrA:99.923896"
>>> gyrA_hit["strand"] = 0
>>> gyrA_hit["perc_coverage"] = 100.0

>>> gyrA = {}
>>> gyrA["mis_matches"] = [
...   [ 'sub', 81, 81, 'D', 'p.G81D', 'GGT', 'GAT', 'G', 'D' ],
...   [ 'sub', 82, 82, 'G', 'p.D82G', 'GAC', 'GGC', 'D', 'G' ] ]
>>> gyrA["hits"] = {}
>>> gyrA["hits"]["gyrA_G81D_GAT_D82G_GGC:1..2628:gyrA:99.923896"] = gyrA_hit

>>> pf_custom_blast = {}
>>> pf_custom_blast["excluded"] = {}
>>> pf_custom_blast["gyrA"] = gyrA
>>> pf_custom_blast["gyrB"] = "No hit found"

```

Create dummy data for KMA results, which are handled slightly differently for
PointFinder.

```python

>>> gyrA_kma_hit = {}
>>> gyrA_kma_hit["sbjct_header"] = "gyrA"
>>> gyrA_kma_hit["perc_ident"] = 99.92
>>> gyrA_kma_hit["HSP_length"] = 2628
>>> gyrA_kma_hit["sbjct_length"] = 2628
>>> gyrA_kma_hit["sbjct_start"] = 1
>>> gyrA_kma_hit["sbjct_end"] = 2628
>>> gyrA_kma_hit["contig_name"] = "NA"
>>> gyrA_kma_hit["query_start"] = "NA"
>>> gyrA_kma_hit["query_end"] = "NA"
>>> gyrA_kma_hit["perc_coverage"] = 100.0
>>> gyrA_kma_hit["depth"] = 21
>>> gyrA_kma_hit["mis_matches"] = [
...   [ 'sub', 81, 81, 'D', 'p.G81D', 'GGT', 'GAT', 'G', 'D' ],
...   [ 'sub', 82, 82, 'G', 'p.D82G', 'GAC', 'GGC', 'D', 'G' ] ]

>>> pf_custom_kma = {}
>>> pf_custom_kma["excluded"] = {}
>>> pf_custom_kma["gyrA"] = gyrA_kma_hit
>>> pf_custom_kma["gyrB"] = "No hit found"

>>> import copy
>>> res_kma_test = copy.deepcopy(res)

```

Create the phenoDB object.

```python

>>> import os
>>> from cge.phenotype2genotype.res_profile import PhenoDB

>>> resfinder_db_path = os.environ["CGE_RESFINDER_RESGENE_DB"]
>>> assert(len(resfinder_db_path) > 0)
>>> pointfinder_db_path = os.environ["CGE_RESFINDER_RESPOINT_DB"]
>>> assert(len(pointfinder_db_path) > 0)

>>> abclassdef_file= "{}/antibiotic_classes.txt".format(resfinder_db_path)
>>> acquired_file= "{}/phenotypes.txt".format(resfinder_db_path)
>>> point_file = ("{}/escherichia_coli/resistens-overview.txt"
...               .format(pointfinder_db_path))
>>> res_pheno_db = PhenoDB(abclassdef_file=abclassdef_file,
...                        acquired_file=acquired_file,
...                        point_file=point_file)

```

## Test ResFinderResultHandler

### standardize\_results(res\_collection, res, ref\_db\_name)

```python

>>> from cge.output.std_results import ResFinderResultHandler
>>> ResFinderResultHandler.standardize_results(res,
...                                            rf_custom_kma,
...                                            "ResFinder")

>>> for k in res["databases"]:
...   print(k)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
ResFinder-...
PointFinder-...

>>> for k in res["seq_regions"]:
...   print(k)
blaOXA-384;;1;;KF986263

```

### load\_res\_profile(res\_collection, isolate)

#### Setup

Create Isolate object

```python

>>> from cge.phenotype2genotype.isolate import Isolate
>>> isolate = Isolate(name="Test sample")

>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_regions")
>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_variations")
>>> isolate.calc_res_profile(res_pheno_db)

>>> import os.path
>>> import inspect
>>> from cgelib.utils.loaders_mixin import LoadersMixin
>>> std_results_file = inspect.getfile(ResFinderResultHandler)
>>> std_results_dir = os.path.dirname(os.path.realpath(std_results_file))
>>> amr_abbreviations_file = ("{}/../../amr_abbreviations.md"
...                           .format(std_results_dir))
>>> amr_abbreviations = LoadersMixin.load_md_table_after_keyword(
...     amr_abbreviations_file, "## Abbreviations")

```

#### Test

```python

>>> from cge.output.std_results import PointFinderResultHandler
>>> ResFinderResultHandler.load_res_profile(res, isolate, amr_abbreviations)
>>> res["phenotypes"]["ampicillin"]["seq_regions"]
[]
>>> res["result_summary"]
'UBE'

```

## Test PointFinderResultHandler

### standardize\_results(res\_collection, res, ref\_db\_name)

```python

>>> PointFinderResultHandler.standardize_results(res,
...                                              pf_custom_blast,
...                                              "PointFinder")

>>> for k in res["databases"]:
...   print(k)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
ResFinder-...
PointFinder-...

>>> for k in res["seq_regions"]:
...   print(k)
blaOXA-384;;1;;KF986263
gyrA

>>> for k in res["seq_variations"]:
...   print(k)
gyrA;;81;;d
gyrA;;82;;g

```

```python

>>> from cge.output.std_results import PointFinderResultHandler
>>> PointFinderResultHandler.standardize_results(res_kma_test,
...                                              pf_custom_kma,
...                                              "PointFinder")

>>> for k in res_kma_test["databases"]:
...   print(k)
... #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
ResFinder-...
PointFinder-...

>>> for k in res_kma_test["seq_regions"]:
...   print(k)
gyrA

>>> for k in res_kma_test["seq_variations"]:
...   print(k)
gyrA;;81;;d
gyrA;;82;;g

```

### create\_amr\_summary\_str(res_collection, amr_abbreviations)

#### Setup

```python

>>> phenotype1 = {
...     "type": "phenotype",
...     "key": "ampicillin+clavulanic acid",
...     "category": "amr",
...     "amr_classes": ['beta-lactam'],
...     "amr_resistance": "ampicillin+clavulanic acid",
...     "amr_resistant": True}
>>> phenotype2 = {
...     "type": "phenotype",
...     "key": "imipenem",
...     "category": "amr",
...     "amr_classes": ['beta-lactam'],
...     "amr_resistance": "imipenem",
...     "amr_resistant": True}
>>> phenotype3 = {
...     "type": "phenotype",
...     "key": "doxycycline",
...     "category": "amr",
...     "amr_classes": ['tetracycline'],
...     "amr_resistance": "doxycycline",
...     "amr_resistant": True}
>>> res.add_class(cl="phenotypes", clobber_warn=False, **phenotype1)
>>> res.add_class(cl="phenotypes", clobber_warn=False, **phenotype2)
>>> res.add_class(cl="phenotypes", clobber_warn=False, **phenotype3)

```

```python

>>> sum_str = ResFinderResultHandler.create_amr_summary_str(res, amr_abbreviations)
>>> sum_str
'AML_IMI_UBE_DOX'

```
