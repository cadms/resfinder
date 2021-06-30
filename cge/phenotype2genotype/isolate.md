# Test Isolate class

**Missing tests***

## Setup

### phenoDB object

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

### Result object

More details about the setup of the result object can be found in the
std_results test documentation.

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

>>> from cge2.output.result import Result
>>> from cge.output.gene_result import GeneResult

>>> res = Result.init_software_result(name="ResFinder", gitdir=".")
>>> res.init_database("ResFinder", ".")
>>> res.init_database("PointFinder", ".")

>>> from cge.output.std_results import ResFinderResultHandler
>>> ResFinderResultHandler.standardize_results(res,
...                                            rf_custom_kma,
...                                            "ResFinder")

>>> from cge.output.std_results import PointFinderResultHandler
>>> PointFinderResultHandler.standardize_results(res,
...                                              pf_custom_kma,
...                                              "PointFinder")

```

## Init Isolate

```python

>>> from cge.phenotype2genotype.isolate import Isolate
>>> isolate = Isolate(name="Test sample")

```

## Methods

### load_finder_results(std_table, phenodb, type)

```python

>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_regions")
>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_variations")

>>> isolate["blaOXA-384_KF986263"][0].ab_class.pop()
'beta-lactam'
>>> isolate["gyrA_81_d"][0].ref_codon
'ggt'
>>> isolate["gyrA_81_d"][0].mut_codon
'gat'
>>> isolate["gyrA_81_d"][0].ab_class.pop()
'fluoroquinolone'
>>> isolate["gyrA_82_g"][0].ref_codon
'gac'
>>> isolate["gyrA_82_g"][0].mut_codon
'ggc'
>>> isolate["gyrA_82_g"][0].ab_class.pop()
'fluoroquinolone'

```

### calc_res_profile(res_pheno_db)

```python

>>> isolate.calc_res_profile(res_pheno_db)

```

### get_phenodb_id(feat_res_dict, type)

```python

>>> feat_res_dict = res["seq_regions"]["blaOXA-384;;1;;KF986263"]
>>> Isolate.get_phenodb_id(feat_res_dict, "seq_regions")
'blaOXA-384_KF986263'

>>> feat_res_dict = res["seq_variations"]["gyrA;;81;;d"]
>>> Isolate.get_phenodb_id(feat_res_dict, "seq_variations")
'gyrA_81_d'

```
