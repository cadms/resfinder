# PhenotypeResult class tests

## Setup

### Create feature (ResGene, ResMutation) and Antibiotics objects

```python

>>> from cge.phenotype2genotype.feature import ResGene, ResMutation
>>> rg = ResGene(unique_id="blaOXA-384_KF986263", start=1, end=90,
...              isolate="isolateA", ab_class=["beta-lactam"])
>>> rm = ResMutation(unique_id="gyrA;;81;;d", seq_region="gyrA", pos=81,
...                  ref_codon="ggt", mut_codon="gat", ref_aa="g",
...                  mut_aa="d", isolate="isolateB", nuc=False,
...                  ab_class=["fluoroquinolone"])

>>> from cge.phenotype2genotype.res_profile import Antibiotics
>>> ab_m1 = Antibiotics(name="ciprofloxacin", classes=["fluoroquinolone"],
...                     feature=rm)
>>> ab_m2 = Antibiotics(name="nalidixic acid", classes=["fluoroquinolone"],
...                     feature=rm)
>>> ab_g1 = Antibiotics(name="carbapenem", classes=["beta-lactam"],
...                     feature=rg)

```

### Create Isolate object

```python

>>> import os
>>> import copy
>>> from cgelib.output.result import Result
>>> from cge.phenotype2genotype.res_profile import PhenoDB
>>> from cge.phenotype2genotype.isolate import Isolate
>>> from cge.output.std_results import ResFinderResultHandler

>>> res = Result.init_software_result(name="ResFinder", gitdir=".")
>>> res.init_database("ResFinder", ".")
>>> res.init_database("PointFinder", ".")
>>> res_empty = copy.deepcopy(res)

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

>>> blaoxa = {
...  "type": "seq_region",
...  "gene": True,
...  "key": "blaOXA-384;;1;;KF986263",
...  "name": "blaOXA",
...  "identity": 98.2,
...  "alignment_length": 90,
...  "ref_seq_lenght": 100,
...  "coverage": 90.0,
...  "depth": 24.8,
...  "ref_id": "blaOXA-384_1_KF986263",
...  "ref_acc": "KF986263",
...  "ref_start_pos": 1,
...  "ref_end_pos": 100,
...  "ref_database": "ResFinder-unknown"
... }

>>> gyra = {
...   'type': 'seq_region',
...   'gene': True,
...   'ref_id': 'gyrA',
...   'name': 'gyrA',
...   'key': 'gyrA',
...   'identity': 99.92,
...   'alignment_length': 2628,
...   'ref_seq_lenght': 2628,
...   'depth': 21,
...   'ref_start_pos': 1,
...   'ref_end_pos': 2628,
...   'query_id': 'gyrA_G81D_GAT_D82G_GGC',
...   'query_start_pos': 1,
...   'query_end_pos': 2628,
...   'ref_database': 'PointFinder-a2b2ce4',
...   'coverage': 100.0
... }

>>> pg81d = {
...  'type': 'seq_variation',
...  'seq_var': 'p.G81D',
...  'ref_codon': 'ggt',
...  'var_codon': 'gat',
...  'codon_change': 'ggt>gat',
...  'ref_aa': 'g',
...  'var_aa': 'd',
...  'ref_start_pos': 81,
...  'ref_end_pos': 81,
...  'substitution': True,
...  'deletion': False,
...  'insertion': False,
...  'ref_id': 'gyrA;;81;;d',
...  'key': 'gyrA;;81;;d',
...  'ref_database': 'PointFinder-a2b2ce4',
...  'seq_regions': ['gyrA']
... }

>>> res.add_class(cl="seq_regions", **blaoxa)
>>> res.add_class(cl="seq_regions", **gyra)
>>> res.add_class(cl="seq_variations", **pg81d)

>>> isolate = Isolate(name="Test sample")
>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_regions")
>>> isolate.load_finder_results(std_table=res, phenodb=res_pheno_db,
...                             type="seq_variations")
>>> isolate.calc_res_profile(res_pheno_db)

```

## Initialise

```python

>>> from cge.output.phenotype_result import PhenotypeResult
>>> pheno = PhenotypeResult(ab_g1)
>>> pheno["key"]
'carbapenem'
>>> pheno["amr_classes"]
['beta-lactam']
>>> pheno["amr_resistance"]
'carbapenem'
>>> pheno["amr_resistant"]
False

```

## Methods

### set_resistant(res)

```python

>>> pheno.set_resistant(True)
>>> pheno["amr_resistant"]
True

```

### get_ref_id_and_type(feature, isolate)

```python

>>> PhenotypeResult.get_ref_id_and_type(rg, isolate)
('blaOXA-384_1_KF986263', 'seq_regions')
>>> PhenotypeResult.get_ref_id_and_type(rm, isolate)
('gyrA;;81;;d', 'seq_variations')

```

### get_keys_matching_ref_id(ref_id, res_collection)

```python

>>> refid, type = PhenotypeResult.get_ref_id_and_type(rg, isolate)
>>> PhenotypeResult.get_keys_matching_ref_id(refid, res[type])
['blaOXA-384;;1;;KF986263']
>>> PhenotypeResult.get_keys_matching_ref_id(refid, res_empty[type])
[]

```

### add_feature(res_collection, isolate, feature)

```python

>>> pheno.add_feature(res, isolate, rg)
>>> res.add_class(cl="phenotypes", **pheno)
>>> res["phenotypes"]["carbapenem"]["seq_regions"]
['blaOXA-384;;1;;KF986263']
>>> res["phenotypes"]["carbapenem"]["amr_classes"]
['beta-lactam']
>>> res["seq_regions"]["blaOXA-384;;1;;KF986263"]["phenotypes"]
['carbapenem']

>>> pheno2 = PhenotypeResult(ab_m1)
>>> pheno2.set_resistant(True)
>>> pheno2.add_feature(res, isolate, rm)
>>> res.add_class(cl="phenotypes", **pheno2)
>>> res["phenotypes"]["ciprofloxacin"]["seq_variations"]
['gyrA;;81;;d']
>>> res["phenotypes"]["ciprofloxacin"]["seq_regions"]
[]
>>> res["phenotypes"]["ciprofloxacin"]["amr_classes"]
['fluoroquinolone']
>>> res["seq_variations"]["gyrA;;81;;d"]["phenotypes"]
['ciprofloxacin']

>>> pheno3 = PhenotypeResult(ab_m2)
>>> pheno3.set_resistant(True)
>>> pheno3.add_feature(res, isolate, rm)
>>> res.add_class(cl="phenotypes", **pheno3)
>>> res["seq_variations"]["gyrA;;81;;d"]["phenotypes"]
['ciprofloxacin', 'nalidixic acid']

```
