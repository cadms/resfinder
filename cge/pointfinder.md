# PointFinder class tests

This is only testing a very small fraction of the PointFinder class.

## Setup

```python

>>> from cge.pointfinder import PointFinder

```

## _remove_insertions_in_seq(sbjct_seq, qry_seq)

```python

>>> sbj_seq = "-GGAGAA-TCCTGAC---GTTAAAATTGATAAGCTAA-"
>>> qry_seq = "XGGAGAAXTCCTGACXXXGTTAAAATTGATAAGCTAAX"
>>> (sbj_seq, qry_seq) = PointFinder._remove_insertions_in_seq(sbj_seq, qry_seq)
>>> assert(sbj_seq == qry_seq)

```
