#! /tools/bin/python3

from __future__ import print_function
import argparse
import os.path
import re
import shutil
from signal import *
import tempfile
import sys
import subprocess


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Phenotype():
    """ A Phenotype object describes the antibiotics a feature/gene causes
        resistance and susceptibility against.
        unique_id: the id is used to locate the specified phenotype.
        phenotype: Tuple of antibiotics that gene causes resistance toward.
        ab_class: Class of resistance (e.g. Beta-Lactamase).
        pub_phenotype: Tuple of published resistance (antibiotics).
        pmid: Tuple of PubMed IDs of publications presenting resistance and
              susceptibility.
        susceptibility: Tuple of antibiotics that the gene is known to be
                        susceptibil towards.
        gene_class: resistance gene class (e.g. class D)
        notes: String containing other information on the resistance gene.
        species: Species exceptions, where resistance is not observed. NOTE:
                 This information is not yet used for anything.
    """
    def __init__(self, unique_id, phenotype, ab_class, sug_phenotype,
                 pub_phenotype, pmid, susceptibile=(), gene_class=None,
                 notes="", species=None, res_mechanics=None):
        self.unique_id = unique_id
        self.phenotype = phenotype
        self.ab_class = ab_class

        self.sug_phenotype = sug_phenotype
        self.pub_phenotype = pub_phenotype
        self.pmid = pmid
        self.susceptibile = susceptibile
        self.gene_class = gene_class
        self.notes = notes
        self.species = species
