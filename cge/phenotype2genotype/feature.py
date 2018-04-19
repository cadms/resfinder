#! /tools/bin/python3

import argparse
import os.path
import re
import shutil
from signal import *
import tempfile
import sys
import subprocess
# import urllib.parse
from itertools import chain


class Feature():
    """ A feature describes a location on a genome/contig.
        The 'type' variable should be used to describe the type of feature. For
        example 'gene', 'promoter' etc. It is suggested that features that only
        describes a part of a gene, promoter etc. is prefixed with "partial_"
        (e.g. 'partial_gene'). It is also suggested that features describing a
        part of the genome without anotations/function is named 'region'.
    """
    def __init__(self, unique_id, seq_region=None, start=None, hit=None,
                 isolate=None):
        self.id = unique_id
        self.unique_id = unique_id
        self.seq_region = seq_region
        if(start):
            self.start = int(start)
        else:
            self.start = None
        self.hit = hit
        self.isolate = isolate


class Gene(Feature):
    """
    """
    def __init__(self, unique_id, seq_region=None, start=None, end=None,
                 hit=None):
        Feature.__init__(self, unique_id, seq_region, start, hit)
        if(end):
            self.end = int(end)
        else:
            self.end = None


class Mutation(Feature):
    """
    """
    def __init__(self, unique_id, seq_region=None, pos=None, hit=None,
                 ref_codon=None, mut_codon=None, ref_aa=None, mut_aa=None):
        Feature.__init__(self, unique_id, seq_region, pos, hit)
        if(pos is not None):
            self.pos = int(pos)
        self.ref_codon = ref_codon
        self.mut_codon = mut_codon
        self.ref_aa = ref_aa
        self.mut_aa = mut_aa
