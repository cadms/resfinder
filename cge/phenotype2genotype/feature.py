#!/usr/bin/env python3

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
        self.seq_region = Feature.na2none(seq_region)
        start = Feature.na2none(start)
        if(start):
            self.start = int(start)
        else:
            self.start = None
        self.hit = Feature.na2none(hit)
        self.isolate = Feature.na2none(isolate)

    @staticmethod
    def na2none(entry):
        if(isinstance(entry, str)):
            if(entry.upper() == "NA"):
                return None
        return entry


class Gene(Feature):
    """
    """
    def __init__(self, unique_id, seq_region=None, start=None, end=None,
                 hit=None, isolate=None):
        Feature.__init__(self, unique_id, seq_region, start, hit, isolate)
        end = Feature.na2none(end)
        if(end):
            self.end = int(end)
        else:
            self.end = None


class ResGene(Gene):
    """
    """
    def __init__(self, unique_id, seq_region=None, start=None, end=None,
                 hit=None, isolate=None, ab_class=None):
        Gene.__init__(self, unique_id, seq_region, start, end, hit, isolate)
        self.ab_class = Feature.na2none(ab_class)


class Mutation(Gene):
    """
    """
    def __init__(self, unique_id, seq_region=None, pos=None, hit=None,
                 ref_codon=None, mut_codon=None, ref_aa=None, mut_aa=None,
                 isolate=None, insertion=None, deletion=None, end=None,
                 nuc=False, premature_stop=0, frameshift=None):
        Gene.__init__(self, unique_id, seq_region, pos, end, hit, isolate)
        if(pos is not None):
            self.pos = int(pos)
        self.ref_codon = Feature.na2none(ref_codon)
        self.mut_codon = Feature.na2none(mut_codon)
        self.ref_aa = Feature.na2none(ref_aa)
        self.mut_aa = Feature.na2none(mut_aa)
        self.nuc = Feature.na2none(nuc)
        self.insertion = Feature.na2none(insertion)
        self.deletion = Feature.na2none(deletion)
        # Indicate how many percent the region was truncated.
        self.premature_stop = Feature.na2none(premature_stop)
        self.frameshift = Feature.na2none(frameshift)

        # Create mut string
        if(insertion and insertion is True):
            self.mut_string_short = (
                "{pos}_{end}ins{codon}"
                .format(pos=pos, end=end, codon=self.mut_codon.upper()))
            self.mut_string = (
                "{region}_{mut_str_short}"
                .format(region=self.seq_region,
                        mut_str_short=self.mut_string_short))
        elif(deletion and deletion is True):
            if(end):
                self.mut_string_short = (
                    "{pos}_{end}del".format(pos=pos, end=end))
            else:
                self.mut_string_short = (
                    "{pos}del{codon}".format(pos=pos, codon=ref_codon.upper()))
            self.mut_string = (
                "{region}_{mut_str_short}"
                .format(region=self.seq_region,
                        mut_str_short=self.mut_string_short))
        elif(nuc and nuc is True):
            self.mut_string_short = (
                "{pos}{ref}>{mut}"
                .format(pos=pos, ref=ref_codon.upper(), mut=mut_codon.upper()))
            self.mut_string = (str(self.seq_region) + "_"
                               + self.mut_string_short)
            self.mut_string = (
                "{region}_{mut_str_short}"
                .format(region=self.seq_region,
                        mut_str_short=self.mut_string_short))
        else:
            self.mut_string_short = (
                "{ref}{pos}{mut}"
                .format(ref=self.ref_aa.upper(), pos=pos,
                        mut=self.mut_aa.upper()))
            self.mut_string = (
                "{region}_{mut_str_short}"
                .format(region=self.seq_region,
                        mut_str_short=self.mut_string_short))


class ResMutation(Mutation):
    """
    """
    def __init__(self, unique_id, seq_region=None, pos=None, hit=None,
                 ref_codon=None, mut_codon=None, ref_aa=None, mut_aa=None,
                 isolate=None, insertion=None, deletion=None, end=None,
                 nuc=False, premature_stop=False, frameshift=None,
                 ab_class=None):
        Mutation.__init__(self, unique_id, seq_region, pos, hit, ref_codon,
                          mut_codon, ref_aa, mut_aa, isolate, insertion,
                          deletion, end, nuc, premature_stop, frameshift)
        self.ab_class = ab_class
