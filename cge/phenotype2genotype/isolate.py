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
from .feature import Feature, ResGene, Mutation
from .phenodbpoint import PhenoDBPoint
from .res_profile import PhenoDB, ResProfile, FeatureGroup
from .dbhit import DBHit


class Isolate(dict):
    """ An isolate class is a dict of Features.
    """
    def __init__(self, name, species=None):
        self.name = name
        self.resprofile = None
        self.species = None

    def load_resfinder_tab(self, tabbed_output):
        with open(tabbed_output, "r") as fh:
            while(True):
                line = fh.readline()
                if(not line):
                    break

                line = line.rstrip()
                if(not line):
                    continue

                ab_class = line
                second_line = fh.readline().rstrip()

                if(second_line == "No hit found"):
                    continue

                # At this point second line must be headers, and are skipped.

                res_hit = fh.readline().rstrip()

                while(res_hit):
                    hit_list = res_hit.split("\t")
                    match_length, ref_length = hit_list[2].split("/")
                    start_ref, end_ref = hit_list[4].split("..")
                    hit = DBHit(name=hit_list[0], identity=hit_list[1],
                                match_length=match_length,
                                ref_length=ref_length, start_ref=start_ref,
                                end_ref=end_ref, acc=hit_list[8],
                                db="resfinder")

                    start_feat, end_feat = hit_list[6].split("..")

                    if(start_feat == "NA"):
                        start_feat = None
                        end_feat = None

                    gene_feat = ResGene(unique_id=hit_list[8],
                                        seq_region=hit_list[5],
                                        start=start_feat, end=end_feat,
                                        hit=hit, ab_class=[ab_class.lower()])

                    if(hit_list[8] in self):
                        temp_list = self[hit_list[8]]
                        temp_list.append(gene_feat)
                        self[hit_list[8]] = temp_list
                    else:
                        self[hit_list[8]] = [gene_feat]

                    res_hit = fh.readline().rstrip()

    def load_pointfinder_tab(self, tabbed_output):
        with open(tabbed_output, "r") as fh:
            for line in fh:

                line = line.rstrip()
                if(not line):
                    continue

                headers = line

                point_hit = fh.readline().rstrip()

                while(point_hit):
                    hit_list = point_hit.split("\t")

                    # Mutation list examples:
                    # [gyrA, p.S83L]
                    # [ampC, promoter, n.-42C>T]
                    mutation_list = hit_list[0].split(" ")
                    mut_name = []
                    ins = None
                    deletion = None
                    mut_end = None
                    protein_mut = None
                    nucleotide_mut = None

                    for m in mutation_list:

                        # Protein / Amino acid mutations
                        # Ex.: "p.T83I"
                        if(m.startswith("p.")):
                            nucleotide_mut = False

                            mut_match = re.search(
                                r"^p.(\D{1})(\d+)(\D{1})$", m)
                            ref_aa = mut_match.group(1).lower()
                            pos = mut_match.group(2)
                            mut_aa = mut_match.group(3).lower()

                        # Nucleotide mutations
                        # Ex. sub: n.-42T>C
                        # Ex. ins: n.-13_-14insG    (TODO: Verify)
                        # Ex. del: n.42delT         (TODO: Verify)
                        # Ex. del: n.42_45del       (TODO: Verify)
                        elif(m.startswith("n.")):
                            nucleotide_mut = True
                            ref_aa = None
                            mut_aa = None

                            sub_match = re.search(
                                r"^n.(-{0,1}\d+)(\D{1})>(\D{1})$", m)
                            ins_match = re.search(
                                r"^n.(-{0,1}\d+)_(-{0,1}\d+)ins([CTGA]+)$", m)
                            del_match = re.search((r"^n.(-{0,1}\d+)_{0,1}(-{0,1}\d*)del[CTGA]*$"), m)
                            if(sub_match):
                                pos = sub_match.group(1)
                            elif(ins_match):
                                ins = True
                                pos = ins_match.group(1)
                                mut_end = ins_match.group(2)
                            elif(del_match):
                                deletion = True
                                pos = del_match.group(1)
                                if(del_match.group(2)):
                                    mut_end = del_match.group(2)

                        else:
                            mut_name.append(m)

                    mut_name_str = "-".join(mut_name)

                    # Codon field looks like: "TCG -> GCG" or "G -> C"
                    # TODO: Camilla, what are possible values of this field?
                    codon_match = re.search(
                        r"^(\S+) -> (\S+)$", hit_list[1].lower())
                    mut_codon = None
                    if(codon_match):
                        ref_codon = codon_match.group(1)
                        mut_codon = codon_match.group(2)
                    else:
                        sys.exit(("EEROR: isolate.py failed to match codon: "
                                  + hit_list[1].lower()))

                    if(not nucleotide_mut):
                        unique_id = mut_name_str + "_" + pos + "_" + mut_aa
                    elif(mut_codon):
                        unique_id = mut_name_str + "_" + pos + "_" + mut_codon
                    elif(deletion):
                        unique_id = (mut_name_str + "_" + pos
                                     + "_del" + ref_codon)
                    else:
                        sys.exit(("ERROR: isolate.py was not able to infer a "
                                  "unique ID. 'mutation_list' contains: "
                                  + str(mutation_list)))

                    mut_feat = Mutation(unique_id=unique_id,
                                        seq_region=mut_name_str,
                                        pos=pos, ref_codon=ref_codon,
                                        mut_codon=mut_codon, ref_aa=ref_aa,
                                        mut_aa=mut_aa,
                                        isolate=self,
                                        insertion=ins,
                                        deletion=deletion,
                                        end=mut_end,
                                        nuc=nucleotide_mut)

                    if(unique_id in self):
                        temp_list = self[unique_id]
                        temp_list.append(mut_feat)
                        self[unique_id] = temp_list
                    else:
                        self[unique_id] = [mut_feat]

                    try:
                        point_hit = fh.readline().rstrip()
                    except StopIteration:
                        point_hit = None

    def calc_res_profile(self, phenodb):
        """
        """
        features = self.values()
        # Flatten features.
        features = list(chain.from_iterable(features))
        self.resprofile = ResProfile(features, phenodb)

    def profile_to_str_table(self, header=False):
        """
        """
        output_str = ""

        if(header):
            output_str = (
                "# ResFinder phenotype results.\n"
                "# Sample: " + self.name + "\n"
                "# \n"
                "# The phenotype 'No resistance' should be interpreted with\n"
                "# caution, as it only means that nothing in the used\n"
                "# database indicate resistance, but resistance could exist\n"
                "# from 'unknown' or not yet implemented sources.\n"
                "# \n"
                "# The 'Match' column stores one of the integers 0, 1, 2, 3.\n"
                "#      0: No match found\n"
                "#      1: Match < 100% ID AND match length < ref length\n"
                "#      2: Match = 100% ID AND match length < ref length\n"
                "#      3: Match = 100% ID AND match length = ref length\n"
                "# If several hits causing the same resistance are found,\n"
                "# the highest number will be stored in the 'Match' column.\n"
                "\n"
            )
            output_str += ("# Antimicrobial\t"
                           "Class\t"
                           "WGS-predicted phenotype\t"
                           "Match\t"
                           "Genetic background\n")

        # For each antibiotic class
        for ab_class in self.resprofile.phenodb.antibiotics.keys():

            # For each antibiotic in current class
            for ab_name in self.resprofile.phenodb.antibiotics[ab_class]:
                output_str += ab_name + "\t" + ab_class

                # Isolate is resistant towards the antibiotic
                if(ab_name in self.resprofile.resistance):
                    ab = self.resprofile.resistance[ab_name]

                    if(ab.published):
                        output_str += "\tResistant"
                    else:
                        output_str += "\tResistant*"

                    # Find the resistance causing gene with the best match
                    # Mutations will always have best match as they are only
                    # either present or absent.
                    best_match = 0
                    for unique_id in ab.features:
                        feature = ab.features[unique_id]
                        if(isinstance(feature, Mutation)
                           or isinstance(feature, FeatureGroup)):
                            best_match = 3
                        # Note: Mutations do not have "hits"
                        elif(feature.hit.match_category > best_match):
                            best_match = feature.hit.match_category

                    output_str += "\t" + str(best_match)

                    gene_list = ab.get_gene_namewacc(tostring=True)
                    mut_list = ab.get_mut_namewannot(tostring=True)

                    if(gene_list and mut_list):
                        gene_mut_str = gene_list + " " + mut_list
                    elif(gene_list):
                        gene_mut_str = gene_list
                    elif(mut_list):
                        gene_mut_str = mut_list
                    else:
                        gene_mut_str = ""

                    output_str += "\t" + gene_mut_str + "\n"

                # Isolate is susceptibile towards the antibiotic
                elif(ab_name in self.resprofile.susceptibile):
                    ab = self.resprofile.susceptibile[ab_name]
                    # Genetic background is not written if susceptibile
                    gene_list = ""
                    # Uncomment next line to write genetic background for susc.
                    # gene_list = ab.get_gene_namewacc(tostring=True)
                    output_str += "\tNo resistance\t0\t" + gene_list + "\n"
                else:
                    output_str += "\tNo resistance\t0\t\n"

        if(self.resprofile.missing_db_features):
            output_str += ("\n# WARNING: Missing features from phenotype "
                           "database:\n")
            output_str += "# Feature_ID\tRegion\tDatabase\tHit\n"

            for feature in self.resprofile.missing_db_features:
                output_str += (feature.unique_id + "\t"
                               + feature.seq_region + "\t")

                if(feature.hit is None):
                    output_str += "\t"
                else:
                    output_str += str(feature.hit.db) + "\t" + feature.hit.name

                output_str += "\n"

        return output_str
