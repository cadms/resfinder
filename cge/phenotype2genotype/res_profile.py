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
# import io
# import urllib.parse
from .feature import Feature, Gene, Mutation


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class PhenoDB(dict):
    """ Loads a text table into dict.
        The dict consists of Phenotype objects. The keys are unique ids.
    """

    def __init__(self, acquired_file=None, point_file=None):

        # Stores non-redundant complete list of antibiotics in DB.
        self.antibiotics = {}

        if(acquired_file is None and point_file is None):
            eprint("ERROR: No pheotype database files where specified.")
            quit(1)

        if(acquired_file):
            self.load_acquired_db(acquired_file)

        if(point_file):
            self.load_point_db(point_file)

    def load_acquired_db(self, txt_file):

        with open(txt_file, "r") as fh:
            # Skip headers
            fh.readline()
            line_counter = 1

            for line in fh:
                try:
                    line_counter += 1

                    # line = line.encode("latin_1")
                    line = line.rstrip()
                    line_list = line.split("\t")
                    line_list = list(map(str.rstrip, line_list))

                    # ID in DB is <gene>_<group>_<acc>. Ex: blaB-2_1_AF189300.
                    # The acc should be unique and is used here.
                    phenodb_id = line_list[0]
                    phenodb_id = phenodb_id.split("_")
                    unique_id = phenodb_id[-1]

                    ab_class = self.get_csv_tuple(line_list[1].lower())

                    pub_phenotype = self.get_csv_tuple(line_list[2].lower())
                    if("unknown" in pub_phenotype or "none" in pub_phenotype):
                        pub_phenotype = ()

                    pmid = self.get_csv_tuple(line_list[3].lower())

                    phenotype = list(pub_phenotype)

                    if(len(line_list) > 4 and line_list[4]):
                        sug_phenotype = self.get_csv_tuple(line_list[4])
                        for p in sug_phenotype:
                            if p not in pub_phenotype:
                                phenotype.append(p)
                    else:
                        sug_phenotype = ()

                    if(len(line_list) > 5 and line_list[5]):
                        gene_class = line_list[5].lower()
                    else:
                        gene_class = None
                    if(len(line_list) > 6 and line_list[6]):
                        susceptibile = self.get_csv_tuple(line_list[6])
                    else:
                        susceptibile = ()
                    if(len(line_list) > 7 and line_list[7]):
                        notes = line_list[7]
                    else:
                        notes = ""
                    if(len(line_list) > 8 and line_list[8]):
                        species = line_list[8].lower()
                    else:
                        species = None

                    # Create non-redundant complete list of antibiotics in DB.
                    for ab in phenotype:
                        for _class in ab_class:
                            if(_class in self.antibiotics):
                                self.antibiotics[_class][ab] = True
                            else:
                                self.antibiotics[_class] = {}
                                self.antibiotics[_class][ab] = True
                    for ab in susceptibile:
                        for _class in ab_class:
                            if(_class in self.antibiotics):
                                self.antibiotics[_class][ab] = True
                            else:
                                self.antibiotics[_class] = {}
                                self.antibiotics[_class][ab] = True

                    pheno = Phenotype(unique_id, phenotype, ab_class,
                                      sug_phenotype, pub_phenotype, pmid,
                                      susceptibile=susceptibile,
                                      gene_class=gene_class, notes=notes,
                                      species=species)

                    self[unique_id] = pheno
                except IndexError:
                    eprint("Error in line " + str(line_counter))
                    eprint("Split line:\n" + str(line_list))

    def load_point_db(self, txt_file):

        with open(txt_file, "r") as fh:
            # Skip headers
            fh.readline()
            line_counter = 1

            for line in fh:
                try:
                    line_counter += 1

                    line_list = line.split("\t")
                    line_list = list(map(str.rstrip, line_list))

                    # ID in DB is just Gene ID and is not unique
                    phenodb_id = line_list[0]
                    codon_pos = line_list[3]
                    res_codon_str = line_list[6].lower()
                    unique_id = (phenodb_id + "_" + codon_pos + "_"
                                 + res_codon_str)

                    ab_class = self.get_csv_tuple(line_list[7].lower())

                    pub_phenotype = self.get_csv_tuple(line_list[8].lower())
                    if("unknown" in pub_phenotype or "none" in pub_phenotype):
                        pub_phenotype = ()

                    pmid = self.get_csv_tuple(line_list[9].lower())

                    phenotype = list(pub_phenotype)

                    if(len(line_list) > 10 and line_list[10]):
                        sug_phenotype = self.get_csv_tuple(line_list[10])
                        for p in sug_phenotype:
                            if p not in pub_phenotype:
                                phenotype.append(p)
                    else:
                        sug_phenotype = ()

                    if(len(line_list) > 11 and line_list[11]):
                        res_mechanics = line_list[11]
                    else:
                        res_mechanics = None

                    if(len(line_list) > 12 and line_list[12]):
                        notes = line_list[12]
                    else:
                        notes = ""

                    # Create non-redundant complete list of antibiotics in DB.
                    for ab in phenotype:
                        for _class in ab_class:
                            if(_class in self.antibiotics):
                                self.antibiotics[_class][ab] = True
                            else:
                                self.antibiotics[_class] = {}
                                self.antibiotics[_class][ab] = True

                    pheno = Phenotype(unique_id, phenotype, ab_class,
                                      sug_phenotype, pub_phenotype, pmid,
                                      notes=notes, res_mechanics=res_mechanics)

                    self[unique_id] = pheno

                    # DEBUG
                    print("Loaded " + unique_id)
                    print("\tPhenotype: " + str(phenotype))

                    # A pointmutation with several different res codons will
                    # never be found using all the res_codons. Instead it will
                    # be found with just one.
                    # Alternative unique ids are therefore made using just a
                    # single res_codon.
                    res_codon = self.get_csv_tuple(line_list[6].lower())
                    if(len(res_codon) > 1):
                        for codon in res_codon:
                            unique_id_alt = (phenodb_id + "_" + codon_pos
                                             + "_" + codon)
                            self[unique_id_alt] = pheno
                except IndexError:
                    eprint("Error in line " + str(line_counter))
                    eprint("Split line:\n" + str(line_list))

    @staticmethod
    def get_csv_tuple(csv_string):
        """ Takes a string containing a comma seperated list, makes it all
            lower case, remove empty entries, remove trailing and preseeding
            whitespaces, and returns the list as a tuple.
        """
        csv_string = csv_string.lower()
        string_list = csv_string.split(",")
        # Remove empty entries.
        out_list = [var.strip() for var in string_list if var]
        return tuple(out_list)

    def print_db_stats(self):
        """ Prints some stats about the database to stdout.
        """
        counter_ab_class = 0
        counter_ab = 0

        ab_output = ""

        print("-------------- LIST OF CLASSES --------------")
        for ab_class in self.antibiotics:
            counter_ab_class += 1
            print(ab_class)
        print("--------------- END OF CLASSES --------------")
        print("|")
        print("|")
        print("------------ LIST OF ANTIBIOTICS ------------")
        for ab_class in self.antibiotics:
            print(ab_class)
            for ab in self.antibiotics[ab_class]:
                counter_ab += 1
                print("\t" + ab)
        print("------------ END OF ANTIBIOTICS -------------")
        print("|")
        print("|")
        print("------------------ SUMMARY ------------------")
        print("No. of classes: " + str(counter_ab_class))
        print("No. of antibiotics: " + str(counter_ab))
        print("-------------------- END --------------------")


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
        self.res_mechanics = res_mechanics


class Antibiotics():
    """ Class is implemented to be key in a dict. The class can be tested
        against isinstances of itself and strings.
    """
    def __init__(self, name, class_, feature, published=None):
        self.name = name
        self.class_ = class_
        self.published = published
        self.features = {}
        self.features[feature.unique_id] = feature

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if isinstance(other, Antibiotics):
            return other.name == self.name
        elif isinstance(other, str):
            return other == self.name
        else:
            return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def add_feature(self, feature):
        self.features[feature.unique_id] = feature

    def get_mut_names(self, _list=False):
        names = {}
        for f in self.features:
            feature = self.features[f]
            if(not isinstance(feature, Mutation)):
                continue

            names[feature.unique_id] = [feature.seq_region,
                                        feature.ref_aa.upper(),
                                        str(feature.pos),
                                        feature.mut_aa.upper()]
        if(_list):
            return names.keys()
        else:
            return names

    def get_mut_namewannot(self, tostring=False):
        names = self.get_mut_names()

        output = []
        for unique_id in names:
            annot = names[unique_id]
            name_w_annot = annot[0] + " (" + "".join(annot[1:]) + ")"
            output.append(name_w_annot)

        if(tostring):
            return ", ".join(output)
        else:
            return output

    def get_gene_names(self, list_=False):
        names = {}
        for f in self.features:
            feature = self.features[f]
            if(not isinstance(feature, Gene)):
                continue

            if(feature.hit.name):
                names[feature.unique_id] = feature.hit.name
            else:
                names[feature.unique_id] = "Not Available"

        if(list_):
            return names.keys()
        else:
            return names

    def get_gene_namewacc(self, tostring=False):
        names = self.get_gene_names()

        output = []
        for acc in names:
            name_w_acc = names[acc] + " (" + acc + ")"
            output.append(name_w_acc)

        if(tostring):
            return ", ".join(output)
        else:
            return output

    def get_pubmed_ids(self, phenodb):
        pmids = {}
        for feature in self.features:
            phenotype = phenodb[feature.unique_id]
            for pmid in phenotype.pmid:
                pmids[pmid] = True
        return tuple(pmids.keys())


class ResProfile():
    """ Given a list of features and a PhenoDB object, an object is created
        that will contain a resistance profile based on the features given and
        the phenotypes described in the PhenoDB object.
        The class contains an add_feature method that makes it possible to add
        features after the object has been created. Per default the add_feature
        method will call the update_profile method, but this can be turned off
        if one needs to add a lot of genes and wish to call the update_profile
        method only once after adding the features.
    """
    def __init__(self, features, phenodb):
        self.phenodb = phenodb
        self.resistance = {}
        self.susceptibile = {}
        self.resistance_classes = {}
        self.missing_db_features = []

        for feature in features:
            if(feature.unique_id in phenodb):
                self.add_feature(feature, update=False)
            else:
                eprint("Not found in PhenoDB: " + feature.unique_id)
                self.missing_db_features.append(feature)
        self.update_profile()

    def add_feature(self, feature, update=True):
        """
        """
        phenotype = self.phenodb[feature.unique_id]

        for antibiotic in phenotype.pub_phenotype:

            # Create collection of features grouped with respect to ab class
            for _class in phenotype.ab_class:
                if(_class not in self.resistance_classes):
                    self.resistance_classes[_class] = set()
                self.resistance_classes[_class].add(feature)

            if(antibiotic not in self.resistance):
                self.resistance[antibiotic] = Antibiotics(antibiotic,
                                                          phenotype.ab_class,
                                                          feature,
                                                          published=True)
            else:
                if(not self.resistance[antibiotic].published):
                    self.resistance[antibiotic].published = True
                    self.resistance[antibiotic].add_feature(feature)

        for antibiotic in phenotype.sug_phenotype:

            # Create collection of features grouped with respect to ab class
            for _class in phenotype.ab_class:
                if(_class not in self.resistance_classes):
                    self.resistance_classes[_class] = set()
                self.resistance_classes[_class].add(feature)

            if(antibiotic not in self.resistance):
                self.resistance[antibiotic] = Antibiotics(antibiotic,
                                                          phenotype.ab_class,
                                                          feature,
                                                          published=False)
            else:
                self.resistance[antibiotic].add_feature(feature)

        for antibiotic in phenotype.susceptibile:
            if(antibiotic not in self.resistance):
                self.susceptibile[antibiotic] = Antibiotics(antibiotic,
                                                            phenotype.ab_class,
                                                            feature)

        if(update):
            self.update_profile()

    def update_profile(self):
        """
        """
        susc_abs = self.susceptibile.keys()
        for ab in susc_abs:
            if ab in self.resistance:
                del self.susceptibile[ab]
