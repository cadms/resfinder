#!/usr/bin/env python3
import random
import string

from .exceptions import DuplicateKeyError
from cge.phenotype2genotype.res_profile import PhenoDB
from .gene_result import GeneResult
from .seq_variation_result import SeqVariationResult
from .phenotype_result import PhenotypeResult

from cge.phenotype2genotype.feature import ResGene, ResMutation

import json


class ResFinderResultHandler():

    @staticmethod
    def standardize_results(res_collection, res, ref_db_name):
        """
            Input:
                res_collection: Result object created by the cge core module.
                res: Custom dictionary of results from ResFinder
                ref_db_name: 'ResFinder' or 'PointFinder'

            Method loads the given res_collection with results from res.
        """
        for db_name, db in res.items():
            if(db_name == "excluded"):
                continue

            if(db == "No hit found"):
                continue

            for unique_id, hit_db in db.items():
                if(unique_id in res["excluded"]):
                    continue

                gene_result = GeneResult(res_collection, hit_db, ref_db_name)

                if gene_result["key"] not in res_collection["seq_regions"]:
                    res_collection.add_class(cl="seq_regions", **gene_result)
                else:
                    raise DuplicateKeyError(
                        "About to overwrite dict entry. This should not be "
                        "happening as all keys are supposed to be unique."
                        "Non-unique key was: {}".format(gene_result["key"]))

    @staticmethod
    def load_res_profile(res_collection, isolate):
        """
            Input:
                res_collection: Result object created by the cge core module.
                isolate: Isolate object

            Method loads the given res_collection with results from res.
        """
        # For each antibiotic class
        for ab_class in isolate.resprofile.phenodb.antibiotics.keys():
            # For each antibiotic in current class
            for phenodb_ab in isolate.resprofile.phenodb.antibiotics[ab_class]:

                phenotype = PhenotypeResult(phenodb_ab)

                # Isolate is resistant towards the antibiotic
                if(phenodb_ab in isolate.resprofile.resistance):
                    phenotype.set_resistant(True)

                    isolate_ab = isolate.resprofile.resistance[phenodb_ab]

                    for unique_id, feature in isolate_ab.features.items():

                        if(isinstance(feature, ResGene)
                           or isinstance(feature, ResMutation)):

                            phenotype.add_feature(res_collection, isolate,
                                                  feature)

                res_collection.add_class(cl="phenotypes", **phenotype)


class PointFinderResultHandler():

    @staticmethod
    def standardize_results(res_collection, res, ref_db_name):
        """
            Input:
                res_collection: Result object created by the cge core module.
                res: Custom dictionary of results from PointFinder
                ref_db_name: 'ResFinder' or 'PointFinder'

            Method loads the given res_collection with results from res.
        """

        for gene_name, db in res.items():

            # Ignore information in excluded dict
            if(gene_name == "excluded"):
                continue

            # Ignore genes found in excluded dict
            if gene_name in res["excluded"]:
                continue

            if(isinstance(db, str)):
                if(db == "No hit found"):
                    continue
                if db.startswith("Gene found with coverage"):
                    continue

            gene_results = []

            # For BLAST results
            db_hits = db.get("hits", {})

            # For KMA results
            if(not db_hits):
                id = db["sbjct_header"]
                db_hits[id] = db

            for unique_id, hit_db in db_hits.items():

                # Ignore genes found in excluded dict
                if(unique_id in res["excluded"]):
                    continue

                gene_result = GeneResult(res_collection, hit_db, ref_db_name)
                res_collection.add_class(cl="seq_regions", **gene_result)
                gene_results.append(gene_result)

            mismatches = db["mis_matches"]

            for mismatch in mismatches:
                seq_var_result = SeqVariationResult(
                    res_collection, mismatch, gene_results, ref_db_name)
                res_collection.add_class(cl="seq_variations", **seq_var_result)
