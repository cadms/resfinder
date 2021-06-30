#!/usr/bin/env python3

from cge.phenotype2genotype.res_profile import PhenoDB
from .gene_result import GeneResult


class SeqVariationResult(dict):
    def __init__(self, res_collection, mismatch, region_results, db_name):
        """
            Input:
                res_collection: Result object created by the cge core module.
                mismatch: Custom list of results from PointFinder
                region_results: List of GeneResult onbjects in which the
                                sequence variation was found.
                db_name: Currently always "PointFinder"

            Method loads the given res_collection with results from mismatch.
        """

        self["type"] = "seq_variation"

        self["seq_var"] = mismatch[4]
        self["ref_codon"] = mismatch[5].lower()
        self["var_codon"] = mismatch[6].lower()

        if(len(self["ref_codon"]) == 3):
            self["codon_change"] = ("{}>{}".format(self["ref_codon"],
                                                   self["var_codon"]))

        if(len(mismatch) > 7):
            self["ref_aa"] = mismatch[7].lower()
            self["var_aa"] = mismatch[8].lower()

        self["ref_start_pos"] = mismatch[1]
        self["ref_end_pos"] = mismatch[2]

        self.load_var_type(mismatch[0])

        region_name = region_results[0]["ref_id"]
        region_name = PhenoDB.if_promoter_rename(region_name)

        if(len(mismatch) > 7):
            self["ref_id"] = ("{id}{deli}{pos}{deli}{var}"
                              .format(id=region_name,
                                      pos=self["ref_start_pos"],
                                      var=self["var_aa"], deli=";;"))
        else:
            self["ref_id"] = ("{id}{deli}{pos}{deli}{var}"
                              .format(id=region_name,
                                      pos=self["ref_start_pos"],
                                      var=self["var_codon"], deli=";;"))

        self["key"] = SeqVariationResult._get_unique_seqvar_key(res_collection,
                                                                self["ref_id"])

        self["ref_database"] = res_collection.get_db_key(db_name)[0]

        region_keys = []
        for result in region_results:
            region_keys.append(result["key"])
        self["seq_regions"] = region_keys

    def load_var_type(self, type):
        """
            Input:
                type: String, must be one of: sub, ins, del

            Sets the correct bool depending on given type.
        """

        self["substitution"] = False
        self["deletion"] = False
        self["insertion"] = False
        if(type == "sub"):
            self["substitution"] = True
        elif(type == "ins"):
            self["insertion"] = True
        elif(type == "del"):
            self["deletion"] = True

    @staticmethod
    def _get_unique_seqvar_key(res_collection, minimum_key,
                               delimiter=";;"):
        """
            Input:
                res_collection: Result object created by the cge core module.
                minimum_key: The smallest key possible. If minimum_key is
                             already unique, minimum_key will be returned by the
                             method.
                delimiter: String used as delimiter inside the returned key.

            Returns a unique seq_variations key. If minimum_key is unique it
            will be returned. Else a random string will be appended after an
            additional delimiter.
        """

        unique_key = minimum_key
        while(unique_key in res_collection["seq_variations"]):
            rnd_str = GeneResult.random_string()
            unique_key = ("{key}{deli}{rnd}"
                          .format(key=minimum_key, deli=delimiter,
                                  rnd=rnd_str))

        return unique_key
