#!/usr/bin/env python3
import random
import string


class GeneResult(dict):
    def __init__(self, res_collection, res):
        self["type"] = "gene"
        self["ref_id"] = res["sbjct_header"]
        self["name"], self.variant, self["ref_acc"] = (
            GeneResult._split_sbjct_header(self["ref_id"]))

        self["ref_start_pos"] = res["sbjct_start"]
        self["ref_end_pos"] = res["sbjct_end"]
        self["key"] = self._get_unique_gene_key(res_collection)

    @staticmethod
    def _split_sbjct_header(header):
        sbjct = header.split("_")
        template = sbjct[0]
        variant = sbjct[1]
        acc = "_".join(sbjct[2:])
        return (template, variant, acc)

    def _get_unique_gene_key(self, res_collection, delimiter=";;"):
        pre_start = "start{}".format(self["ref_start_pos"])
        pre_end = "end{}".format(self["ref_end_pos"])

        gene_key = (
            "{name}{deli}{var}{deli}{ref_acc}{deli}"
            "{ref_start_pos}{deli}{ref_end_pos}{deli}"
            .format(deli=delimiter, var=self.variant, **self))

        # Attach random string if key already exists
        minimum_gene_key = gene_key
        while(gene_key in res_collection["genes"]):
            rnd_str = GeneResult.randomString()
            gene_key = ("{key}{deli}{rnd}{deli}"
                        .format(key=minimum_gene_key, deli=delimiter,
                                rnd=rnd_str))

        return gene_key

    @staticmethod
    def randomString(stringLength=4):
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(stringLength))


class ResFinderResultHandler():

    @staticmethod
    def standardize_results(res_collection, res):
        for db_name, db in res.items():
            if(db_name == "excluded"):
                continue

            if(db == "No hit found"):
                continue

            for unique_id, hit_db in db.items():
                if(unique_id in res["excluded"]):
                    continue

                # unique_key = GeneResult._get_unique_gene_key(prefix,
                #                                             res_collection)
