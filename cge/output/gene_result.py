#!/usr/bin/env python3
import random
import string

from cge.phenotype2genotype.res_profile import PhenoDB


class GeneResult(dict):
    def __init__(self, res_collection, res, db_name):
        """
            Input:
                res_collection: Result object created by the cgelib package.
                res: Custom dictionary containing information about a single hit
                     from ResFinder.
                ref_db_name: 'ResFinder' or 'PointFinder'

            Method creates a seq_region dict as defined in the BeOne template:
            https://bitbucket.org/genomicepidemiology/cgelib/src/master/src/
            cgelib/output/templates_json/beone/
            from res.
        """

        self.db_name = db_name
        self["type"] = "seq_region"
        self["gene"] = True

        self["ref_id"] = res["sbjct_header"]
        self["ref_id"] = PhenoDB.if_promoter_rename(self["ref_id"])

        if(db_name == "ResFinder"):
            self["name"], self.variant, self["ref_acc"] = (
                GeneResult._split_sbjct_header(self["ref_id"]))
        elif(db_name == "PointFinder"):
            self["name"] = self["ref_id"]

        self["key"] = self._get_unique_gene_key(res_collection)

        self["identity"] = res["perc_ident"]
        self["alignment_length"] = res["HSP_length"]
        self["ref_seq_lenght"] = res["sbjct_length"]
        self["depth"] = res.get("depth", None)
        self["ref_start_pos"] = res["sbjct_start"]
        self["ref_end_pos"] = res["sbjct_end"]
        self["query_id"] = res["contig_name"]
        self["query_start_pos"] = res["query_start"]
        self["query_end_pos"] = res["query_end"]
        self["ref_database"] = res_collection.get_db_key(db_name)[0]

        # BLAST coverage formatted results
        coverage = res.get("coverage", None)
        if(coverage is None):
            # KMA coverage formatted results
            coverage = res["perc_coverage"]
        else:
            coverage = float(coverage) * 100
        self["coverage"] = coverage

        self.remove_NAs()

    def remove_NAs(self):
        """
            Remove all entries containing NA og None as values.

            Removing None is not necessary as the Result object will ignore all
            entries with None values.
        """
        na_keys = []
        for key, val in self.items():
            if(val == "NA" or val is None):
                na_keys.append(key)
        for key in na_keys:
            del self[key]

    @staticmethod
    def get_rnd_unique_gene_key(gene_key, res_collection,
                                minimum_gene_key, delimiter):
        """
            Input:
                gene_key: None-unique key
                res_collection: Result object created by the cgelib package.
                minimum_key: Key prefix
                delimiter: String used as delimiter inside the returned key.
            Output:
                gene_key: Unique key (string)

            If gene_key is found in res_collection. Creates a unique key by
            appending a random string ton minimum_gene_key.
        """
        while(gene_key in res_collection["seq_regions"]):
            rnd_str = GeneResult.random_string(str_len=4)
            gene_key = ("{key}{deli}{rnd}"
                        .format(key=minimum_gene_key, deli=delimiter,
                                rnd=rnd_str))
        return gene_key

    @staticmethod
    def random_string(str_len=4):
        """
            Output:
                random string of length 'str_len'

            Return a random string of the provided length. The string will only
            consist of lowercase ascii letters.
        """
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(str_len))

    @staticmethod
    def _split_sbjct_header(header):
        """
            Input:
                header: database entry header (ref_header/subject_header)
            Output:
                template: name of entry (string)
                variant: Variant interger (string) or None
                acc: Accession number given by sequnce database (string) or None

            Splits the input header by underscores and returns first list item
            as template. If list is > 1 then variant and acc will also return
            strings. If not they return None.
        """
        sbjct = header.split("_")
        template = sbjct[0]

        if(len(sbjct) > 1):
            variant = sbjct[1]
            acc = "_".join(sbjct[2:])
        else:
            variant = None
            acc = None

        return (template, variant, acc)

    def _get_unique_gene_key(self, res_collection, delimiter=";;"):
        """
            Input:
                res_collection: Result object created by the cgelib package.
                delimiter: String used as delimiter inside the returned key.
            Output:
                key: Unique key for hit

            Creates a unique key for GeneResult instance. Key format depends on
            database. If gene result is considered indentical to an existing
            gene result in the provided res_collection, it will not create a new
            key. Two restults are considered identical if they have the same
            query_id, query_start_pos and query_end_pos. If it is mapping
            results the query_* doesn't exist, and results will never be
            considered identical.
        """
        if(self.db_name == "ResFinder"):
            gene_key = ("{name}{deli}{var}{deli}{ref_acc}"
                        .format(deli=delimiter, var=self.variant, **self))

        if(self.db_name == "PointFinder"):
            gene_key = self["name"]

        # Attach random string if key already exists
        minimum_gene_key = gene_key
        if gene_key in res_collection["seq_regions"]:

            query_id = self.get("query_id", "NA")

            if(query_id == "NA"):
                gene_key = GeneResult.get_rnd_unique_gene_key(
                    gene_key, res_collection, minimum_gene_key, delimiter)

            elif (query_id
                    != res_collection["seq_regions"][gene_key]["query_id"]
                  or self["query_start_pos"]
                    != res_collection["seq_regions"][gene_key]["query_start_pos"]
                  or self["query_end_pos"]
                    != res_collection["seq_regions"][gene_key]["query_end_pos"]):
                gene_key = GeneResult.get_rnd_unique_gene_key(
                    gene_key, res_collection, minimum_gene_key, delimiter)

        return gene_key
