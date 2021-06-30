#!/usr/bin/env python3


from cge.phenotype2genotype.feature import ResGene, ResMutation


class PhenotypeResult(dict):
    def __init__(self, antibiotic):
        self["type"] = "phenotype"
        self["category"] = "amr"
        self["key"] = antibiotic.name
        self["amr_classes"] = antibiotic.classes
        self["amr_resistance"] = antibiotic.name
        self["amr_resistant"] = False

    def set_resistant(self, res):
        self["amr_resistant"] = res

    def add_feature(self, res_collection, isolate, feature):
        """
            Input:
                res_collection: Result object created by the cgelib package.
                isolate: Isolate object. Must have been loaded with at
                         resistance profile using Isolate.load_finder_results
                         and then Isolate.calc_res_profile.
                feature: Either ResGene or ResMutation object (inherit feature)
        """
        # Get all keys in the result that matches the feature in question.
        # Most of the time this will be a one to one relationship.
        # However if several identical features has been found in a sample,
        # they will all have different keys, but identical ref ids.

        ref_id, type = PhenotypeResult.get_ref_id_and_type(feature, isolate)
        feature_keys = PhenotypeResult.get_keys_matching_ref_id(
            ref_id, res_collection[type])
        # Add keys to phenotype results
        pheno_feat_keys = self.get(type, [])
        pheno_feat_keys = pheno_feat_keys + feature_keys
        self[type] = pheno_feat_keys

        # Add phenotype keys to feature results
        features = res_collection[type]
        for feat_key in feature_keys:
            feat_result = features[feat_key]
            pheno_keys = feat_result.get("phenotypes", [])
            pheno_keys.append(self["key"])
            feat_result["phenotypes"] = pheno_keys

        if(type == "seq_regions"):
            db_key = res_collection.get_db_key("ResFinder")[0]
        elif(type == "seq_variations"):
            db_key = res_collection.get_db_key("PointFinder")[0]

        self["ref_database"] = db_key

    @staticmethod
    def get_ref_id_and_type(feature, isolate):
        """
            Input:
                feature: Either ResGene or ResMutation object (inherit feature)
                isolate: Isolate object. Must have been loaded with at
                         resistance profile using Isolate.load_finder_results
                         and then Isolate.calc_res_profile.
            Output (tuple):
                ref_id: id to identity the feature in relevant database
                type: 'seq_regions' or 'seq_variations'
        """
        type = None
        ref_id = None
        if(isinstance(feature, ResGene)):
            type = "seq_regions"
            ref_id = isolate.resprofile.phenodb.id_to_idwithvar[
                feature.unique_id]
        elif(isinstance(feature, ResMutation)):
            type = "seq_variations"
            ref_id = feature.unique_id
        return (ref_id, type)

    @staticmethod
    def get_keys_matching_ref_id(ref_id, res_collection):
        """
            Input:
                ref_id: Key/ID to search for
                res_collection: Result object created by the cgelib package.
                                The Result object used with this method found is
                                within the main Result object and is accessed
                                either by main_result["genes"] or by
                                main_result["seq_variations"].
            Output (list):
                out_keys: Returns a list of keys where the
                          input ref_id == res_collection[key]["ref_id"].
                          No matches will return an empty list.
        """
        out_keys = []
        for key, results in res_collection.items():
            if(ref_id == results["ref_id"]):
                out_keys.append(key)

        return out_keys
