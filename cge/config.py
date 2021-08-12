#!/usr/bin/env python3

import os.path
import os
import sys
import subprocess

from cgelib.utils.loaders_mixin import LoadersMixin

from cge.pointfinder import PointFinder


class Config():

    DEFAULT_VALS = {
        "inputfasta": None,
        "inputfastq": None,
        "outputPath": None,
        "blastPath": "blastn",
        "kmaPath": "kma",
        "species": None,
        "ignore_missing_species": False,
        "db_path_res": None,
        "db_path_res_kma": None,
        "databases": None,
        "acquired": None,
        "acq_overlap": 30,
        "min_cov": 0.6,
        "threshold": 0.8,
        "point": None,
        "db_path_point": None,
        "db_path_point_kma": "kma",
        "specific_gene": None,
        "unknown_mut": False,
        "min_cov_point": 0.01,
        "threshold_point": 0.8,
        "ignore_indels": False,
        "ignore_stop_codons": False,
        "pickle": False
    }

    def __init__(self, args):

        # Directoy of config.py substracted the last dir 'cge'
        self.resfinder_root = os.path.dirname(os.path.realpath(__file__))[:-3]
        self.readme_file = "{}{}".format(self.resfinder_root, "README.md")
        Config.set_default_and_env_vals(args, self.readme_file)

        self.set_general_opts(args)
        if(self.acquired):
            self.set_resfinder_opts(args)
        if(self.point):
            self.set_pointfinder_opts(args)
        self.set_phenotype_opts(args)

        if self.acquired is False and self.point is False:
            sys.exit("Please specify to look for acquired resistance genes, "
                     "chromosomal mutaitons or both!\n")

    def set_general_opts(self, args):
        self.outputPath = os.path.abspath(args.outputPath)
        os.makedirs(self.outputPath, exist_ok=True)

        self.acquired = bool(args.acquired)
        self.point = bool(args.point)
        self.species = Config.get_species(args.species, self.readme_file)

        if(args.inputfasta):
            self.set_fasta_related_opts(args)

        if(args.inputfastq):
            self.set_fastq_related_opts(args)

        amr_abbreviations_file = ("{}/amr_abbreviations.md"
                                  .format(self.resfinder_root))
        self.amr_abbreviations = LoadersMixin.load_md_table_after_keyword(
            amr_abbreviations_file, "## Abbreviations")

        self.pickle = args.pickle

    @staticmethod
    def get_species(in_species, species_def_filepath):
        out_species = in_species
        if(in_species is not None and in_species.lower() == "other"):
            out_species = None
        elif(in_species is not None):
            out_species = in_species.lower()

        species_transl = LoadersMixin.load_md_table_after_keyword(
            species_def_filepath, "#### Species Abbreviations Table",
            header_key="Abbreviation")

        fixed_species = species_transl.get(out_species, None)
        if(fixed_species):
            out_species = fixed_species[0]

        return out_species

    def set_fasta_related_opts(self, args):
        self.inputfastq_1 = None
        self.inputfasta = self.get_abs_path_and_check(args.inputfasta)

        self.outPath_res_blast = "{}/resfinder_blast".format(self.outputPath)
        os.makedirs(self.outPath_res_blast, exist_ok=True)

        self.outPath_point_blast = ("{}/pointfinder_blast"
                                    .format(self.outputPath))
        os.makedirs(self.outPath_point_blast, exist_ok=True)

        self.sample_name = os.path.basename(self.inputfasta)
        self.method = PointFinder.TYPE_BLAST
        self.blast = self.get_prg_path(args.blastPath)
        self.kma = None

    def set_fastq_related_opts(self, args):
        self.inputfasta = None
        self.inputfastq_1 = self.get_abs_path_and_check(
            args.inputfastq[0])
        if(len(args.inputfastq) == 2):
            self.inputfastq_2 = self.get_abs_path_and_check(
                args.inputfastq[1])
        elif(len(args.inputfastq) > 2):
            sys.exit("ERROR: More than 2 files were provided to inputfastq: {}."
                     .format(inputfastq))
        else:
            self.inputfastq_2 = None

        self.outPath_res_kma = "{}/resfinder_kma".format(self.outputPath)
        os.makedirs(self.outPath_res_kma, exist_ok=True)

        self.outPath_point_kma = "{}/pointfinder_kma".format(self.outputPath)
        os.makedirs(self.outPath_point_kma, exist_ok=True)

        self.sample_name = os.path.basename(args.inputfastq[0])
        self.method = PointFinder.TYPE_KMA
        self.kma = self.get_prg_path(args.kmaPath)

    def set_resfinder_opts(self, args):
        self.databases = args.databases
        self.set_path_resfinderdb(args)
        self.db_config_file = "{}/config".format(self.db_path_res)
        self.db_notes_file = "{}/notes.txt".format(self.db_path_res)
        if not os.path.exists(self.db_config_file or self.db_notes_fil):
            sys.exit("Input Error: The database config or notes.txt file could "
                     "not be found in the ResFinder database directory.")

        args.min_cov = float(args.min_cov)
        args.threshold = float(args.threshold)

        # Check if coverage/identity parameters are valid
        if(args.min_cov > 1.0 or args.min_cov < 0.0):
            sys.exit("ERROR: Minimum coverage above 1 or below 0 is not "
                     "allowed. Please select a minimum coverage within the "
                     "range 0-1 with the flag -l. Given value: {}."
                     .format(args.min_cov))
        self.rf_gene_cov = args.min_cov

        if(args.threshold > 1.0 or args.threshold < 0.0):
            sys.exit("ERROR: Threshold for identity of ResFinder above 1 or "
                     "below 0 is not allowed. Please select a threshold for "
                     "identity within the range 0-1 with the flag -t. Given "
                     "value: {}.".format(args.threshold))
        self.rf_gene_id = args.threshold

        self.rf_overlap = int(args.acq_overlap)

    def set_pointfinder_opts(self, args):
        if(not self.species and not args.ignore_missing_species):
            sys.exit("ERROR: Chromosomal point mutations cannot be located if "
                     "no species has been provided. Please provide species "
                     "using the --species option.")
        elif(not self.species and args.ignore_missing_species):
            self.point = False
            return

        self.set_path_pointdb(args)
        if(self.db_path_point is None):
            self.point = False
            return

        self.specific_gene = args.specific_gene

        args.min_cov_point = float(args.min_cov_point)
        args.threshold_point = float(args.threshold_point)

        # Check if coverage/identity parameters are valid
        if(args.min_cov_point > 1.0 or args.min_cov_point < 0.0):
            sys.exit("ERROR: Minimum coverage above 1 or below 0 is not "
                     "allowed. Please select a minimum coverage within the "
                     "range 0-1 with the flag -l. Given value: {}."
                     .format(args.min_cov_point))
        self.pf_gene_cov = args.min_cov_point

        if(args.threshold_point > 1.0 or args.threshold_point < 0.0):
            sys.exit("ERROR: Threshold for identity of ResFinder above 1 or "
                     "below 0 is not allowed. Please select a threshold for "
                     "identity within the range 0-1 with the flag -t. Given "
                     "value: {}.".format(args.threshold_point))
        self.pf_gene_id = args.threshold_point

        self.unknown_mut = args.unknown_mut
        self.ignore_indels = args.ignore_indels
        self.ignore_stop_codons = args.ignore_stop_codons

    def set_phenotype_opts(self, args):
        self.point_file = None
        if(self.point):
            self.point_file = ("{}/resistens-overview.txt"
                               .format(self.db_path_point))
            _ = self.get_abs_path_and_check(self.point_file)

        self.abclassdef_file = ("{}/antibiotic_classes.txt"
                                .format(self.db_path_res))
        _ = self.get_abs_path_and_check(self.abclassdef_file)

        self.phenotype_file = ("{}/phenotypes.txt".format(self.db_path_res))
        _ = self.get_abs_path_and_check(self.phenotype_file)

    @staticmethod
    def get_abs_path_and_check(path, allow_exit=True):
        abs_path = os.path.abspath(path)
        if(not os.path.isfile(abs_path) and not os.path.isdir(abs_path)):
            if(allow_exit):
                sys.exit("ERROR: Path not found: {}".format(path))
            else:
                raise FileNotFoundError
        return abs_path

    def set_path_pointdb(self, args):
        tmp_list = self.species.split()
        if(len(tmp_list) != 1 and len(tmp_list) != 2):
            sys.exit("ERROR: Species name must contain 1 or 2 names. Given "
                     "value: {}".format(self.species))

        if(len(tmp_list) == 2):
            tmp_species_dir = "_".join(tmp_list)
        else:
            tmp_species_dir = tmp_list[0]

        if(args.db_path_point is None):
            args.db_path_point = "{}{}".format(self.resfinder_root,
                                               "/db_pointfinder")

        path_pointdb = os.path.abspath(args.db_path_point)

        self.species_dir = self._parse_species_dir(path_pointdb,
                                                   tmp_species_dir,
                                                   args.ignore_missing_species)

        if(self.species_dir is not None):
            self.db_path_point_root = path_pointdb
            self.db_path_point = "{}/{}".format(path_pointdb, self.species_dir)

    def set_path_resfinderdb(self, args):
        self.db_path_res = args.db_path_res
        if(self.db_path_res is None):
            self.db_path_res = "{}{}".format(self.resfinder_root,
                                             "/db_resfinder")
        self.db_path_res_kma = args.db_path_res_kma
        if(self.db_path_res_kma is None):
            self.db_path_res_kma = self.db_path_res

        try:
            self.db_path_res = Config.get_abs_path_and_check(
                self.db_path_res, allow_exit=False)
        except FileNotFoundError:
            sys.exit("Could not locate ResFinder database path: {}"
                     .format(self.db_path_res))

        try:
            self.db_path_res_kma = Config.get_abs_path_and_check(
                self.db_path_res_kma, allow_exit=False)
        except FileNotFoundError:
            if(self.acquired and self.inputfastq_1):
                sys.exit("Could not locate ResFinder database index path: {}"
                         .format(self.db_path_res_kma))
            else:
                pass

    def _parse_species_dir(self, path_pointdb, species_dir,
                           ignore_missing_species):
        # Check if a database for species exists
        point_dbs = PointFinder.get_db_names(path_pointdb)
        if(species_dir not in point_dbs):
            # If no db for species is found check if db for genus is found
            # and use that instead
            tmp_list = self.species.split()
            if(tmp_list[0] in point_dbs):
                species_dir = tmp_list[0]
            elif(ignore_missing_species):
                self.species = None
                self.db_path_point_root = None
                self.db_path_point = None
                species_dir = None
            else:
                sys.exit("ERROR: species '{}' ({}) does not seem to exist"
                         " as a PointFinder database."
                         .format(self.species, species_dir))
        return species_dir

    @staticmethod
    def get_prg_path(prg_path):
        try:
            prg_path = Config.get_abs_path_and_check(prg_path,
                                                     allow_exit=False)
        except FileNotFoundError:
            pass

        try:
            _ = subprocess.check_output([prg_path, "-h"])
        except PermissionError:
            sys.exit("ERROR: Missing permission. Unable to execute blastn from "
                     "the path: {}".format(prg_path))
        return prg_path

    @staticmethod
    def set_default_and_env_vals(args, env_def_filepath):

        known_envs = LoadersMixin.load_md_table_after_keyword(
            env_def_filepath, "#### Environment Variables Table")

        # Set flag values defined in environment variables
        for var, entries in known_envs.items():

            try:
                cli_val = getattr(args, entries[0])
                # Flags set by user will not be None, default vals will be None.
                if(cli_val is not None):
                    continue

                var_val = os.environ.get(var, None)
                if(var_val is not None):
                    setattr(args, entries[0], var_val)

            except AttributeError:
                sys.exit("ERROR: A flag set in the Environment Variables Table "
                         "in the README file did not match any valid flags in "
                         "ResFinder. Flag not recognized: {}."
                         .format(entries[0]))

        Config._set_default_values(args)

    @staticmethod
    def _set_default_values(args):
        for flag, def_val in Config.DEFAULT_VALS.items():
            val = getattr(args, flag)
            if(val is None):
                setattr(args, flag, def_val)
