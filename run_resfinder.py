#!/usr/bin/env python3
import sys
import os
import subprocess
from argparse import ArgumentParser
import pickle
import json

from cgelib.output.result import Result
from cgelib.utils.loaders_mixin import LoadersMixin

from cge.config import Config
from cge.resfinder import ResFinder
from cge.pointfinder import PointFinder
from cge.output.std_results import ResFinderResultHandler
from cge.output.std_results import PointFinderResultHandler

#  Modules used to create the extended ResFinder output (phenotype output)
from cge.phenotype2genotype.isolate import Isolate
from cge.phenotype2genotype.res_profile import PhenoDB
from cge.phenotype2genotype.res_sumtable import ResSumTable
from cge.phenotype2genotype.res_sumtable import PanelNameError

# TODO list:
# TODO: JSON output summary is not species dependent


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser(allow_abbrev=False)

# General options
parser.add_argument("-ifa", "--inputfasta",
                    help="Input fasta file.",
                    default=None)
parser.add_argument("-ifq", "--inputfastq",
                    help="Input fastq file(s). Assumed to be single-end fastq \
                          if only one file is provided, and assumed to be \
                          paired-end data if two files are provided.",
                    nargs="+",
                    default=None)
parser.add_argument("-o", "--outputPath",
                    help=("Output directoy. If it doesn't exist, it will be "
                          "created."),
                    required=True,
                    default=None)
parser.add_argument("-b", "--blastPath",
                    help="Path to blastn",
                    default=None)
parser.add_argument("-k", "--kmaPath",
                    help="Path to KMA",
                    default=None)
parser.add_argument("-s", "--species",
                    help="Species in the sample",
                    default=None)
parser.add_argument("--ignore_missing_species",
                    action="store_true",
                    help="If set, species is provided and --point flag is set, "
                         "will not throw an error if no database is found for "
                         "the provided species. If species is not found. Point "
                         "mutations will silently be ignored.",
                    default=False)

# Acquired resistance options
parser.add_argument("-db_res", "--db_path_res",
                    help=("Path to the databases for ResFinder. Defaults to "
                          "'db_resfinder' in the ResFinder application "
                          "directory."),
                    default=None)
parser.add_argument("-db_res_kma", "--db_path_res_kma",
                    help=("Path to the ResFinder databases indexed with KMA. "
                          "Defaults to the value of the --db_res flag."),
                    default=None)
parser.add_argument("-d", "--databases",
                    help="Databases chosen to search in - if none is specified\
                          all is used",
                    default=None)
parser.add_argument("-acq", "--acquired",
                    action="store_true",
                    help="Run resfinder for acquired resistance genes",
                    default=None)
parser.add_argument("-ao", "--acq_overlap",
                    help="Genes are allowed to overlap this number of\
                          nucleotides. Default: {}.".format(
                        Config.DEFAULT_VALS["acq_overlap"]),
                    type=int,
                    default=None)
parser.add_argument("-l", "--min_cov",
                    help=("Minimum (breadth-of) coverage of ResFinder within "
                          "the range 0-1."),
                    type=float,
                    default=None)
parser.add_argument("-t", "--threshold",
                    help=("Threshold for identity of ResFinder within the "
                          "range 0-1."),
                    type=float,
                    default=None)

# Point resistance option
parser.add_argument("-c", "--point",
                    action="store_true",
                    help="Run pointfinder for chromosomal mutations",
                    default=None)
parser.add_argument("-db_point", "--db_path_point",
                    help="Path to the databases for PointFinder",
                    default=None)
parser.add_argument("-db_point_kma", "--db_path_point_kma",
                    help="Path to the PointFinder databases indexed with KMA. \
                          Defaults to the 'kma_indexing' directory inside the \
                          given database directory.",
                    default=None)
parser.add_argument("-g", "--specific_gene",
                    nargs='+',
                    help="Specify genes existing in the database to \
                          search for - if none is specified all genes are \
                          included in the search.",
                    default=None)
parser.add_argument("-u", "--unknown_mut",
                    action="store_true",
                    help="Show all mutations found even if in unknown to the\
                          resistance database",
                    default=None)
parser.add_argument("-l_p", "--min_cov_point",
                    help=("Minimum (breadth-of) coverage of Pointfinder within "
                          "the range 0-1. If None is selected, the minimum "
                          "coverage of ResFinder will be used."),
                    type=float,
                    default=None)
parser.add_argument("-t_p", "--threshold_point",
                    help=("Threshold for identity of Pointfinder within the "
                          "range 0-1. If None is selected, the minimum "
                          "coverage of ResFinder will be used."),
                    type=float,
                    default=None)

# Temporary option only available temporary
parser.add_argument("--pickle",
                    action="store_true",
                    help="Create a pickle dump of the Isolate object. \
                          Currently needed in the CGE webserver. Dependency \
                          and this option is being removed.",
                    default=False)

args = parser.parse_args()

# Parse and check all arguments and expected files.
conf = Config(args)

# Initialise result dict
std_result = Result.init_software_result(name="ResFinder",
                                         gitdir=conf.resfinder_root)
if(conf.acquired):
    std_result.init_database("ResFinder", conf.db_path_res)
if(conf.point):
    std_result.init_database("PointFinder", conf.db_path_point_root)

##########################################################################
# ResFinder
##########################################################################

if(conf.acquired is True):

    blast_results = None
    kma_run = None

    # Actually running ResFinder (for acquired resistance)
    acquired_finder = ResFinder(db_conf_file=conf.db_config_file,
                                databases=conf.databases,
                                db_path=conf.db_path_res,
                                notes=conf.db_notes_file,
                                db_path_kma=conf.db_path_res_kma)

    if(conf.inputfasta):
        blast_results = acquired_finder.blast(inputfile=conf.inputfasta,
                                              out_path=conf.outPath_res_blast,
                                              min_cov=conf.rf_gene_cov,
                                              threshold=conf.rf_gene_id,
                                              blast=conf.blast,
                                              allowed_overlap=conf.rf_overlap)

        # DEPRECATED
        # TODO: make a write method that depends on the json output
        acquired_finder.write_results(out_path=conf.outputPath,
                                      result=blast_results,
                                      res_type=ResFinder.TYPE_BLAST)

        ResFinderResultHandler.standardize_results(std_result,
                                                   blast_results.results,
                                                   "ResFinder")

    else:
        kma_run = acquired_finder.kma(inputfile_1=conf.inputfastq_1,
                                      inputfile_2=conf.inputfastq_2,
                                      out_path=conf.outPath_res_kma,
                                      db_path_kma=conf.db_path_res_kma,
                                      databases=acquired_finder.databases,
                                      min_cov=conf.rf_gene_cov,
                                      threshold=conf.rf_gene_id,
                                      kma_path=conf.kma,
                                      sample_name="",
                                      kma_cge=True,
                                      kma_apm="p",
                                      kma_1t1=True)

        # DEPRECATED
        # TODO: make a write method that depends on the json output
        acquired_finder.write_results(out_path=conf.outputPath,
                                      result=kma_run.results,
                                      res_type=ResFinder.TYPE_KMA)

        ResFinderResultHandler.standardize_results(std_result,
                                                   kma_run.results,
                                                   "ResFinder")

##########################################################################
# PointFinder
##########################################################################

if(conf.point):

    finder = PointFinder(db_path=conf.db_path_point, species=conf.species_dir,
                         gene_list=conf.specific_gene)

    if(conf.inputfasta):

        method = PointFinder.TYPE_BLAST

        blast_run = finder.blast(inputfile=conf.inputfasta,
                                 out_path=conf.outPath_point_blast,
                                 min_cov=0.01,  # Sorts on coverage later
                                 threshold=conf.pf_gene_id,
                                 blast=conf.blast,
                                 cut_off=False)
        results = blast_run.results

    else:

        method = PointFinder.TYPE_KMA

        kma_run = finder.kma(inputfile_1=conf.inputfastq_1,
                             inputfile_2=conf.inputfastq_2,
                             out_path=conf.outPath_point_kma,
                             db_path_kma=conf.db_path_point,
                             databases=[conf.species_dir],
                             min_cov=0.01,  # Sorts on coverage later
                             threshold=conf.pf_gene_id,
                             kma_path=conf.kma,
                             sample_name=conf.sample_name,
                             kma_cge=True,
                             kma_apm="p",
                             kma_1t1=True)

        results = kma_run.results

    if(conf.specific_gene):
        results = PointFinder.discard_unwanted_results(
            results=results, wanted=conf.specific_gene)

    if(method == PointFinder.TYPE_BLAST):
        results_pnt = finder.find_best_seqs(results, conf.rf_gene_cov)
    else:
        results_pnt = results[finder.species]
        if(results_pnt == "No hit found"):
            results_pnt = {}
        else:
            results_pnt["excluded"] = results["excluded"]

    # DEPRECATED
    # TODO: make a write method that depends on the json output
    finder.write_results(out_path=conf.outputPath, result=results,
                         res_type=method, unknown_flag=conf.unknown_mut,
                         min_cov=conf.pf_gene_cov, perc_iden=conf.pf_gene_id)

    PointFinderResultHandler.standardize_results(std_result,
                                                 results_pnt,
                                                 "PointFinder")

##########################################################################
# Phenotype to genotype
##########################################################################

# Load genotype to phenotype database
res_pheno_db = PhenoDB(
    abclassdef_file=conf.abclassdef_file, acquired_file=conf.phenotype_file,
    point_file=conf.point_file)

# Isolate object store results
isolate = Isolate(name=conf.sample_name)

if(conf.acquired):
    isolate.load_finder_results(std_table=std_result,
                                phenodb=res_pheno_db,
                                type="seq_regions")
if(conf.point):
    isolate.load_finder_results(std_table=std_result,
                                phenodb=res_pheno_db,
                                type="seq_variations")

isolate.calc_res_profile(res_pheno_db)
ResFinderResultHandler.load_res_profile(std_result, isolate,
                                        conf.amr_abbreviations)

std_result_file = "{}/std_format.json".format(conf.outputPath)

with open(std_result_file, 'w') as fh:
    fh.write(std_result.json_dumps())

# Create and write the downloadable tab file
pheno_profile_str = isolate.profile_to_str_table(header=True)

# TODO: REMOVE THE NEED FOR THE PICKLED FILE
if(conf.pickle):
    isolate_pickle = open("{}/isolate.p".format(conf.outputPath), "wb")
    pickle.dump(isolate, isolate_pickle, protocol=2)

pheno_table_file = "{}/pheno_table.txt".format(conf.outputPath)
with open(pheno_table_file, 'w') as fh:
    fh.write(pheno_profile_str)

if(conf.species is not None):
    # Apply AMR panel
    input_amr_panels = "{}/phenotype_panels.txt".format(conf.db_path_res)
    res_sum_table = ResSumTable(pheno_profile_str)
    res_sum_table.load_amr_panels(input_amr_panels)

    try:
        panel_profile_str = res_sum_table.get_amr_panel_str(
            panel_name_raw=conf.species, header=True)
    # If specified species does not have an associated panel, just ignore it
    # and exit.
    except PanelNameError:
        eprint("Warning: No panel was detected for the species: {}"
               .format(conf.species))
        sys.exit()

    amr_panel_filename = conf.species.replace(" ", "_")

    panel_tabel_file = ("{tablename}_{species}.txt"
                        .format(tablename=pheno_table_file[:-4],
                                species=amr_panel_filename))

    with open(panel_tabel_file, "w") as fh:
        fh.write(panel_profile_str)

sys.exit()
