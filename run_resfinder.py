#!/home/data1/tools/bin/anaconda/bin/python
import sys
import os
import subprocess
from argparse import ArgumentParser

from cge.resfinder import ResFinder
from cge.pointfinder import PointFinder

#  Modules used to create the extended ResFinder output (phenotype output)
from cge.phenotype2genotype.isolate import Isolate
from cge.phenotype2genotype.res_profile import PhenoDB
from cge.phenotype2genotype.res_sumtable import ResSumTable

# TODO list:
# TODO: Add input data check


# ########################################################################### #
# #########                         FUNCTIONS                       ######### #
# ########################################################################### #

def create_tab_acquired(isolate, phenodb):
   """ Alternative method to create the downloadeable tabbed result file. This
       method will include the additional information from the phenotype
       database.
   """
   output_str = ("Resistance gene\tIdentity\tAlignment Length/Gene Length\t"
                 "Position in reference\tContig\tPosition in contig\tPhenotype"
                 "\tClass\tPMID\tAccession no.\tNotes\n")

   for unique_id in isolate:
      for feature in isolate[unique_id]:

         # Extract phenotypes
         phenotype_out_list = []
         phenotype = phenodb[feature.unique_id]

         # Append stars to phenotypes that are suggested by the curators and
         # not published
         for antibiotic in phenotype.phenotype:
            if(antibiotic in phenotype.sug_phenotype):
               antibiotic = antibiotic + "*"
            phenotype_out_list.append(antibiotic)

         phenotype_out_str = ",".join(phenotype_out_list)

         output_str += (feature.hit.name + "\t"
                        + str(feature.hit.identity) + "\t"
                        + str(feature.hit.match_length)
                        + "/" + str(feature.hit.ref_length) + "\t"
                        + str(feature.hit.start_ref)
                        + ".." + str(feature.hit.end_ref) + "\t"
                        + feature.seq_region + "\t"
                        + str(feature.start)
                        + ".." + str(feature.end) + "\t"
                        + phenotype_out_str + "\t"
                        + ",".join(phenotype.ab_class) + "\t"
                        + ",".join(phenotype.pmid) + "\t"
                        + feature.hit.acc + "\t"
                        + phenotype.notes + "\n")

   # Find AMR classes with no hits
   no_class_hits = []
   for ab_class in phenodb.antibiotics:
      if(ab_class not in isolate.resprofile.resistance_classes):
         no_class_hits.append(ab_class)

   if(no_class_hits):
      output_str += ("\nNo hits found in the classes: "
                     + ",".join(no_class_hits) + "\n")

   return output_str


# TODO: Add fix species choice
species_transl = {"C. jejuni": "campylobacter",
                  "C. jejuni": "C. jejuni",
                  "E. coli": "e.coli",
                  "E. coli": "E. coli",
                  }

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()

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

parser.add_argument("-scripts", "--scrtips",
                    dest="scripts",
                    help="Path to ResFinder and PointFinder scritps. Defaults\
                          to the directory of run_resfinder.py",
                    default=None)
parser.add_argument("-o", "--outputPath",
                    dest="out_path",
                    help="Path to blast output",
                    default='')
parser.add_argument("-b", "--blastPath",
                    dest="blast_path",
                    help="Path to blast",
                    default='blastn')
parser.add_argument("-k", "--kmaPath",
                    dest="kma_path",
                    help="Path to KMA",
                    default=None)
parser.add_argument("-s", "--species",
                    dest="species",
                    help="Species in the sample")
parser.add_argument("-l", "--min_cov",
                    dest="min_cov",
                    help="Minimum coverage",
                    type=float,
                    default=0.60)
parser.add_argument("-t", "--threshold",
                    dest="threshold",
                    help="Threshold for identity",
                    type=float,
                    default=0.80)

# Acquired resistance options
parser.add_argument("-db_res", "--databasePath_res",
                    dest="db_path_res",
                    help="Path to the databases for ResFinder",
                    default=None)
parser.add_argument("-db_res_kma", "--databasePath_res_kma",
                    dest="db_path_kma",
                    help="Path to the ResFinder databases indexed with KMA. \
                          Defaults to the 'kma_indexing' directory inside the \
                          given database directory.",
                    default=None)
parser.add_argument("-d", "--databases",
                    dest="databases",
                    help="Databases chosen to search in - if none is specified\
                          all is used",
                    default=None)
parser.add_argument("-acq", "--acquired",
                    action="store_true",
                    dest="acquired",
                    help="Run resfinder for acquired resistance genes",
                    default=False)

# Point resistance option
parser.add_argument("-c", "--point",
                    action="store_true",
                    dest="point",
                    help="Run pointfinder for chromosomal mutations",
                    default=False)
parser.add_argument("-db_point", "--databasePath_point",
                    dest="db_path_point",
                    help="Path to the databases for PointFinder",
                    default=None)
parser.add_argument("-g",
                    dest="specific_gene",
                    nargs='+',
                    help="Specify genes existing in the database to \
                          search for - if none is specified all genes are \
                          included in the search.",
                    default=None)
parser.add_argument("-u", "--unknown_mut",
                    dest="unknown_mutations",
                    action="store_true",
                    help="Show all mutations found even if in unknown to the\
                          resistance database",
                    default=False)

# Phemotype2genotype options pheno_db_path
parser.add_argument("-db_pheno", "--databasePath_pheno",
                    dest="pheno_db_path",
                    help="Path to phenotype database.",
                    default=None)
parser.add_argument("--amr_panel",
                    help="Name of AMR panel to apply. If None is specified all\
                          possible phenotypes are reported.",
                    default=None)

args = parser.parse_args()

# Create a "sample" name
if(args.inputfasta):
   sample_name = os.path.basename(args.inputfasta)
   method = PointFinder.TYPE_BLAST
else:
   sample_name = os.path.basename(args.inputfastq[0])
   method = PointFinder.TYPE_KMA

# TODO: Add input data check
scripts = args.scripts
if(args.inputfastq):
   inputfastq_1 = args.inputfastq[0]
   if(len(args.inputfastq) == 2):
      inputfastq_2 = args.inputfastq[1]
   else:
      inputfastq_2 = None

blast = args.blast_path
species = args.species

# Check KMA path cge/kma/kma
if(args.inputfastq and args.kma_path is None):
   kma = (os.path.dirname(
       os.path.realpath(__file__)) + "/cge/kma/kma")
   kma = os.path.abspath(kma)

# Check output directory
args.out_path = os.path.abspath(args.out_path)
os.makedirs(args.out_path, exist_ok=True)

# Check script directory.
if(not args.scripts):
    script_resfinder = os.path.dirname(
        os.path.realpath(__file__)) + "/ResFinder.py"
    script_pointfinder = os.path.dirname(
        os.path.realpath(__file__)) + "/PointFinder.py"
else:
   script_resfinder = scritps + "/ResFinder.py"
   script_pointfinder = scritps + "/PointFinder.py"

if args.acquired is False and args.point is False:
   sys.exit("Please specify to look for acquired resistance genes, "
            "chromosomal mutaitons or both!\n")

# Check ResFinder KMA database
if(args.db_path_kma is None and args.acquired):
   db_path_kma = (os.path.dirname(
       os.path.realpath(__file__)) + "/database/kma_indexing/")
   db_path_kma = os.path.abspath(db_path_kma)

# Check Poinfinder database
if(args.db_path_point is None and args.point):
   db_path_point = (os.path.dirname(
       os.path.realpath(__file__)) + "/database_pointfinder/"
       + args.species)
   db_path_point = os.path.abspath(db_path_point)

# Check phenotype database
if(args.pheno_db_path is None):
   pheno_db_path = os.path.dirname(
       os.path.realpath(__file__)) + "/database_phenotype"
   pheno_db_path = os.path.abspath(pheno_db_path)

if not os.path.exists(pheno_db_path):
   sys.exit("Input Error: The specified phenotype database directory does not "
            "exist!\nProvided path: " + str(pheno_db_path))

##########################################################################
# ResFinder
##########################################################################

if args.acquired is True:

   databases = args.databases
   min_cov = float(args.min_cov)
   threshold = float(args.threshold)

   if(args.inputfasta):
      out_res_blast = args.out_path + "/resfinder_blast"
      os.makedirs(out_res_blast, exist_ok=True)
   if(args.inputfastq):
      out_res_kma = args.out_path + "/resfinder_kma"
      os.makedirs(out_res_kma, exist_ok=True)

   db_path_res = args.db_path_res

   # Check if valid database is provided
   if(db_path_res is None):
      db_path_res = os.path.dirname(os.path.realpath(__file__)) + "/database"

   if not os.path.exists(db_path_res):
      sys.exit("Input Error: The specified database directory does not "
               "exist!\nProvided path: " + str(db_path_res))
   else:
      # Check existence of config file
      db_config_file = '%s/config' % (db_path_res)
      if not os.path.exists(db_config_file):
         sys.exit("Input Error: The database config file could not be found!")

   # Check existence of notes file
   notes_path = "%s/notes.txt" % (db_path_res)
   if not os.path.exists(notes_path):
      sys.exit('Input Error: notes.txt not found! (%s)' % (notes_path))

   # Actually running ResFinder (for acquired resistance)
   acquired_finder = ResFinder(db_conf_file=db_config_file,
                               databases=args.databases, db_path=db_path_res,
                               notes=notes_path, db_path_kma=db_path_kma)

   blast_results = None
   kma_results = None

   if(args.inputfasta):
      blast_results = acquired_finder.blast(inputfile=args.inputfasta,
                                            out_path=out_res_blast,
                                            min_cov=min_cov,
                                            threshold=threshold,
                                            blast=blast)

      acquired_finder.write_results(out_path=out_res_blast,
                                    result=blast_results,
                                    res_type=ResFinder.TYPE_BLAST)

   if(args.inputfastq):
      kma_results = acquired_finder.kma(inputfile_1=inputfastq_1,
                                        inputfile_2=inputfastq_2,
                                        out_path=out_res_kma,
                                        db_path_kma=db_path_kma,
                                        databases=acquired_finder.databases,
                                        min_cov=min_cov,
                                        threshold=args.threshold,
                                        kma_path=kma,
                                        sample_name="",
                                        kma_mrs=0.5, kma_gapopen=-3,
                                        kma_gapextend=-1, kma_penalty=-2,
                                        kma_reward=1,
                                        kma_pm="p",
                                        kma_fpm="p")

      acquired_finder.write_results(out_path=out_res_kma, result=kma_results,
                                    res_type=ResFinder.TYPE_KMA)

##########################################################################
# PointFinder
##########################################################################

if args.point is True:

   if(args.inputfasta):
      out_point = os.path.abspath(args.out_path + "/pointfinder_blast")
      os.makedirs(out_point, exist_ok=True)
   if(args.inputfastq):
      out_point = os.path.abspath(args.out_path + "/pointfinder_kma")
      os.makedirs(out_point, exist_ok=True)

   finder = PointFinder(db_path=db_path_point, species=args.species,
                        gene_list=args.specific_gene)

   if(args.inputfasta):
      blast_run = finder.blast(inputfile=args.inputfasta,
                               out_path=out_point,
                               min_cov=args.min_cov,
                               threshold=args.threshold,
                               blast=blast,
                               cut_off=False)
      results = blast_run.results

   # Note: ResFinder is able to do a fasta and a fastq call, hence its
   #       two if statements. PointFinder can only handle eiter fasta
   #       or fastq, hence the if-else statement.
   else:

      method = PointFinder.TYPE_KMA

      results = finder.kma(inputfile_1=inputfastq_1,
                           inputfile_2=inputfastq_2,
                           out_path=out_point,
                           db_path_kma=db_path_point,
                           databases=[args.species],
                           min_cov=args.min_cov,
                           threshold=args.threshold,
                           kma_path=kma,
                           sample_name="",
                           kma_mrs=0.5, kma_gapopen=-5, kma_gapextend=-2,
                           kma_penalty=-3, kma_reward=1, kma_pm="p",
                           kma_fpm="p")

   if(args.specific_gene):
      results = PointFinder.discard_unwanted_results(results=results,
                                                     wanted=args.specific_gene)

   finder.write_results(out_path=args.out_path, result=results,
                        res_type=method, unknown_flag=args.unknown_mutations)

##########################################################################
# Phenotype to genotype
##########################################################################

# Load genotype to phenotype database
res_pheno_db = PhenoDB(acquired_file=pheno_db_path + "/acquired_db.txt",
                       point_file=db_path_point + "/resistens-overview.txt")

# Isolate object stores results
isolate = Isolate(name=sample_name)

if(args.acquired):
   if(method == PointFinder.TYPE_BLAST):
      isolate.load_resfinder_tab(out_res_blast + "/results_table.txt")
   else:
      isolate.load_resfinder_tab(out_res_kma + "/results_table.txt")
if(args.point):
   isolate.load_pointfinder_tab(args.out_path + "/PointFinder_results.txt")

isolate.calc_res_profile(res_pheno_db)

# Create and write the downloadable tab file
pheno_profile_str = isolate.profile_to_str_table(header=True)

pheno_table_file = args.out_path + '/pheno_table.txt'
with open(pheno_table_file, 'w') as fh:
   fh.write(pheno_profile_str)

if(args.amr_panel is not None):
   # Apply AMR panel
   input_amr_panels = pheno_db_path + "/amr_panels.txt"
   res_sum_table = ResSumTable(pheno_profile_str)
   res_sum_table.load_amr_panels(input_amr_panels)
   panel_profile_str = res_sum_table.get_amr_panel_str(
       panel_name=args.amr_panel, header=True)

   amr_panel_filename = args.amr_panel.replace(".", "")
   amr_panel_filename = amr_panel_filename.replace(" ", "")

   panel_tabel_file = pheno_table_file[:-4] + "_" + amr_panel_filename + ".txt"
   with open(panel_tabel_file, "w") as fh:
      fh.write(panel_profile_str)

sys.exit()
