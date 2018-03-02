#!/home/data1/tools/bin/anaconda/bin/python
import sys
import os
import subprocess
from argparse import ArgumentParser

from cge.resfinder import ResFinder

#  Modules used to create the extended ResFinder output (phenotype output)
from phenotype2genotype.isolate import Isolate
from phenotype2genotype.res_profile import PhenoDB
from phenotype2genotype.res_sumtable import ResSumTable

# TODO list:
# TODO: Add input check

# python = "/home/data1/tools/bin/anaconda/bin/python"
python = "/services/tools/anaconda3/4.0.0/bin/python3"


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
                    default="cge/kma/kma")
parser.add_argument("-s", "--species",
                    dest="species",
                    help="Species in the sample")
parser.add_argument("-u", "--unknown_mut",
                    dest="unknown_mutations",
                    action="store_true",
                    help="Show all mutations found even if in unknown to the\
                          resistance database",
                    default=False)

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
parser.add_argument("-l", "--min_cov",
                    dest="min_cov",
                    help="Minimum coverage",
                    default=0.60)
parser.add_argument("-t", "--threshold",
                    dest="threshold",
                    help="Blast threshold for identity",
                    default=0.90)
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
                    default='')

args = parser.parse_args()


##########################################################################
# MAIN
##########################################################################

# TODO: Add more input checks?
scripts = args.scripts
if(args.inputfastq):
   inputfastq_1 = args.inputfastq[0]
   if(len(args.inputfastq) == 2):
      inputfastq_2 = args.inputfastq[1]
   else:
      inputfastq_2 = None

blast = args.blast_path
kma = args.kma_path
species = args.species

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

if args.acquired is False and args.pont is False:
   sys.exit("Please specify to look for acquired resistance genes, "
            "chromosomal mutaitons or both!\n")

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
   if db_path_res is None:
         sys.exit("Input Error: No database directory was provided!\n")
   elif not os.path.exists(db_path_res):
      sys.exit("Input Error: The specified database directory does not "
               "exist!\n")
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
                               notes=notes_path, db_path_kma=args.db_path_kma)

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
                                        out_path=out_res_kma, min_cov=min_cov,
                                        kma_path=kma)

      acquired_finder.write_results(out_path=out_res_kma, result=kma_results,
                                    res_type=ResFinder.TYPE_BLAST)

if args.point is True:
   db_path_point = args.db_path_point

   out_point = args.out_path + "/pointfinder_out"
   os.makedirs(out_point, exist_ok=True)

   # Run PointFinder
   # TODO: Python path
   cmd = ("%s %s "
          "-i %s -o %s -s %s -p %s -b %s" % (python, script_pointfinder,
                                             inputfile, out_point, species,
                                             db_path_point, blast))
   # Add -u if unknown mutations should be included in the output
   if args.unknown_mutations is True:
      cmd = cmd + ' -u'

   print(cmd)
   process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
   out, err = process.communicate()
