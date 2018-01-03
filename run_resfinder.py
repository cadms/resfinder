#!/home/data1/tools/bin/anaconda/bin/python
import sys
import os
import subprocess
from argparse import ArgumentParser

#  Modules used to create the extended ResFinder output (phenotype output)
from phenotype2genotype.isolate import Isolate
from phenotype2genotype.res_profile import PhenoDB
from phenotype2genotype.res_sumtable import ResSumTable

# TODO list:
# TODO: Python path
# TODO: Add input check

python = "/home/data1/tools/bin/anaconda/bin/python"


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
parser.add_argument("-i", "--inputfile",
                    dest="inputfile",
                    help="Input file",
                    default='')
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
                    default='')
parser.add_argument("-d", "--databases",
                    dest="databases",
                    help="Databases chosen to search in - if non is specified\
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

# TODO: Add input check
scripts = args.scripts
inputfile = args.inputfile
# path = args.out_path
blast = args.blast_path
species = args.species

# Check output directory
args.out_path = os.path.abspath(args.out_path)
if(not os.path.isfile(args.out_path)):
    print("Output directory not found:", args.config_queue)
    quit(1)

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

   out_res = args.out_path + "/resfinder_out"
   os.makedirs(out_res, exist_ok=True)

   db_path_res = args.db_path_res

   # Run ResFinder
   # TODO: Python path
   cmd = ("%s %s -i "
          "%s -o %s -b %s -p %s -d %s -l %f -t %f" % (python, script_resfinder,
                                                      inputfile,
                                                      out_res, blast,
                                                      db_path_res,
                                                      databases, min_cov,
                                                      threshold))
   process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
   out, err = process.communicate()

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
