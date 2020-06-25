#!/usr/bin/env python3
from __future__ import division
from argparse import ArgumentParser
from cgecore.blaster import Blaster
from cgecore.cgefinder import CGEFinder
from distutils.spawn import find_executable
from tabulate import tabulate
import collections
import gzip
import json
import os
import pprint
import random
import re
import subprocess
import sys
import time


class ResFinder(CGEFinder):

   def __init__(self, db_conf_file, notes, db_path, db_path_kma=None, databases=None):
      self.db_path = db_path

      if(db_path_kma is None):
         self.db_path_kma = db_path + "/kma_indexing"
      else:
         self.db_path_kma = db_path_kma

      self.configured_dbs = dict()
      self.description_dbs = dict()
      self.kma_db_files = None
      self.load_db_config(db_conf_file=db_conf_file)

      self.databases = []
      self.load_databases(databases=databases)

      self.phenos = dict()
      self.load_notes(notes=notes)
      self.blast_results = None


   def blast(self, inputfile, out_path, min_cov=0.9, threshold=0.6, blast="blastn", allowed_overlap=0):
      blast_run = Blaster(inputfile=inputfile, databases=self.databases,
                          db_path=self.db_path, out_path=out_path,
                          min_cov=min_cov, threshold=threshold, blast=blast,
                          allowed_overlap=allowed_overlap)
      self.blast_results = blast_run.results
      return blast_run

   def create_results(self, results_method, outdir, json_out=False):

      results = results_method.results
      query_aligns = results_method.gene_align_query
      homo_aligns = results_method.gene_align_homo
      sbjct_aligns = results_method.gene_align_sbjct
      json_results = dict()

      titles = dict()
      rows = dict()
      headers = dict()
      txt_file_seq_text = dict()
      split_print = collections.defaultdict(list)

      hits = []
      for db in results:
          contig_res = {}
          if db == 'excluded':
              continue
          db_name = str(db).capitalize()
          if db_name not in json_results:
              json_results[db_name] = {}
          if db not in json_results[db_name]:
              json_results[db_name][db] = {}
          if results[db] == "No hit found":
              json_results[db_name][db] = "No hit found"
          else:
              for contig_id, hit in results[db].items():
                  identity = float(hit["perc_ident"])
                  coverage = float(hit["perc_coverage"])

                   # Skip hits below coverage
                  if coverage < (min_cov * 100) or identity < (threshold * 100):
                      continue

                  bit_score = identity * coverage
                  if contig_id not in contig_res:
                      contig_res[contig_id] = []
                  contig_res[contig_id].append([hit["query_start"], hit["query_end"],
                                               bit_score, hit])

          if not contig_res:
              json_results[db_name][db] = "No hit found"

          # Check for overlapping hits, only report the best
          for contig_id, hit_lsts in contig_res.items():

              hit_lsts.sort(key=lambda x: x[0])
              hits = [hit[3] for hit in hit_lsts]

              for hit in hits:

                  header = hit["sbjct_header"]
                  tmp = header.split("_")
                  try:
                      gene = tmp[0]
                      note = tmp[1]
                      acc = "_".join(tmp[2:])
                  except IndexError:
                      gene = ":".join(tmp)
                      note = ""
                      acc = ""
                  try:
                      variant = tmp[3]
                  except IndexError:
                      variant = ""

                  identity = hit["perc_ident"]
                  coverage = hit["perc_coverage"]
                  sbj_length = hit["sbjct_length"]
                  HSP = hit["HSP_length"]
                  positions_contig = "%s..%s" % (hit["query_start"],
                                                 hit["query_end"])
                  positions_ref = "%s..%s" % (hit["sbjct_start"], hit["sbjct_end"])
                  contig_name = hit["contig_name"]
                  pheno = self.phenos.get(gene, "Warning: gene is missing from "
                                       "Notes file. Please inform curator.")
                  pheno = pheno.strip()

                  # Write JSON results dict
                  json_results[db_name][db].update({contig_id: {}})
                  json_results[db_name][db][contig_id] = {
                      "resistance_gene": gene,
                      "identity": round(identity, 2),
                      "HSP_length": HSP,
                      "template_length": sbj_length,
                      "position_in_ref": positions_ref,
                      "contig_name": contig_name,
                      "positions_in_contig": positions_contig,
                      "note": note,
                      "accession": acc,
                      "predicted_phenotype": pheno,
                      "coverage": round(coverage, 2),
                      "hit_id": contig_id,
                      "subject_header": header}

      # Get run info for JSON file
      service = os.path.basename(__file__).replace(".py", "")
      date = time.strftime("%d.%m.%Y")
      time_ = time.strftime("%H:%M:%S")

      # Make JSON output file
      data = {service: {}}

      userinput = {"filename(s)": args.inputfile,
                   "method": method,
                   "file_format": file_format}
      run_info = {"date": date, "time": time_}

      data[service]["user_input"] = userinput
      data[service]["run_info"] = run_info
      data[service]["results"] = json_results

      if(json_out):
          print(json.dumps(data))
      else:
          pprint.pprint(data)

      # Save json output
      result_file = "{}/data_resfinder.json".format(tmp_dir)
      with open(result_file, "w") as outfile:
          json.dump(data, outfile)

      # Getting and writing out the results
      header = ["Resistance gene", "Identity", "Query / Template length",
                "Contig", "Position in contig", "Predicted phenotype",
                "Accession number"]

      if args.extended_output:
          # Define extented output
          table_filename = "{}/results_tab.tsv".format(outdir)
          query_filename = "{}/Hit_in_genome_seq.fsa".format(outdir)
          sbjct_filename = "{}/Resistance_genes.fsa".format(outdir)
          result_filename = "{}/results.txt".format(outdir)
          table_file = open(table_filename, "w")
          query_file = open(query_filename, "w")
          sbjct_file = open(sbjct_filename, "w")
          result_file = open(result_filename, "w")
          # Make results file
          result_file.write("{} Results\n\nAntibiotic(s): {}\n\n"
                            .format(service, ", ".join(self.configured_dbs)))
          # Write tsv table
          ## TODO##
          rows = [["Database"] + header]
          for species, dbs_info in sorted(json_results.items()):
              for db_name, db_hits in dbs_info.items():
                  result_file.write("*" * len("\t".join(header)) + "\n")
                  result_file.write(db_name.capitalize() + "\n")
                  db_rows = []

                  # Check it hits are found
                  if isinstance(db_hits, str):
                      content = [''] * len(header)
                      content[int(len(header) / 2)] = db_hits
                      result_file.write(ResFinder.text_table(header, [content]) + "\n")
                      continue

                  for gene_id, gene_info in sorted(
                          db_hits.items(),
                          key=lambda x: (x[1]['resistance_gene'],
                                         x[1]['accession'])):

                      res_gene = gene_info["resistance_gene"]
                      identity = str(gene_info["identity"])
                      coverage = str(gene_info["coverage"])

                      template_HSP = (
                          "{hsp_len} / {template_len}".
                          format(hsp_len=gene_info["HSP_length"],
                                 template_len=gene_info["template_length"]))

                      position_in_ref = gene_info["position_in_ref"]
                      position_in_contig = gene_info["positions_in_contig"]
                      predicted_phenotype = gene_info["predicted_phenotype"]
                      acc = gene_info["accession"]
                      contig_name = gene_info["contig_name"]

                      # Add rows to result tables
                      db_rows.append([res_gene, identity, template_HSP, contig_name,
                                      position_in_contig, predicted_phenotype, acc])
                      rows.append([db_name, res_gene, identity, template_HSP,
                                   contig_name, position_in_contig, predicted_phenotype,
                                   acc])

                      # Write query fasta output
                      hit_name = gene_info["hit_id"]
                      query_seq = query_aligns[db_name][hit_name]
                      sbjct_seq = sbjct_aligns[db_name][hit_name]

                      if coverage == 100 and identity == 100:
                          match = "PERFECT MATCH"
                      else:
                          match = "WARNING"
                      qry_header = (">{}:{} ID:{}% COV:{}% Best_match:{}\n"
                                    .format(res_gene, match, identity, coverage,
                                            gene_id))
                      query_file.write(qry_header)
                      for i in range(0, len(query_seq), 60):
                          query_file.write(query_seq[i:i + 60] + "\n")

                      # Write template fasta output
                      sbj_header = ">{}\n".format(gene_id)
                      sbjct_file.write(sbj_header)
                      for i in range(0, len(sbjct_seq), 60):
                          sbjct_file.write(sbjct_seq[i:i + 60] + "\n")

                  # Write db results tables in results file and table file
                  result_file.write(ResFinder.text_table(header, db_rows) + "\n")

              result_file.write("\n")

          for row in rows:
              table_file.write("\t".join(row) + "\n")

          # Write allignment output
          result_file.write("\n\nExtended Output:\n\n")
          self.make_aln(result_file, json_results, query_aligns, homo_aligns,
                   sbjct_aligns)

          # Close all files

          query_file.close()
          sbjct_file.close()
          table_file.close()
          result_file.close()
      if args.quiet:
          f.close()

   def make_aln(self,file_handle, json_data, query_aligns, homol_aligns, sbjct_aligns):
       for dbs_info in json_data.values():
           for db_name, db_info in dbs_info.items():
               if isinstance(db_info, str):
                   continue

               for gene_id, gene_info in sorted(
                       db_info.items(),
                       key=lambda x: (x[1]['resistance_gene'],
                                      x[1]['accession'])):

                   seq_name = ("{gene}_{acc}"
                               .format(gene=gene_info["resistance_gene"],
                                       acc=gene_info["accession"]))
                   hit_name = gene_info["hit_id"]

                   seqs = ["", "", ""]
                   seqs[0] = sbjct_aligns[db_name][hit_name]
                   seqs[1] = homol_aligns[db_name][hit_name]
                   seqs[2] = query_aligns[db_name][hit_name]

                   self.write_align(seqs, seq_name, file_handle)

   def write_align(self,seq, seq_name, file_handle):
       file_handle.write("# {}".format(seq_name) + "\n")
       sbjct_seq = seq[0]
       homol_seq = seq[1]
       query_seq = seq[2]
       for i in range(0, len(sbjct_seq), 60):
           file_handle.write("%-10s\t%s\n" % ("template:", sbjct_seq[i:i + 60]))
           file_handle.write("%-10s\t%s\n" % ("", homol_seq[i:i + 60]))
           file_handle.write("%-10s\t%s\n\n" % ("query:", query_seq[i:i + 60]))


   @staticmethod
   def text_table(headers, rows, empty_replace='-'):
       ''' Create text table

       USAGE:
           >>> from tabulate import tabulate
           >>> headers = ['A','B']
           >>> rows = [[1,2],[3,4]]
           >>> print(text_table(headers, rows))
           **********
             A     B
           **********
             1     2
             3     4
           ==========
       '''
       # Replace empty cells with placeholder
       rows = map(lambda row: map(lambda x: x if x else empty_replace, row), rows)
       # Create table
       table = tabulate(rows, headers, tablefmt='simple').split('\n')
       # Prepare title injection
       width = len(table[0])
       # Switch horisontal line
       table[1] = '*' * (width + 2)
       # Update table with title
       table = (("%s\n" * 3)
                % ('*' * (width + 2), '\n'.join(table), '=' * (width + 2)))
       return table

   def load_notes(self, notes):
      with open(notes, 'r') as f:
         for line in f:
            line = line.strip()
            if line.startswith("#"):
               continue
            else:
               tmp = line.split(":")

               self.phenos[tmp[0]] = "%s %s" % (tmp[1], tmp[2])

               if(tmp[2].startswith("Alternate name; ")):
                  self.phenos[tmp[2][16:]] = "%s %s" % (tmp[1], tmp[2])

   def load_databases(self, databases):
      # Check if databases and config file are correct/correponds
      if databases == '':
            sys.exit("Input Error: No database was specified!\n")
      elif databases is None:
         # Choose all available databases from the config file
         self.databases = self.configured_dbs.keys()
      else:
         # Handle multiple databases
         databases = databases.split(',')
         # Check that the ResFinder DBs are valid
         for db_prefix in databases:
            if db_prefix in self.configured_dbs:
               self.databases.append(db_prefix)
            else:
               sys.exit("Input Error: Provided database was not "
                        "recognised! (%s)\n" % db_prefix)

   def load_db_config(self, db_conf_file):
      extensions = []
      with open(db_conf_file) as f:
         for line in f:
            line = line.strip()

            if not line:
               continue

            if line[0] == '#':
               if 'extensions:' in line:
                  extensions = [s.strip() for s
                                in line.split('extensions:')[-1].split(',')]
               continue

            tmp = line.split('\t')
            if len(tmp) != 3:
               sys.exit(("Input Error: Invalid line in the database"
                         " config file!\nA proper entry requires 3 tab "
                         "separated columns!\n%s") % (line))

            db_prefix = tmp[0].strip()
            name = tmp[1].split('#')[0].strip()
            description = tmp[2]

            # Check if all db files are present
            for ext in extensions:
               db = "%s/%s.%s" % (self.db_path, db_prefix, ext)
               if not os.path.exists(db):
                  sys.exit(("Input Error: The database file (%s) "
                            "could not be found!") % (db_path))

            if db_prefix not in self.configured_dbs:
               self.configured_dbs[db_prefix] = []
            self.configured_dbs[db_prefix].append(name)
            self.description_dbs[db_prefix] = description

      if len(self.configured_dbs) == 0:
         sys.exit("Input Error: No databases were found in the "
                  "database config file!")

      # Loading paths for KMA databases.
      for drug in self.configured_dbs:
         kma_db = self.db_path_kma + drug
         self.kma_db_files = [kma_db + ".b", kma_db + ".length.b",
                              kma_db + ".name.b", kma_db + ".align.b"]

def is_gzipped(file_path):
    ''' Returns True if file is gzipped and False otherwise.
         The result is inferred from the first two bits in the file read
         from the input path.
         On unix systems this should be: 1f 8b
         Theoretically there could be exceptions to this test but it is
         unlikely and impossible if the input files are otherwise expected
         to be encoded in utf-8.
    '''
    with open(file_path, mode='rb') as fh:
        bit_start = fh.read(2)
    if(bit_start == b'\x1f\x8b'):
        return True
    else:
        return False

def get_file_format(input_files):
    """
    Takes all input files and checks their first character to assess
    the file format. Returns one of the following strings; fasta, fastq,
    other or mixed. fasta and fastq indicates that all input files are
    of the same format, either fasta or fastq. other indiates that all
    files are not fasta nor fastq files. mixed indicates that the inputfiles
    are a mix of different file formats.
    """
    # Open all input files and get the first character
    file_format = []
    invalid_files = []
    for infile in input_files:
        if is_gzipped(infile):
            f = gzip.open(infile, "rb")
            fst_char = f.read(1)
        else:
            f = open(infile, "rb")
            fst_char = f.read(1)
        f.close()
        # Assess the first character
        if fst_char == b"@":
            file_format.append("fastq")
        elif fst_char == b">":
            file_format.append("fasta")
        else:
            invalid_files.append("other")
    if len(set(file_format)) != 1:
        return "mixed"
    return ",".join(set(file_format))


if __name__ == '__main__':

   ##########################################################################
   # PARSE COMMAND LINE OPTIONS
   ##########################################################################

   parser = ArgumentParser()
   parser.add_argument("-i", "--inputfile",
                       dest="inputfile",
                       help="FASTA or FASTQ input files.",
                       nargs="+",
                       required=True)
   parser.add_argument("-o", "--outputPath",
                       dest="out_path",
                       help="Path to blast output.",
                       default='.')
   parser.add_argument("-tmp", "--tmp_dir",
                       help=("Temporary directory for storage of the results from the external software."))
   parser.add_argument("-mp", "--methodPath",
                       dest="method_path",
                       help="Path to method to use (kma or blastn).")

   parser.add_argument("-p", "--databasePath",
                       dest="db_path",
                       help="Path to the databases.",
                       default='')

   parser.add_argument("-d", "--databases",
                       dest="databases",
                       help="Databases chosen to search in - if none are specified all are used.",
                       default=None)
   parser.add_argument("-l", "--min_cov",
                       dest="min_cov",
                       help="Minimum coverage.",
                       default=0.60)
   parser.add_argument("-t", "--threshold",
                       dest="threshold",
                       help="Blast threshold for identity.",
                       default=0.90)
   parser.add_argument("-ao", "--acq_overlap",
                       help="Genes are allowed to overlap this number of nucleotides. Default: 30.",
                       type=int,
                       default=30)
   parser.add_argument("-matrix", "--matrix",
                       help=("Gives the counts all all called bases at each "
                             "position in each mapped template. Columns are: "
                             "reference base, A count, C count, G count, T "
                             "count, N count,- count."),
                       dest="kma_matrix",
                       action='store_true',
                       default=False)
   parser.add_argument("-x", "--extended_output",
                       help=("Give extented output with allignment files, "
                             "template and query hits in fasta and a tab "
                             "seperated file with gene profile results"),
                       action="store_true")
   parser.add_argument("-j", "--json",
                       help=("JSON is program's default output"),
                       action="store_true",
                       default=True)
   parser.add_argument("-q", "--quiet",
                       action="store_true",
                       default=False)
   args = parser.parse_args()

   ##########################################################################
   # MAIN
   ##########################################################################

   if args.quiet:
       f = open('/dev/null', 'w')
       sys.stdout = f

   # Defining varibales

   min_cov = float(args.min_cov)
   threshold = float(args.threshold)
   method_path = args.method_path

   # Check if valid database is provided
   if args.db_path is None:
         sys.exit("Input Error: Missing database directory!\n")
   elif not os.path.exists(args.db_path):
      sys.exit(f"Input Error: Database directory '{args.db_path}' does not exist!\n")
   else:
      # Check existence of config file
      db_config_file = '%s/config' % (args.db_path)
      if not os.path.exists(db_config_file):
         sys.exit("Input Error: The database config file could not be found!")
      # Save path
      db_path = args.db_path

   # Check existence of notes file
   notes_path = "%s/notes.txt" % (args.db_path)
   if not os.path.exists(notes_path):
      sys.exit('Input Error: notes.txt not found! (%s)' % (notes_path))

   # Check if valid input files are provided
   if args.inputfile is None:
       sys.exit("Input Error: No input file was provided!\n")
   elif not os.path.exists(args.inputfile[0]):
       sys.exit("Input Error: Input file does not exist!\n")
   elif len(args.inputfile) > 1:
       if not os.path.exists(args.inputfile[1]):
           sys.exit("Input Error: Input file does not exist!\n")
       infile = args.inputfile
   else:
       infile = args.inputfile

   # Check if valid output directory exist
   if not os.path.exists(args.out_path):
      # sys.exit(f"Input Error: Output dirctory '{args.out_path}' does not exists!\n")
      os.system('mkdir -m 775 ' + args.out_path)
      out_path = args.out_path
   else:
      out_path = args.out_path

   # Add kma_matrix
   if (args.kma_matrix):
      extra_args = "-matrix"
   else:
      extra_args = None

   # Check if valid tmp directory is provided
   if args.tmp_dir:
       if not os.path.exists(args.tmp_dir):
           sys.exit(f"Input Error: Tmp dirctory, '{args.tmp_dir}', does not exist!\n")
       else:
           tmp_dir = os.path.abspath(args.tmp_dir)
   else:
       tmp_dir = out_path

   file_format = get_file_format(infile)

   if file_format == "fastq":
       if not method_path:
           method_path = "kma"
       if find_executable(method_path) is None:
           sys.exit("No valid path to a kma program was provided. Use the -mp "
                    "flag to provide the path.")

    # Check the number of files
       if len(infile) == 1:
           infile_1 = infile[0]
           infile_2 = None
       elif len(infile) == 2:
           infile_1 = infile[0]
           infile_2 = infile[1]
       else:
           sys.exit("Only 2 input file accepted for raw read data,\
                    if data from more runs is avaliable for the same\
                    sample, please concatinate the reads into two files")

       sample_name = os.path.basename(sorted(args.inputfile)[0])
       method = "kma"
       finder = ResFinder(db_conf_file=db_config_file, databases=args.databases,
                          db_path=db_path, db_path_kma=args.db_path,
                          notes=notes_path)
       # if input_fastq2 is None, it is ignored by the kma method.
       kma_run = finder.kma(inputfile_1=infile_1, inputfile_2=infile_2,
                            out_path=tmp_dir, databases=finder.databases,
                            db_path_kma=finder.db_path_kma,
                            min_cov=min_cov, threshold=threshold,
                            kma_path=method_path,
                            kma_add_args=extra_args)

       finder.create_results(results_method=kma_run, outdir=out_path,
                             json_out=args.json)
   elif file_format == "fasta":
       if not method_path:
           method_path = "blastn"
       if find_executable(method_path) is None:
           sys.exit("No valid path to a blastn program was provided. Use the "
                    "-mp flag to provide the path.")

       # Assert that only one fasta file is inputted
       assert len(infile) == 1, "Only one input file accepted for assembled data"
       infile = infile[0]
       method = "blast"
       finder = ResFinder(db_conf_file=db_config_file, databases=args.databases,
                          db_path=db_path, notes=notes_path)

       blast_run = finder.blast(inputfile=infile, out_path=tmp_dir,
                                min_cov=min_cov, threshold=threshold,
                                blast=method_path,
                                allowed_overlap=args.acq_overlap)
       finder.create_results(results_method=blast_run, outdir=out_path,
                             json_out=args.json)
