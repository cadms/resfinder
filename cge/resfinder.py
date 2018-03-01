#!/home/data1/tools/bin/anaconda/bin/python
from __future__ import division
import sys
import os
import time
import random
import re
import subprocess
from argparse import ArgumentParser
from tabulate import tabulate
import collections

from cge.blaster.blaster import Blaster


class ResFinder():

   def __init__(self, db_conf_file, notes, db_path, databases=None):
      """
      """
      self.db_path = db_path
      self.configured_dbs = dict()
      self.load_db_config(db_conf_file=db_conf_file)
      self.databases = []
      self.load_databases(databases=databases)
      self.phenos = dict()
      self.load_notes(notes=notes)
      self.blast_results = None
      self.kma_results = None
      self.results = None

   def KMA(self, inputfile_1, databases, db_path, out_path, sample_name,
           min_cov, mapping_path):
      """
         I expect that there will only be one hit pr gene, but if there are
         more, I assume that the sequence of the hits are the same in the res
         file and the aln file.
      """

      self.kma_results = dict()

      for drug in databases:
         # TODO: kma database, should probalble be included as an argument
         kma_db = db_path + "/kma_indexing/" + drug
         kma_outfile = out_path + "/kma_" + drug + "_" + sample_name
         kma_cmd = ("%s -i %s -t_db %s -SW -o %s -e 1.0" % (mapping_path,
                    inputfile_1, kma_db, kma_outfile))

         print("KMA_cmd: " + kma_cmd)
         # Call KMA
         process = subprocess.Popen(kma_cmd, shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
         out, err = process.communicate()

         self.kma_results[drug] = 'No hit found'

         # Fetch kma output files
         align_filename = kma_outfile + ".aln"
         res_filename = kma_outfile + ".res"

         # Open res file, find coverage and the gene names of genes found
         with open(res_filename, "r") as res_file:
            header = res_file.readline()
            for line in res_file:
               if self.kma_results[drug] == 'No hit found':
                  self.kma_results[drug] = dict()
               data = [data.strip() for data in line.split("\t")]
               gene = data[0]
               # Check if gene one of the user specified genes
   #                if gene not in gene_list:
   #                    continue
               sbjct_len = int(data[3])
               coverage = float(data[4])
               print(gene)
               if gene not in self.kma_results[drug]:
                  hit = gene
               else:
                  hit = gene + "_" + str(len(self.kma_results[drug][gene]) + 1)

               self.kma_results[drug][hit] = dict()
               self.kma_results[drug][hit]['sbjct_length'] = sbjct_len
               self.kma_results[drug][hit]['coverage'] = coverage / 100
               self.kma_results[drug][hit]["sbjct_string"] = []
               self.kma_results[drug][hit]["query_string"] = []
               self.kma_results[drug][hit]["homology"] = []
               self.kma_results[drug][hit]["sbjct_header"] = gene
               self.kma_results[drug][hit]["split_length"] = 'Not given'
               self.kma_results[drug][hit]["perc_ident"] = coverage
               self.kma_results[drug][hit]["query_start"] = 'Not given'
               self.kma_results[drug][hit]["query_end"] = 'Not given'
               self.kma_results[drug][hit]["contig_name"] = 'Not given'

         if self.kma_results[drug] == 'No hit found':
            continue

         # Open align file
         with open(align_filename, "r") as align_file:
            hit_no = dict()
            gene = ""
            for line in align_file:
               if line.startswith("#"):
                  gene = line[1:].strip()

                  if gene not in hit_no:
                     hit_no[gene] = str(1)
                  else:
                     hit_no[gene] += str(int(hit_no[gene]) + 1)

               else:
                  # Check if gene one of the user specified genes
                  if hit_no[gene] == '1':
                     hit = gene
                  else:
                     hit = gene + "_" + hit_no[gene]

                  if hit in self.kma_results[drug]:
                     line_data = line.split("\t")[-1].strip()
                     if line.startswith("template"):
                        self.kma_results[drug][hit]["sbjct_string"] += [line_data]
                     elif line.startswith("query"):
                        self.kma_results[drug][hit]["query_string"] += [line_data]
                     else:
                        self.kma_results[drug][hit]["homology"] += [line_data]
                  else:
                     print(hit + " not in results: ", self.kma_results)

         # concatinate all sequences lists and find subject start and subject
         # end
         seq_start_search_str = re.compile("^-*(\w+)")

         for hit in self.kma_results[drug]:
            self.kma_results[drug][hit]['sbjct_string'] = "".join(
                self.kma_results[drug][hit]['sbjct_string'])
            self.kma_results[drug][hit]['query_string'] = "".join(
                self.kma_results[drug][hit]['query_string'])
            self.kma_results[drug][hit]['homology'] = "".join(
                self.kma_results[drug][hit]['homology'])

            seq_start_object = seq_start_search_str.search(
                self.kma_results[drug][hit]['query_string'])
            sbjct_start = seq_start_object.start() + 1
            self.kma_results[drug][hit]['sbjct_start'] = sbjct_start
            self.kma_results[drug][hit]["sbjct_end"] = (
                self.kma_results[drug][hit]["sbjct_length"] - sbjct_start + 1)

   def write_results(self, out_path):
      """
      """
      if(self.blast_results is None):
         sys.exit("The blast method needs to be called before calling this "
                  "method.")

      with open(out_path + "/results_tab.txt", "w") as fh:
         fh.write(self.blast_results[0])
      with open(out_path + "/results_table.txt", "w") as fh:
         fh.write(self.blast_results[1])
      with open(out_path + "/results.txt", "w") as fh:
         fh.write(self.blast_results[2])
      with open(out_path + "/Resistance_gene_seq.fsa", "w") as fh:
         fh.write(self.blast_results[3])
      with open(out_path + "/Hit_in_genome_seq.fsa", "w") as fh:
         fh.write(self.blast_results[4])

   def blast(self, inputfile, out_path, min_cov=0.9, threshold=0.6,
             blast="blastn"):
      """
      """
      blast_run = Blaster(inputfile=inputfile, databases=self.databases,
                          db_path=self.db_path, out_path=out_path,
                          min_cov=min_cov, threshold=threshold, blast=blast)
      self.blast_results = blast_run.results

      self.results_to_str(query_align=blast_run.gene_align_query,
                          homo_align=blast_run.gene_align_homo,
                          sbjct_align=blast_run.gene_align_sbjct)

      # self.results = (tab_str, table_str, txt_str, ref_str, hit_str)
      self.write_results(out_path=out_path)

   def results_to_str(self, query_align=None, homo_align=None,
                      sbjct_align=None):

      # TODO: Do not use results variable.
      results = self.blast_results

      # Write the header for the tab file
      tab_str = ("Resistance gene\tIdentity\tAlignment Length/Gene Length\t"
                 "Position in reference\tContig\tPosition in contig\t"
                 "Phenotype\tAccession no.\n")

      table_str = ""
      txt_str = ""
      ref_str = ""
      hit_str = ""

      # Getting and writing out the results
      titles = dict()
      rows = dict()
      headers = dict()
      txt_file_seq_text = dict()
      split_print = collections.defaultdict(list)

      for db in results:
         profile = str(self.configured_dbs[db][0])
         if results[db] == "No hit found":
            table_str += ("%s\n%s\n\n" % (profile, results[db]))
         else:
            titles[db] = "%s" % (profile)
            headers[db] = ["Resistance gene", "Identity",
                           "Alignment Length/Gene Length",
                           "Position in reference", "Contig",
                           "Position in contig", "Phenotype",
                           "Accession no."]
            table_str += ("%s\n" % (profile))
            table_str += ("Resistance gene\tIdentity\t"
                          "Alignment Length/Gene Length\t"
                          "Position in reference\tContig\tPosition in contig\t"
                          "Phenotype\tAccession no.\n")

            rows[db] = list()
            txt_file_seq_text[db] = list()

            for hit in results[db]:
               res_header = results[db][hit]["sbjct_header"]
               tmp = res_header.split("_")
               gene = tmp[0]
               acc = tmp[2]
               ID = results[db][hit]["perc_ident"]
               sbjt_length = results[db][hit]["sbjct_length"]
               HSP = results[db][hit]["HSP_length"]
               positions_contig = "%s..%s" % (results[db][hit]["query_start"],
                                              results[db][hit]["query_end"])
               positions_ref = "%s..%s" % (results[db][hit]["sbjct_start"],
                                           results[db][hit]["sbjct_end"])
               contig_name = results[db][hit]["contig_name"]
               pheno = self.phenos[gene]
               pheno = pheno.strip()

               if "split_length" in results[db][hit]:
                  total_HSP = results[db][hit]["split_length"]
                  split_print[res_header].append([gene, ID, total_HSP,
                                                  sbjt_length, positions_ref,
                                                  contig_name,
                                                  positions_contig, pheno,
                                                  acc])
                  tab_str += ("%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\n"
                              % (gene, ID, HSP, sbjt_length, positions_ref,
                                 contig_name, positions_contig, pheno, acc)
                              )
               else:
                  # Write tabels
                  table_str += ("%s\t%.2f\t%s/%s\t%s\t%s\t%s\t%s\t%s\n"
                                % (gene, ID, HSP, sbjt_length, positions_ref,
                                   contig_name, positions_contig, pheno, acc)
                                )
                  tab_str += ("%s\t%.2f\t%s/%s\t%s\t%s\t%s\t%s\t%s\n"
                              % (gene, ID, HSP, sbjt_length, positions_ref,
                                 contig_name, positions_contig, pheno, acc)
                              )

                  # Saving the output to write the txt result table
                  hsp_length = "%s/%s" % (HSP, sbjt_length)
                  rows[db].append([gene, ID, hsp_length, positions_ref,
                                   contig_name, positions_contig, pheno, acc])

               # Writing subjet/ref sequence
               ref_seq = sbjct_align[db][hit]
               ref_str += (">%s_%s\n" % (gene, acc))
               for i in range(0, len(ref_seq), 60):
                  ref_str += ("%s\n" % (ref_seq[i:i + 60]))

               # Getting the header and text for the txt file output
               sbjct_start = results[db][hit]["sbjct_start"]
               sbjct_end = results[db][hit]["sbjct_end"]
               text = ("%s, ID: %.2f %%, Alignment Length/Gene Length: %s/%s, "
                       "Positions in reference: %s..%s, Contig name: %s, "
                       "Position: %s" % (gene, ID, HSP, sbjt_length,
                                         sbjct_start, sbjct_end, contig_name,
                                         positions_contig))
               hit_str += (">%s\n" % text)

               # Writing query/hit sequence
               hit_seq = query_align[db][hit]
               for i in range(0, len(hit_seq), 60):
                  hit_str += ("%s\n" % (hit_seq[i:i + 60]))

               # Saving the output to print the txt result file allignemts
               txt_file_seq_text[db].append((text, ref_seq,
                                             homo_align[db][hit], hit_seq))

            for res in split_print:
               gene = split_print[res][0][0]
               ID = split_print[res][0][1]
               HSP = split_print[res][0][2]
               sbjt_length = split_print[res][0][3]
               positions_ref = split_print[res][0][4]
               contig_name = split_print[res][0][5]
               positions_contig = split_print[res][0][6]
               pheno = split_print[res][0][7]
               acc = split_print[res][0][8]

               for i in range(1, len(split_print[res])):
                  ID = "%s, %.2f" % (ID, split_print[res][i][1])
                  positions_ref = positions_ref + ", " + split_print[res][i][4]
                  contig_name = contig_name + ", " + split_print[res][i][5]
                  positions_contig = (positions_contig + ", "
                                      + split_print[res][i][6])

               table_str += ("%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\n"
                             % (gene, ID, HSP, sbjt_length, positions_ref,
                                contig_name, positions_contig, pheno, acc)
                             )

               hsp_length = "%s/%s" % (HSP, sbjt_length)

               rows[db].append([gene, ID, hsp_length, positions_ref,
                                contig_name, positions_contig, pheno, acc])

            table_str += ("\n")

      # Writing table txt for all hits
      for db in titles:
         # Txt file table
         table = ResFinder.text_table(titles[db], headers[db], rows[db])
         txt_str += table

      # Writing alignment txt for all hits
      for db in titles:
         # Txt file alignments
         txt_str += ("##################### %s #####################\n"
                     % (db))
         for text in txt_file_seq_text[db]:
            txt_str += ("%s\n\n" % (text[0]))
            for i in range(0, len(text[1]), 60):
               txt_str += ("%s\n" % (text[1][i:i + 60]))
               txt_str += ("%s\n" % (text[2][i:i + 60]))
               txt_str += ("%s\n\n" % (text[3][i:i + 60]))
            txt_str += ("\n")

   self.results = (tab_str, table_str, txt_str, ref_str, hit_str)

   @staticmethod
   def text_table(title, headers, rows, table_format='psql'):
      ''' Create text table

      USAGE:
         >>> from tabulate import tabulate
         >>> title = 'My Title'
         >>> headers = ['A','B']
         >>> rows = [[1,2],[3,4]]
         >>> print(text_table(title, headers, rows))
         +-----------+
         | My Title  |
         +-----+-----+
         |   A |   B |
         +=====+=====+
         |   1 |   2 |
         |   3 |   4 |
         +-----+-----+
      '''
      # Create table
      table = tabulate(rows, headers, tablefmt=table_format)
      # Prepare title injection
      width = len(table.split('\n')[0])
      tlen = len(title)
      if tlen + 4 > width:
         # Truncate oversized titles
         tlen = width - 4
         title = title[:tlen]
      spaces = width - 2 - tlen
      left_spacer = ' ' * int(spaces / 2)
      right_spacer = ' ' * (spaces - len(left_spacer))
      # Update table with title
      table = '\n'.join(['+%s+' % ('-' * (width - 2)),
                         '|%s%s%s|' % (left_spacer, title, right_spacer),
                         table, '\n'])
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

   def load_databases(self, databases):
      """
      """
      # Check if databases and config file are correct/correponds
      if databases is '':
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
      """
      """
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

      if len(self.configured_dbs) == 0:
         sys.exit("Input Error: No databases were found in the "
                  "database config file!")


if __name__ == '__main__':

   ##########################################################################
   # PARSE COMMAND LINE OPTIONS
   ##########################################################################

   parser = ArgumentParser()
   parser.add_argument("-i", "--inputfile",
                       dest="inputfile",
                       help="Input file",
                       default=None)
   parser.add_argument("-1", "--fastq1",
                       help="Raw read data file 1.",
                       default=None)
   parser.add_argument("-2", "--fastq2",
                       help="Raw read data file 2 (only required if data is \
                             paired-end).",
                       default=None)
   parser.add_argument("-o", "--outputPath",
                       dest="out_path",
                       help="Path to blast output",
                       default='')

   parser.add_argument("-b", "--blastPath",
                       dest="blast_path",
                       help="Path to blast",
                       default='blastn')
   parser.add_argument("-p", "--databasePath",
                       dest="db_path",
                       help="Path to the databases",
                       default='')
   parser.add_argument("-d", "--databases",
                       dest="databases",
                       help="Databases chosen to search in - if none are \
                             specified all are used",
                       default=None)
   parser.add_argument("-l", "--min_cov",
                       dest="min_cov",
                       help="Minimum coverage",
                       default=0.60)
   parser.add_argument("-t", "--threshold",
                       dest="threshold",
                       help="Blast threshold for identity",
                       default=0.90)
   args = parser.parse_args()

   ##########################################################################
   # MAIN
   ##########################################################################

   # Defining varibales

   min_cov = args.min_cov
   threshold = args.threshold

   # Check if valid database is provided
   if args.db_path is None:
         sys.exit("Input Error: No database directory was provided!\n")
   elif not os.path.exists(args.db_path):
      sys.exit("Input Error: The specified database directory does not "
               "exist!\n")
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

   # Check for input
   if args.inputfile is None and args.fastq1 is None:
      sys.exit("Input Error: No Input were provided!\n")

   # Check if valid input file for assembly is provided
   if args.inputfile:
      if not os.path.exists(args.inputfile):
         sys.exit("Input Error: Input file does not exist!\n")
      else:
         inputfile = args.inputfile
   else:
      inputfile = None

   # Check if valid input files for raw data
   if args.fastq1:

      if not os.path.exists(args.fastq1):
         sys.exit("Input Error: fastq1 file does not exist!\n")
      else:
         input_fastq1 = args.fastq1

      if args.fastq2:
         if not os.path.exists(args.fastq2):
            sys.exit("Input Error: fastq2 file does not exist!\n")
         else:
            input_fastq2 = args.fastq2
      else:
         input_fastq2 = None
   else:
      input_fastq1 = None
      input_fastq2 = None

   # Check if valid output directory is provided
   if not os.path.exists(args.out_path):
      sys.exit("Input Error: Output dirctory does not exists!\n")
   else:
      out_path = args.out_path

   if(inputfile is not None):
      finder = ResFinder(db_conf_file=db_config_file, databases=args.databases,
                         db_path=db_path, notes=notes_path)

      finder.blast(inputfile=inputfile, out_path=out_path, min_cov=min_cov,
                   threshold=threshold, blast=args.blast_path)
