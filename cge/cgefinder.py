#!/home/data1/tools/bin/anaconda/bin/python
import subprocess
import re


class CGEFinder():

    # Variables used by methods to distinguish results created by different
    # methods.
    TYPE_BLAST = "blast"
    TYPE_KMA = "kma"

    @staticmethod
    def kma(inputfile_1, out_path, databases, db_path_kma, min_cov=0.9,
            threshold=0.6, kma_path="cge/kma/kma", sample_name="",
            inputfile_2=None, kma_mrs=None, kma_gapopen=None,
            kma_gapextend=None, kma_penalty=None, kma_reward=None):
        """
           I expect that there will only be one hit pr gene, but if there are
           more, I assume that the sequence of the hits are the same in the res
           file and the aln file.
        """

        kma_results = dict()
        kma_results["excluded"] = dict()

        if(sample_name):
           sample_name = "_" + sample_name

        for db in databases:
           kma_db = db_path_kma + "/" + db
           kma_outfile = out_path + "/kma_" + db + sample_name
           kma_cmd = ("%s -t_db %s -SW -o %s -e 1.0 -i %s" % (kma_path, kma_db,
                      kma_outfile, inputfile_1))
           if(inputfile_2 is not None):
              kma_cmd += " " + inputfile_2
           if(kma_mrs is not None):
               kma_cmd += " -mrs " + str(kma_mrs)
           if(kma_gapopen is not None):
               kma_cmd += " -gapopen " + str(kma_gapopen)
           if(kma_gapextend is not None):
               kma_cmd += " -gapextend " + str(kma_gapextend)
           if(kma_gapextend is not None):
               kma_cmd += " -gapextend " + str(kma_gapextend)
           if(kma_penalty is not None):
               kma_cmd += " -penalty " + str(kma_penalty)
           if(kma_reward is not None):
               kma_cmd += " -reward " + str(kma_reward)

           # Call KMA
           process = subprocess.Popen(kma_cmd, shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
           out, err = process.communicate()

           kma_results[db] = 'No hit found'

           # Fetch kma output files
           align_filename = kma_outfile + ".aln"
           res_filename = kma_outfile + ".res"

           # Open res file, find coverage and the gene names of genes found
           with open(res_filename, "r") as res_file:
              header = res_file.readline()

              for line in res_file:

                 if kma_results[db] == 'No hit found':
                    kma_results[db] = dict()
                    # kma_results[db]["excluded"] = dict()
                    # continue

                 data = [data.strip() for data in line.split("\t")]
                 gene = data[0]

                 sbjct_len = int(data[3])
                 sbjct_ident = float(data[4])
                 coverage = float(data[5])

                 if gene not in kma_results[db]:
                    hit = gene
                 else:
                    hit = gene + "_" + str(len(kma_results[db][gene]) + 1)

                 exclude_reasons = []
                 if(coverage < min_cov):
                     exclude_reasons.append("coverage: " + str(coverage))
                 elif(sbjct_ident < threshold):
                     exclude_reasons.append("identity: " + str(sbjct_ident))

                 if(exclude_reasons):
                     # kma_results[db]["excluded"][hit] = exclude_reasons
                     kma_results["excluded"][hit] = exclude_reasons

                 kma_results[db][hit] = dict()
                 kma_results[db][hit]['sbjct_length'] = sbjct_len
                 kma_results[db][hit]["perc_coverage"] = coverage
                 kma_results[db][hit]["sbjct_string"] = []
                 kma_results[db][hit]["query_string"] = []
                 kma_results[db][hit]["homology"] = []
                 kma_results[db][hit]["sbjct_header"] = gene
                 kma_results[db][hit]["perc_ident"] = sbjct_ident
                 kma_results[db][hit]["query_start"] = "NA"
                 kma_results[db][hit]["query_end"] = "NA"
                 kma_results[db][hit]["contig_name"] = "NA"
                 kma_results[db][hit]["HSP_length"] = "NA"

           if kma_results[db] == 'No hit found':
              continue

           # Open align file
           with open(align_filename, "r") as align_file:
              hit_no = dict()
              gene = ""
              # Parse through alignments
              for line in align_file:
                 # Skip empty lines
                 if(not line.strip()):
                      continue

                 # Check when a new gene alignment start
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

                    if hit in kma_results[db]:
                       line_data = line.split("\t")[-1].strip()
                       if line.startswith("template"):
                          kma_results[db][hit]["sbjct_string"] += [line_data]
                       elif line.startswith("query"):
                          kma_results[db][hit]["query_string"] += [line_data]
                       else:
                          kma_results[db][hit]["homology"] += [line_data]
                    else:
                       print(hit + " not in results: ", kma_results)

           # concatinate all sequences lists and find subject start and subject
           # end
           seq_start_search_str = re.compile("^-*(\w+)")

           for hit in kma_results[db]:
              # if(hit == "excluded"):
              # continue
              kma_results[db][hit]['sbjct_string'] = "".join(
                  kma_results[db][hit]['sbjct_string'])
              kma_results[db][hit]['query_string'] = "".join(
                  kma_results[db][hit]['query_string'])
              kma_results[db][hit]['homology'] = "".join(
                  kma_results[db][hit]['homology'])

              seq_start_object = seq_start_search_str.search(
                  kma_results[db][hit]['query_string'])
              sbjct_start = seq_start_object.start() + 1
              kma_results[db][hit]['sbjct_start'] = sbjct_start
              kma_results[db][hit]["sbjct_end"] = (
                  kma_results[db][hit]["sbjct_length"] - sbjct_start + 1)

        return kma_results
