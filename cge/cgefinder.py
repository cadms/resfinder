#!/usr/bin/env python3
import subprocess
import re
import os.path

# TODO import blaster and make blaster function in CGEFinder
# from cge.blaster.blaster import Blaster


class FinderResult():
    def __init__(self, results, align_sbjct=None, align_query=None,
                 align_homo=None):
        self.results = results  # Results
        self.gene_align_query = align_query  # Sequence alignment lines
        self.gene_align_homo = align_homo  # Sequence alignment homolog string
        self.gene_align_sbjct = align_sbjct  # Sequence alignment allele string


class CGEFinder():

    # Variables used by methods to distinguish results created by different
    # methods.
    TYPE_BLAST = "blast"
    TYPE_KMA = "kma"

    @staticmethod
    def kma(inputfile_1, out_path, databases, db_path_kma, min_cov=0.9,
            threshold=0.6, kma_path="cge/kma/kma", sample_name="",
            inputfile_2=None, kma_mrs=None, kma_gapopen=None,
            kma_gapextend=None, kma_penalty=None, kma_reward=None, kma_pm=None,
            kma_fpm=None):
        """
           I expect that there will only be one hit pr gene, but if there are
           more, I assume that the sequence of the hits are the same in the res
           file and the aln file.
        """
        threshold = threshold * 100
        min_cov = min_cov * 100

        kma_results = dict()
        kma_results["excluded"] = dict()

        if(sample_name):
           sample_name = "_" + sample_name

        for db in databases:
            kma_db = db_path_kma + "/" + db
            kma_outfile = out_path + "/kma_" + db + sample_name
            kma_cmd = ("%s -t_db %s -o %s -e 1.0" % (kma_path,
                       kma_db, kma_outfile))
            if(inputfile_2 is not None):
                kma_cmd += " -ipe " + inputfile_1 + " " + inputfile_2
            else:
                kma_cmd += " -i " + inputfile_1
            if(kma_mrs is not None):
                kma_cmd += " -mrs " + str(kma_mrs)
            if(kma_gapopen is not None):
                kma_cmd += " -gapopen " + str(kma_gapopen)
            if(kma_gapextend is not None):
                kma_cmd += " -gapextend " + str(kma_gapextend)
            if(kma_penalty is not None):
                kma_cmd += " -penalty " + str(kma_penalty)
            if(kma_reward is not None):
                kma_cmd += " -reward " + str(kma_reward)
            if(kma_pm is not None):
                kma_cmd += " -pm " + kma_pm
            if(kma_fpm is not None):
                kma_cmd += " -fpm " + kma_fpm

            # kma output files
            align_filename = kma_outfile + ".aln"
            res_filename = kma_outfile + ".res"

            # If .res file exists then skip mapping
            if(os.path.exists(res_filename)):
                print("Found " + res_filename + " skipping DB.")
            else:
                # Call KMA
                process = subprocess.Popen(kma_cmd, shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                out, err = process.communicate()

            kma_results[db] = 'No hit found'

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
                    q_value = float(data[-2])

                    if gene not in kma_results[db]:
                        hit = gene
                    else:
                        hit = gene + "_" + str(len(kma_results[db][gene]) + 1)

                    exclude_reasons = []

                    if(coverage < min_cov or sbjct_ident < threshold):
                        exclude_reasons.append(coverage)
                        exclude_reasons.append(sbjct_ident)

                    if(exclude_reasons):
                        # kma_results[db]["excluded"][hit] = exclude_reasons
                        kma_results["excluded"][hit] = exclude_reasons

                    kma_results[db][hit] = dict()
                    kma_results[db][hit]['sbjct_length'] = sbjct_len
                    kma_results[db][hit]["perc_coverage"] = coverage
                    kma_results[db][hit]["sbjct_string"] = []
                    kma_results[db][hit]["query_string"] = []
                    kma_results[db][hit]["homo_string"] = []
                    kma_results[db][hit]["sbjct_header"] = gene
                    kma_results[db][hit]["perc_ident"] = sbjct_ident
                    kma_results[db][hit]["query_start"] = "NA"
                    kma_results[db][hit]["query_end"] = "NA"
                    kma_results[db][hit]["contig_name"] = "NA"
                    kma_results[db][hit]["HSP_length"] = ""
                    kma_results[db][hit]["cal_score"] = q_value

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
                                kma_results[db][hit]["sbjct_string"] += (
                                    [line_data])
                            elif line.startswith("query"):
                                kma_results[db][hit]["query_string"] += (
                                    [line_data])
                            else:
                                kma_results[db][hit]["homo_string"] += (
                                    [line_data])
                        else:
                            print(hit + " not in results: ", kma_results)

            # concatinate all sequences lists and find subject start
            # and subject end

            gene_align_sbjct = {db: {}}
            gene_align_query = {db: {}}
            gene_align_homo = {db: {}}

            for hit in kma_results[db]:
                # if(hit == "excluded"):
                # continue
                align_sbjct = "".join(kma_results[db][hit]['sbjct_string'])
                align_query = "".join(kma_results[db][hit]['query_string'])
                align_homo = "".join(kma_results[db][hit]['homo_string'])

                # Extract only aligned sequences
                start = re.search("^-*(\w+)", align_query).start(1)
                end = re.search("\w+(-*)$", align_query).start(1)

                kma_results[db][hit]['sbjct_string'] = align_sbjct[start:end]
                kma_results[db][hit]['query_string'] = align_query[start:end]
                kma_results[db][hit]['homo_string'] = align_homo[start:end]

                # Save align start and stop positions relative to
                # subject sequence
                kma_results[db][hit]['sbjct_start'] = start + 1
                kma_results[db][hit]["sbjct_end"] = end + 1
                kma_results[db][hit]["HSP_length"] = end - start

                # Count gaps in the alignment
                kma_results[db][hit]["gaps"] = (
                    kma_results[db][hit]['sbjct_string'].count("-") +
                    kma_results[db][hit]['query_string'].count("-"))

                # Save sequences covering the entire subject sequence
                # in seperate variables
                gene_align_sbjct[db][hit] = align_sbjct
                gene_align_query[db][hit] = align_query
                gene_align_homo[db][hit] = align_homo

        return FinderResult(kma_results, gene_align_sbjct, gene_align_query,
                            gene_align_homo)
