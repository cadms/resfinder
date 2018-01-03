#!/home/data1/tools/bin/anaconda/bin/python
import sys
import os
import re
import math
from argparse import ArgumentParser

from cge.blaster import *


##########################################################################
# FUNCTIONS
##########################################################################


def get_db_mutations(mut_db_path, gene_list):
    """
    This function opens the file resistenss-overview.txt, and reads the
    content into a dict of dicts. The dict will contain information about
    all known mutations given in the database. This dict is returned.
    """

    # Open resistens-overview.txt
    try:
        drugfile = open(mut_db_path, "r")
    except:
        sys.exit("wrong path", mut_db_path)

    # Initiate dict
    known_mutations = dict()

    indelflag = False
    # Go throug the file line by line
    for line in drugfile:
        # Ignore headers and check when the indel section starts
        if line.startswith("#"):
            if "indel" in line.lower():
                indelflag = True
            continue

        # Assert that all lines have the correct set of columns
        mutation = [data.strip() for data in line.strip().split("\t")]
        assert len(mutation) == 9, "mutation overview file (%s) must have 9 columns, %s" % (mut_db_path, mutation)

        # Extract all info on the line (even though it is not all used)
        gene_ID = mutation[0]

        # Only consider mutations in genes found in the gene list
        if gene_ID in gene_list:
            gene_name = mutation[1]
            no_of_mut = int(mutation[2])
            mut_pos = int(mutation[3])
            ref_codon = mutation[4]
            ref_aa = mutation[5]
            alt_aa = mutation[6].split(",")
            res_drug = mutation[7].replace("\t", " ")
            pmid = mutation[8].split(",")

            # Initiate empty dict to store relevant mutation information
            mut_info = dict()

            # Save need mutation info with pmid cooresponding to the amino
            # acid change
            for i in range(len(alt_aa)):
                try:
                    mut_info[alt_aa[i]] = {"gene_name": gene_name,
                                           "drug": res_drug,
                                           "pmid": pmid[i]}

                except IndexError:
                    print("pmid not found, %s") % line
                    mut_info[alt_aa[i]] = {"gene_name": gene_name,
                                           "drug": res_drug,
                                           "pmid": "-"}

            # Check if more than one mutations is needed for resistance
            if no_of_mut != 1:
                print("More than one mutation is needed, this is not "
                      + "implemented", mutation)

            # Add all possible types of mutations to the dict
            if gene_ID not in known_mutations:
                known_mutations[gene_ID] = {"sub": dict(), "ins": dict(),
                                            "del": dict()}

            # Check for the type of mutation
            if indelflag is False:
                mutation_type = "sub"
            else:
                mutation_type = ref_aa

            # Save mutations positions with required information given in
            # mut_info
            if mut_pos not in known_mutations[gene_ID][mutation_type]:
                known_mutations[gene_ID][mutation_type][mut_pos] = dict()
            for aa in alt_aa:
                known_mutations[gene_ID][mutation_type][mut_pos][aa] = mut_info[aa]

    drugfile.close()

    for gene in gene_list:
        if gene not in known_mutations:
            known_mutations[gene] = {"sub": dict(), "ins": dict(),
                                     "del": dict()}

    return known_mutations


def find_best_sequence(hits_found):
    """
    This function takes the list hits_found as argument. This contains all
    hits found for the blast search of one gene. A hit includes the subjct
    sequence, the query, and the start and stop position of the allignment
    corresponding to the subject sequence. This function finds the best
    hit by concatinating sequences of found hits. If different overlap
    sequences occurr these are saved in the list alternative_overlaps. The
    subject and query sequence of the concatinated sequence to gether with
    alternative overlaps and the corresponding start stop
    positions are returned.
    """

    # Get information from the fisrt hit found
    all_start = hits_found[0][0]
    current_end = hits_found[0][1]
    final_sbjct = hits_found[0][2]
    final_qry = hits_found[0][3]
    alternative_overlaps = []

    # Check if more then one hit was found within the same gene
    for i in range(len(hits_found) - 1):
        # Save information from previous hit
        pre_block_start = hits_found[i][0]
        pre_block_end = hits_found[i][1]
        pre_sbjct = hits_found[i][2]
        pre_qry = hits_found[i][3]
        # Save information from next hit
        next_block_start = hits_found[i + 1][0]
        next_block_end = hits_found[i + 1][1]
        next_sbjct = hits_found[i + 1][2]
        next_qry = hits_found[i + 1][3]

        # Check for overlapping sequences, collaps them and save alternative
        # overlaps if any
        if next_block_start <= current_end:
            # Find overlap start and take gaps into account
            pos_count = 0
            overlap_pos = pre_block_start
            for i in range(len(pre_sbjct)):
                # Stop loop if overlap_start position is reached
                if overlap_pos == next_block_start:
                    overlap_start = pos_count
                    break
                if pre_sbjct[i] != "-":
                    overlap_pos += 1
                pos_count += 1

            # Find overlap length and add next sequence to final sequence
            if len(pre_sbjct[overlap_start:]) > len(next_sbjct):
                #  <--------->
                #     <--->
                overlap_len = len(next_sbjct)
                overlap_end_pos = next_block_end
            else:
                #  <--------->
                #        <--------->
                overlap_len = len(pre_sbjct[overlap_start:])
                overlap_end_pos = pre_block_end
                # Update current end
                current_end = next_block_end
                # Use the entire pre sequence and add the last part of the
                # next sequence
                final_sbjct += next_sbjct[overlap_len:]
                final_qry += next_qry[overlap_len:]

            # Find query overlap sequences
            pre_qry_overlap = pre_qry[overlap_start: (overlap_start + overlap_len)]  # can work for both types of overlap
            next_qry_overlap = next_qry[:overlap_len]
            sbjct_overlap = next_sbjct[:overlap_len]
            # If alternative query overlap excist save it
            if pre_qry_overlap != next_qry_overlap:
                print("OVERLAP WARNING:")
                print(pre_qry_overlap, "\n", next_qry_overlap)
                # save alternative overlaps
                alternative_overlaps += [(next_block_start, overlap_end_pos,
                                          sbjct_overlap, next_qry_overlap)]

        elif next_block_start > current_end:
            #  <------->
            #              <------->
            gap_size = next_block_start - current_end - 1
            for i in range(gap_size):
                final_sbjct += "N"
                final_qry += "N"
            current_end = next_block_end
            final_sbjct += next_sbjct
            final_qry += next_qry

    return final_sbjct, final_qry, all_start, current_end, alternative_overlaps


def find_mismatches(gene, sbjct_start, sbjct_seq, qry_seq,
                    alternative_overlaps=[]):
    """
    This function finds mis matches between two sequeces. Depending on the
    the sequence type either the function find_codon_mismatches or
    find_nucleotid_mismatches are called, if the sequences contains both
    a promoter and a coding region both functions are called. The function
    can also call it self if alternative overlaps is give. All found mis
    matches are returned
    """

    # Initiate the mis_matches list that will store all found mis matcehs
    mis_matches = []

    # Find mis matches in RNA genes
    if gene in RNA_gene_list:
        mis_matches += find_nucleotid_mismatches(sbjct_start, sbjct_seq,
                                                 qry_seq)
    else:
        # Check if the gene sequence is with a promoter
        regex = r"promoter_size_(\d+)(?:bp)"
        promtr_gene_objt = re.search(regex, gene)

        # Check for promoter sequences
        if promtr_gene_objt:
            # Get promoter length
            promtr_len = int(promtr_gene_objt.group(1))

            # Extract promoter sequence, while considering gaps
            # --------agt-->----
            #    ---->?
            if sbjct_start <= promtr_len:
                # Find position in sbjct sequence where promoter ends
                promtr_end = 0
                nuc_count = sbjct_start - 1
                for i in range(len(sbjct_seq)):
                    promtr_end += 1
                    if sbjct_seq[i] != "-":
                        nuc_count += 1
                    if nuc_count == promtr_len:
                        break

                # Check if only a part of the promoter is found
                # --------agt-->----
                # ----
                promtr_sbjct_start = -1
                if nuc_count < promtr_len:
                    promtr_sbjct_start = nuc_count - promtr_len

                # Get promoter part of subject and query
                sbjct_promtr_seq = sbjct_seq[:promtr_end]
                qry_promtr_seq = qry_seq[:promtr_end]

                # For promoter part find nucleotide mis matches
                mis_matches += find_nucleotid_mismatches(promtr_sbjct_start,
                                                         sbjct_promtr_seq,
                                                         qry_promtr_seq,
                                                         promoter=True)

                # Check if gene is also found
                # --------agt-->----
                #     -----------
                if (sbjct_start + len(sbjct_seq.replace("-", ""))) > promtr_len:
                    sbjct_gene_seq = sbjct_seq[promtr_end:]
                    qry_gene_seq = qry_seq[promtr_end:]
                    sbjct_gene_start = 1

                    # Find mismatches in gene part
                    mis_matches += find_codon_mismatches(sbjct_gene_start,
                                                         sbjct_gene_seq,
                                                         qry_gene_seq)

            # No promoter, only gene is found
            # --------agt-->----
            #            -----
            else:
                sbjct_gene_start = sbjct_start - promtr_len

                # Find mismatches in gene part
                mis_matches += find_codon_mismatches(sbjct_gene_start,
                                                     sbjct_seq, qry_seq)

        else:
            # Find mismatches in gene
            mis_matches += find_codon_mismatches(sbjct_start, sbjct_seq,
                                                 qry_seq)

    # Find mismatches in alternative overlaps if any
    for overlap in alternative_overlaps:
        mis_matches += find_mismatches(gene, overlap[0], overlap[2],
                                       overlap[3])

    return mis_matches


def find_nucleotid_mismatches(sbjct_start, sbjct_seq, qry_seq, promoter=False):
    """
    This function takes two alligned sequence (subject and query), and the
    position on the subject where the alignment starts. The sequences are
    compared one nucleotide at a time. If mis matches are found they are
    saved. If a gap is found the function find_nuc_indel is called to find
    the entire indel and it is also saved into the list mis_matches. If
    promoter sequences are given as arguments, these are reversed the and
    the absolut value of the sequence position  used, but when mutations
    are saved the negative value and det reverse sequences are saved in
    mis_mathces.
    """

    # Initiate the mis_matches list that will store all found mis matcehs
    mis_matches = []

    sbjct_start = abs(sbjct_start)
    seq_pos = sbjct_start

    # Set variables depending on promoter status
    factor = 1
    mut_prefix = "r."
    if promoter is True:
        factor = (-1)
        mut_prefix = "n."
        # Reverse promoter sequences
        sbjct_seq = sbjct_seq[::-1]
        qry_seq = qry_seq[::-1]

    # Go through sequences one nucleotide at a time
    shift = 0
    for index in range(sbjct_start - 1, len(sbjct_seq)):
        mut_name = mut_prefix
        mut = ""
        # Shift index according to gaps
        i = index + shift

        # If the end of the sequence is reached, stop
        if i == len(sbjct_seq):
            break

        sbjct_nuc = sbjct_seq[i]
        qry_nuc = qry_seq[i]

        # Check for mis matches
        if sbjct_nuc != qry_nuc:

            # check for insertions and deletions
            if sbjct_nuc == "-" or qry_nuc == "-":
                if sbjct_nuc == "-":
                    mut = "ins"
                    indel_start_pos = (seq_pos - 1) * factor
                    indel_end_pos = seq_pos * factor
                    indel = find_nuc_indel(sbjct_seq[i:], qry_seq[i:])
                else:
                    mut = "del"
                    indel_start_pos = seq_pos * factor
                    indel = find_nuc_indel(qry_seq[i:], sbjct_seq[i:])
                    indel_end_pos = (seq_pos + len(indel) - 1) * factor
                    seq_pos += len(indel) - 1

                # Shift the index to the end of the indel
                shift += len(indel) - 1

                # Write mutation name, depending on sequnce
                if len(indel) == 1 and mut == "del":
                    mut_name += str(indel_start_pos) + mut + indel
                else:
                    if promoter is True:
                        # Reverse the sequence and the start and end positions
                        indel = indel[::-1]
                        temp = indel_start_pos
                        indel_start_pos = indel_end_pos
                        indel_end_pos = temp

                    mut_name += (str(indel_start_pos) + "_"
                                 + str(indel_end_pos) + mut + indel)

                mis_matches += [(mut, seq_pos * factor, seq_pos * factor,
                                 indel, mut_name, mut, indel)]

            # Check for substitutions mutations
            else:
                mut = "sub"
                mut_name += str(seq_pos * factor) + sbjct_nuc + ">" + qry_nuc
                mis_matches += [(mut, seq_pos * factor, seq_pos * factor,
                                 qry_nuc, mut_name, sbjct_nuc, qry_nuc)]
        # Increment sequence position
        if mut != "ins":
            seq_pos += 1

    return mis_matches


def find_nuc_indel(gapped_seq, indel_seq):
    """
    This function finds the entire indel missing in from a gapped sequence
    compared to the indel_seqeunce. It is assumes that the sequences start
    with the first position of the gap.
    """
    ref_indel = indel_seq[0]
    for j in range(1, len(gapped_seq)):
        if gapped_seq[j] == "-":
            ref_indel += indel_seq[j]
        else:
            break
    return ref_indel


def aa(codon, mut="?"):
    """
    This function converts a codon to an amino acid. If the codon is not
    valid an error message is given, or else, the amino acid is returned.
    """
    codon = codon.upper()
    aa = {"ATT": "I", "ATC": "I", "ATA": "I",
          "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L",
          "TTG": "L",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
          "TTT": "F", "TTC": "F",
          "ATG": "M",
          "TGT": "C", "TGC": "C",
          "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
          "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S",
          "AGC": "S",
          "TAT": "Y", "TAC": "Y",
          "TGG": "W",
          "CAA": "Q", "CAG": "Q",
          "AAT": "N", "AAC": "N",
          "CAT": "H", "CAC": "H",
          "GAA": "E", "GAG": "E",
          "GAT": "D", "GAC": "D",
          "AAA": "K", "AAG": "K",
          "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
          "AGG": "R",
          "TAA": "*", "TAG": "*", "TGA": "*"}

    try:
        amino_a = aa[codon]
    except KeyError:
        if len(codon) != 3:
            amino_a = "?"
            print("ERROR, this is not a codon", codon)
        elif "---" in codon and mut == "ins":
            amino_a = mut
        else:
            amino_a = "?"
    return amino_a


def get_codon(seq, codon_no):
    """
    This function takes a sequece and a codon number and returns the codon
    found in the sequence at that position
    """
    seq = seq.replace("-", "")
    codon = seq[int(codon_no - 1) * 3:int(codon_no - 1) * 3 + 3]
    return codon


def translate_indel(sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel,
                    codon_no, mut):
    """
    This function serves to give a name to insertion and deletion
    mutations found in the coding region of a gene. The function detects
    if a frameshift occurred, and this will be shown in the name if the
    mutation. If an insert occurr the insert is named after the two amino
    acids inbetween the mutation. If it is a deletion the first and the
    last amino acid deleted is used in the name of the mutation, unless
    only one amino acid is affected
    """

    # Get the indel sequences without gaps
    sbjct_nucs = sbjct_rf_indel.replace("-", "")
    qry_nucs = qry_rf_indel.replace("-", "")
    if qry_rf_indel == "-":
        print("HERE", sbjct_nucs, qry_nucs)


    if mut == "ins" and (sbjct_rf_indel.count("-") > 2 or len(sbjct_rf_indel) == sbjct_rf_indel.count("-")):
        start_codon_no = codon_no - 1
        start_codon = get_codon(sbjct_seq, start_codon_no)
        end_codon = get_codon(sbjct_seq, codon_no)
        pos_name = (aa(start_codon, mut) + str(start_codon_no) + "_"
                    + aa(end_codon, mut) + str(codon_no) + mut)

    else:
        pos_name = aa(sbjct_nucs[:3]) + str(codon_no)
        if (mut == "ins" and len(sbjct_nucs) == 3) or (mut == "del" and len(qry_nucs) == 3):  # mut == "ins" or len(sbjct_nucs) == 3:
            pos_name += aa(qry_nucs, mut)
        elif mut == "del" and len(sbjct_nucs) > 3:
            end_codon_no = codon_no + math.ceil(len(sbjct_nucs) / 3) - 1
            end_codon = get_codon(sbjct_seq, end_codon_no)
            pos_name += "_" + aa(end_codon, mut) + str(int(end_codon_no)) + mut
        elif mut == "del":
            pos_name += mut

    # Set mut_name
    mut_name = "p." + pos_name

    # Translate indel to amino acids
    aa_ref = ""
    aa_alt = ""
    for i in range(0, len(sbjct_nucs), 3):
        aa_ref += aa(sbjct_nucs[i:i + 3])
    for i in range(0, len(qry_nucs), 3):
        aa_alt += aa(qry_nucs[i:i + 3])
    if mut == "ins" and sbjct_nucs == "":
        mut_name += aa_alt
        aa_ref = mut
    elif mut == "del" and qry_nucs == "" and len(sbjct_nucs) > 3:
        mut_name += aa_ref
        aa_alt = mut

    # Check for frameshift mutations
    if (mut == "ins" and sbjct_rf_indel.count("-") % 3 != 0) or (mut == "del" and qry_rf_indel.count("-") % 3 != 0):
        mut_name += " - Frameshift"

    return mut_name, aa_ref, aa_alt


def get_inframe_gap(seq, nucs_needed=3):
    """
    This funtion takes a sequnece starting with a gap or the complementary
    seqeuence to the gap, and the number of nucleotides that the seqeunce
    should contain in order to maintain the correct reading frame. The
    sequence is gone through and the number of non-gap characters are
    counted. When the number has reach the number of needed nucleotides
    the indel is returned. If the indel is a 'clean' insert or deletion
    that starts in the start of a codon and can be divided by 3, then only
    the gap is returned.
    """

    nuc_count = 0
    gap_indel = ""
    for i in range(len(seq)):
        # Check if the character is not a gap
        if seq[i] != "-":
            # Check in the indel is a 'clean' insert or deletion that starts
            # in the start of a codon and can be divided by 3
            if gap_indel.count("-") == len(gap_indel) and len(gap_indel) % 3 == 0 and len(gap_indel) != 0:
                return gap_indel
            nuc_count += 1

        gap_indel += seq[i]
        # if the number of nucleotides in the indel equals the amount needed
        # the indel is returned
        if nuc_count == nucs_needed:
            return gap_indel
    # This will only happen if the gap is in the very end of a sequence
    return gap_indel


def get_indels(sbjct_seq, qry_seq, start_pos):
    """
    This function uses regex to find inserts and deletions in sequences
    given as arguments. A list of these indels are returned. The list
    includes, type of mutations(ins/del), subject codon no of found
    mutation, subject sequence position, insert/deletions nucleotide
    sequence, and the affected qry codon no.
    """

    seqs = [sbjct_seq, qry_seq]
    indels = []
    gap_obj = re.compile(r"-+")
    for i in range(len(seqs)):
        for match in gap_obj.finditer(seqs[i]):
            pos = int(match.start())
            gap = match.group()

            # find position of the mutation corresponding to the subject
            # sequence
            sbj_pos = len(sbjct_seq[:pos].replace("-", "")) + start_pos - 1

            # get indel sequence and the affected sequences in sbjct and qry
            # in the reading frame
            indel = seqs[abs(i - 1)][pos:pos + len(gap)]

            # find codon number for mutation
            codon_no = int(math.ceil((sbj_pos) / 3))
            qry_pos = len(qry_seq[:pos].replace("-", "")) + start_pos - 1
            qry_codon = int(math.ceil((qry_pos) / 3))
            if i == 0:
                mut = "ins"
            else:
                mut = "del"
            indels.append([mut, codon_no, sbj_pos, indel, qry_codon])
    indels = sorted(indels, key=lambda x: (x[1], x[2]))

    return indels


##########################################################################
def find_codon_mismatches(sbjct_start, sbjct_seq, qry_seq):
    """
    This function takes two alligned sequence (subject and query), and
    the position on the subject where the alignment starts. The sequences
    are compared codon by codon. If a mis matches is found it is saved in
    'mis_matches'. If a gap is found the function get_inframe_gap is used
    to find the indel sequence and keep the sequence in the correct
    reading frame. The function translate_indel is used to name indel
    mutations and translate the indels to amino acids
    The function returns a list of tuples containing all needed informations
    about the mutation in order to look it up in the database dict known
    mutation and the with the output files the the user.
    """
    mis_matches = []

    # Find start pos of first codon in frame, i_start
    codon_offset = (sbjct_start - 1) % 3
    i_start = 0
    if codon_offset != 0:
        i_start = 3 - codon_offset
    sbjct_start = sbjct_start + i_start

    # Set sequences in frame
    sbjct_seq = sbjct_seq[i_start:]
    qry_seq = qry_seq[i_start:]

    # Find codon number of the first codon in the sequence, start at 0
    codon_no = int((sbjct_start - 1) / 3)  # 1,2,3 start on 0

    # s_shift and q_shift are used when gaps appears
    q_shift = 0
    s_shift = 0
    s_gaps = 0
    q_gaps = 0
    mut_no = 0

    # Find inserts and deletions in sequence
    indel_no = 0
    indels = get_indels(sbjct_seq, qry_seq, sbjct_start)

    # Go through sequence and save mutations when found
    for index in range(0, len(sbjct_seq), 3):
        # Count codon number
        codon_no += 1

        # Shift index according to gaps
        s_i = index + s_shift
        q_i = index + q_shift

        # get codons
        sbjct_codon = sbjct_seq[s_i:s_i + 3]
        qry_codon = qry_seq[q_i:q_i + 3]

        # print(codon_no, sbjct_codon, qry_codon)#, aa(qry_codon))
        if len(sbjct_seq[s_i:].replace("-", "")) + len(qry_codon[q_i:].replace("-", "")) < 6:
            break
        # print(sbjct_codon, qry_codon)
        # Check for mutations
        if sbjct_codon != qry_codon:
            # print(codon_no, int(codon_no*3) - 2, sbjct_codon, aa(sbjct_codon), qry_codon, aa(qry_codon))
            # Check for codon insertions and deletions and frameshift mutations
            if "-" in sbjct_codon or "-" in qry_codon:
                print(sbjct_codon, qry_codon)
                # print(indels)
                # get indel info
                try:
                    indel_data = indels[indel_no]
                except IndexError:
                    print(sbjct_codon, qry_codon)
                    print(indels)
                    print(gene, indel_data, indel_no)
                    sys.exit("indel_data list is out of range, bug!")
                mut = indel_data[0]
                codon_no_indel = indel_data[1]
                seq_pos = indel_data[2] + sbjct_start - 1
                indel = indel_data[3]
                indel_no += 1

                # Get the affected sequence in frame for both for sbjct and qry
                if mut == "ins":
                    sbjct_rf_indel = get_inframe_gap(sbjct_seq[s_i:])
                    qry_rf_indel = get_inframe_gap(qry_seq[q_i:], int(math.floor(len(sbjct_rf_indel) / 3) * 3))
                    N_mut = sbjct_rf_indel.upper().replace("A", "N").replace("T", "N").replace("C", "N").replace("G", "N")
                else:
                    qry_rf_indel = get_inframe_gap(qry_seq[q_i:])
                    sbjct_rf_indel = get_inframe_gap(sbjct_seq[s_i:], int(math.floor(len(qry_rf_indel) / 3) * 3))
                    N_mut = qry_rf_indel.upper().replace("A", "N").replace("T", "N").replace("C", "N").replace("G", "N")
                # print(sbjct_rf_indel, qry_rf_indel)
                no_of_indels = N_mut.count("-N")
                # count gaps
                s_gaps += sbjct_rf_indel.count("-")
                q_gaps += qry_rf_indel.count("-")

                mut_name, aa_ref, aa_alt = translate_indel(sbjct_seq, indel,
                                                           sbjct_rf_indel,
                                                           qry_rf_indel,
                                                           codon_no, mut)

                # set index to the correct reading frame after the indel gap
                shift_diff_before = abs(s_shift - q_shift)
                s_shift += len(sbjct_rf_indel) - 3
                q_shift += len(qry_rf_indel) - 3
                shift_diff = abs(s_shift - q_shift)

                # print("I AM HERE!!", sbjct_rf_indel, qry_rf_indel)
                # Check if the next mutation in the indels list is in the
                # current codon
                if no_of_indels > 1:
                    mut_no += 1
                    if mut_no == no_of_indels:
                        mut_no = 0
                    else:
                        s_shift -= len(sbjct_rf_indel)
                        q_shift -= len(qry_rf_indel)
                        codon_no -= 1
                        s_gaps -= sbjct_rf_indel.count("-")
                        q_gaps -= qry_rf_indel.count("-")
                        mis_matches += [(mut, codon_no_indel, seq_pos, indel,
                                         mut_name, sbjct_rf_indel,
                                         qry_rf_indel, aa_ref, aa_alt)]
                        continue

                # Check if the indel restores the reading frame, if yes,
                # change the last mutation
                if shift_diff_before != 0 and shift_diff % 3 == 0:
                    print(shift_diff_before, shift_diff)
                    if mut == "ins" and shift_diff >= 3:
                        sbjct_rf_indel = "-" * sbjct_rf_indel.count("-")
                        aa_ref = mut
                        s_shift = max(q_shift, s_shift)
                        q_shift = max(q_shift, s_shift)
                    elif mut == "del" and shift_diff >= 3:

                        qry_rf_indel = "-" * qry_rf_indel.count("-")
                        aa_alt = mut
                        s_shift = max(q_shift, s_shift)
                        q_shift = max(q_shift, s_shift)
                    print(mut_name)
                    mut_name, dont_use, dont_use2 = translate_indel(
                        sbjct_seq, indel, sbjct_rf_indel, qry_rf_indel,
                        codon_no, mut)
                    print(mut_name)

                    if "Frameshift" in mut_name:
                        mut_name = mut_name.split("-")[0] + "- Frame restored"

                mis_matches += [(mut, codon_no_indel, seq_pos, indel, mut_name,
                                 sbjct_rf_indel, qry_rf_indel, aa_ref, aa_alt)]

                # set codon number, and save nucleotides from out of frame
                # mutations
                if mut == "del":
                    codon_no += int((len(sbjct_rf_indel) - 3) / 3)
                # if evaluated insert is only gaps codon_no should not
                # increment
                elif sbjct_rf_indel.count("-") == len(sbjct_rf_indel):
                    codon_no -= 1

            # Check of point mutations
            else:
                mut = "sub"
                aa_ref = aa(sbjct_codon)
                aa_alt = aa(qry_codon)

                if aa_ref != aa_alt:
                    # End search for mutation if a premature stop codon is
                    # found
                    mut_name = "p." + aa_ref + str(codon_no) + aa_alt

                    # if a stop codon is found in the query, stop searching for
                    # additional mutations
                    if aa_alt == "*":
                        mut_name += " - Premature stop codon"
                        mis_matches += [(mut, codon_no, codon_no, aa_alt,
                                         mut_name, sbjct_codon, qry_codon,
                                         aa_ref, aa_alt)]
                        break
                    mis_matches += [(mut, codon_no, codon_no, aa_alt, mut_name,
                                     sbjct_codon, qry_codon, aa_ref, aa_alt)]

    # sort mutations on position
    mis_matches = sorted(mis_matches, key=lambda x: x[1])

    return mis_matches


def write_output(gene, gene_name, mis_matches, known_mutations, unknown_flag):
    """
    This function takes a gene name a list of mis matches found betreewn subject and query of
    this gene, the dictionary of known mutation in the point finder database, and the flag telling
    weather the user wants unknown mutations to be reported.
    All mis matches are looked up in the known mutation dict to se if the mutation is known,
    and in this case what drug resistence it causes.
    The funtions returns a 3 strings that are used as output to the users.
    One string is only tab seperated and contains the mutations listed line by line.
    If the unknown flag is set to true it will contain both known and unknown mutations.
    The next string contains only known mutation and are given in in a format that is easy to
    convert to HTML. The last string is the HTML tab sting from the unknown mutations.
    """
    RNA = False
    known_header = "Mutation\tNucleotide change\tAmino acid change\tResistance\tPMID\n"
    unknown_header = "Mutation\tNucleotide change\tAmino acid change\n"
    if gene in RNA_gene_list:
        RNA = True
        known_header = "Mutation\tNucleotide change\tResistance\tPMID\n"
        unknown_header = "Mutation\tNucleotide change\n"

    known_lst = []
    unknown_lst = []
    all_results_lst = []
    output_mut = []

    # Go through each mutation
    for i in range(len(mis_matches)):
        m_type = mis_matches[i][0]
        pos = mis_matches[i][1]  # sort on pos?
        look_up_pos = mis_matches[i][2]
        look_up_mut = mis_matches[i][3]
        mut_name = mis_matches[i][4]
        nuc_ref = mis_matches[i][5]
        nuc_alt = mis_matches[i][6]
        ref = mis_matches[i][-2]
        alt = mis_matches[i][-1]

        # First index in list indicates if mutation is known
        output_mut += [[]]
        output_mut[i] = [0]

        # Define output vaiables
        gene_mut_name = gene_name + " " + mut_name
        print(nuc_ref, nuc_alt)
        codon_change = nuc_ref + " -> " + nuc_alt
        aa_change = ref + " -> " + alt
        if RNA is True:
            aa_change = "RNA mutations"
        elif pos < 0:
            aa_change = "Promoter mutations"
        resistence = "Unknown"
        pmid = "-"

        # Check if mutation is known
        if look_up_pos in known_mutations[gene][m_type]:
            if look_up_mut in known_mutations[gene][m_type][look_up_pos]:
                # First index set to 1 indicates that mutation is known
                output_mut[i][0] = 1
                gene_mut_name = known_mutations[gene][m_type][look_up_pos][look_up_mut]['gene_name'] + " " + mut_name
                resistence = known_mutations[gene][m_type][look_up_pos][look_up_mut]['drug']
                pmid = known_mutations[gene][m_type][look_up_pos][look_up_mut]['pmid']

        output_mut[i] += [gene_mut_name, codon_change, aa_change, resistence,
                          pmid]

        # Add mutation to output strings for known mutations
        if output_mut[i][0] == 1:
            if RNA is True:
                # don't include the amino acid change field for RNA mutations
                known_lst += ["\t".join(output_mut[i][1:3]) + "\t" + "\t".join(output_mut[i][4:])]
            else:
                known_lst += ["\t".join(output_mut[i][1:])]
            all_results_lst += ["\t".join(output_mut[i][1:])]

        # Add mutation to output strings for unknown mutations
        else:
            if RNA is True:
                unknown_lst += ["\t".join(output_mut[i][1:3])]
            else:
                unknown_lst += ["\t".join(output_mut[i][1:4])]
            if unknown_flag is True:
                all_results_lst += ["\t".join(output_mut[i][1:])]

        # Check that you do not print two equal lines (can happen it to indels
        # occure in the same codon)
        if len(output_mut) > 1:
            if output_mut[i] == output_mut[i - 1]:

                if output_mut[i][0] == 1:
                    known_lst = known_lst[:-1]
                    all_results_lst = all_results_lst[:-1]
                else:
                    unknown_lst = unknown_lst[:-1]
                    if unknown_flag is True:
                        all_results_lst = all_results_lst[:-1]

    # Creat final strings
    all_results = "\n".join(all_results_lst)
    total_known_str = ""
    total_unknown_str = ""

    # Check if there are only unknown mutations
    known_index = [mut[0] for mut in output_mut]
    if sum(known_index) > 0:
        total_known_str = known_header + "\n".join(known_lst)
    else:
        total_known_str = "No known mutations found in %s" % gene_name

    if sum(known_index) < len(known_index):
        total_unknown_str = unknown_header + "\n".join(unknown_lst)
    else:
        total_unknown_str = "No unknown mutations found in %s" % gene_name

    return all_results, total_known_str, total_unknown_str


##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################


parser = ArgumentParser(description="This program...")

# positional arguments
parser.add_argument("-i", "--inputfiles", nargs='+', help="Input fasta file(s)")  # Don't konw if default should be sat to anything
parser.add_argument("-o", "--out_path", help="Path to blast output")  # er det runroot eller er det tmp directiory
parser.add_argument("-s", "--species", choices=['e.coli', 'gonorrhoeae', 'campylobacter', 'salmonella', 'tuberculosis', 'malaria'], help="Specie chosen for point mutation detetion")

# optional arguments
parser.add_argument("-p", "--databasePath", dest="db_path",help="Path to the databases", default='')  #/home/data1/services/ResFinder/database_archive/database-3.0_Test_indels
parser.add_argument("-b", "--blastPath", dest="blast_path",help="Path to blast", default='blastn')
parser.add_argument("-t", "--threshold", dest="threshold",help="Blast threshold for identity", default=0.75)
parser.add_argument("-l", "--min_cov", dest="min_cov",help="Minimum coverage", default=0.2)
parser.add_argument("-u", "--unknown_mut", dest="unknown_mutations", action="store_true", help="Show all mutations found even if in unknown to the resistance database", default=False)
parser.add_argument("-g", "--specific_genes", nargs='+', dest="specific_genes",help="Specifie genes existing in the database to search for only - if non is specified all genes are used", default=None)


args = parser.parse_args()

# Check if valid input file is provided
for infile in args.inputfiles:
    if not os.path.exists(infile):
        sys.exit("Input Error: Input file, {:s} does not exist!\n".format(infile))
inputfiles = args.inputfiles

# Check if valid output directory is provided
if not os.path.exists(args.out_path):
    # Try  to create the given output path
    out_path = os.path.abspath(args.out_path)
    try:
        os.makedirs(out_path)
    except:
        sys.exit("Input Error: Output dirctory does not exists!\n")
else:
    out_path = args.out_path

# Check if valid database path is provided
if not os.path.exists(args.db_path):
    sys.exit("Input Error: The specified database directory, %s, does not exist!\n" % args.db_path)
else:
    db_path = args.db_path

    # file configuartion
    # Check existence of config file
"""
    db_config_file = '%s/config'%(args.db_path)
    if not os.path.exists(db_config_file):
        sys.exit("Input Error: The database config file could not be "
                          "found!")
    # Check existence of notes file
    notes_path = "%s/notes.txt"%(args.db_path)
    if not os.path.exists(notes_path):
        sys.exit('Input Error: notes.txt not found! (%s)'%(notes_path))
    # Save path
    db_path = args.db_path
"""

# Save speices name and path to database
specie = args.species
specie_path = db_path + "/" + specie
# Check if the species database exist
if not os.path.exists(specie_path):
    sys.exit("database error: the provided database, {:s}, is not correctly installed").format(specie_path)
else:
    # Get names of genes that database consist of by reading the file genes.txt
    try:
        with open(specie_path + "/genes.txt", "r") as gene_file:
            gene_list = []
            for line in gene_file:
                gene_list.append(line.strip())
    except IOError:
        sys.exit("The gene list \"genes.txt\" doesn't exist in database: %s" % (specie_path))
    # Get names of RNA genes
    try:
        with open(specie_path + "/RNA_genes.txt", "r") as RNA_gene_file:
            RNA_gene_list = []
            for line in RNA_gene_file:
                RNA_gene_list.append(line.strip())
    except IOError:
        sys.exit("The gene list \"RNA_genes.txt\" doesn't exist in database: %s" % (specie_path))

# Check if valid blast path is provided
if not os.path.exists(args.blast_path):
    sys.exit("Input Error: The specified blast path, %s, does not exist!\n" % args.blast_path)
else:
    blast = args.blast_path

# Create user defined gene_list if applied
if args.specific_genes:

    # Check that the genes are valid
    genes_specified = []
    for gene in args.specific_genes:
        if gene in gene_list:
            genes_specified.append(gene)
        else:
            print("Input Error: Provided database was not recognised! (%s)\
                   \nYou can between the following genes:\n" % (gene))
            for gene in gene_list:
                print(gene)
            sys.exit()
    # Change the gene_list to the user defined gene_list
    gene_list = genes_specified

min_cov = args.min_cov
threshold = args.threshold

# Set a flag for specifing if unknown mutations are reported
unknown_flag = False
if args.unknown_mutations is True:
    unknown_flag = True


###############################################################################
# MAIN
###############################################################################


# Open resistens-overview file and extract mutation information
# known_mutations = get_known_muts(specie_path) #, databases)

known_mutations = get_db_mutations(specie_path + "/resistens-overview.txt",
                                   gene_list)


# Blast input(s) file against db
for inputfile in inputfiles:
    filename = inputfile.split("/")[-1]
    sample_name = ("").join(filename.split(".")[:-1])
    if sample_name == "":
        sample_name = filename
    # Calling blast and parsing output
    min_cov_blast = 0.01
    results, qry_align, homo_align, sbjct_align = Blaster(
        inputfile, gene_list, specie_path, out_path, min_cov_blast, threshold,
        blast)

    GENES = dict()

    # Find gene hit with the largest coverage
    for gene, hits in results.items():
        # save all hits in the list 'hits_found'
        hits_found = []

        GENES[gene] = dict()

        # Check for hits in blast results
        if type(hits) is dict:
            GENES[gene]['found'] = 'partially'

            # check coverage for each hit, if coverage is 100% save gene
            # directly else save hit info. if only one hit, save gene directly
            for hit in hits.keys():
                coverage = results[gene][hit]['coverage']

                # append tuble with subject start and end positions to the
                # list 'hits_found'
                hits_found += [(results[gene][hit]['sbjct_start'],results[gene][hit]['sbjct_end'], results[gene][hit]['sbjct_string'], results[gene][hit]['query_string'], coverage)]

                # if coverage is 100% change found to yes and stop loop
                if coverage == 1.0:
                    GENES[gene]['found'] = 'yes'

            # Sort positions found
            hits_found = sorted(hits_found, key=lambda x: x[0])

            # Find best hit by concatenating sequences if more hits exist
            final_sbjct, final_qry, all_start, all_end, alternative_overlaps = find_best_sequence(hits_found)

            # Find gene coverage if sequences has been concatinated
            if len(hits_found) > 1:
                no_call = final_sbjct.count("N")
                sbjct_len = results[gene][hit]['sbjct_length']
                coverage = ((all_end - all_start - no_call) / float(sbjct_len))

            # Save blast output of gene in GENES
            if coverage >= min_cov:
                GENES[gene]['coverage'] = coverage
                GENES[gene]['sbjct_string'] = final_sbjct
                GENES[gene]['query_string'] = final_qry
                GENES[gene]['sbjct_start'] = all_start
                GENES[gene]['sbjct_end'] = all_end
                # GENES[gene]['sbjct_len'] = sbjct_len
                GENES[gene]['alternative_overlaps'] = alternative_overlaps
                GENES[gene]['mis_matches'] = []

            else:
                # Gene not found above given coverage
                GENES[gene]['found'] = 'Gene found with coverage (%f) below  minimum coverage threshold: %s' %(coverage, min_cov)
        else:
            # gene not found!
            GENES[gene]['found'] = 'Gene not found'

    # Find known mutations and write output files
    # Output filenames
    output_files = ["PointFinder_results.txt", "PointFinder_table.txt"]

    # Initiate output stings with header
    output_strings = ["Mutation\tNucleotide change\tAmino acid change\tResistance\tPMID",
                      "Chromosomal point mutations - Results\nSpecies: %s\n\n\nKnown Mutations\n"%(specie)]
    total_unknown_str = ""

    # Find mutation in gene if gene is found
    for gene in GENES:
        # Start writing output string (to HTML tab file)
        gene_name = gene
        regex = r"promoter_size_(\d+)(?:bp)"
        promtr_gene_objt = re.search(regex, gene)
        if promtr_gene_objt:
            gene_name = gene.split("_")[0]
        output_strings[1] += "\n%s\n" % (gene_name)

        # Check if gene is found
        if GENES[gene]['found'] == 'yes' or GENES[gene]['found'] == 'partially':
            sbjct_start = GENES[gene]['sbjct_start']
            sbjct_seq = GENES[gene]['sbjct_string']
            qry_seq = GENES[gene]['query_string']
            alternative_overlaps = GENES[gene]['alternative_overlaps']
            # Find and save mis_matches in gene
            GENES[gene]['mis_matches'] = find_mismatches(gene, sbjct_start, sbjct_seq, qry_seq, alternative_overlaps)

        # If gene isn't found write the reason saved in GENES[gene]['found']
        if GENES[gene]['found'] != 'yes' and GENES[gene]['found'] != 'partially':
            output_strings[1] += GENES[gene]['found'] + "\n"
        else:
            # Check if any mutations was found
            if len(GENES[gene]['mis_matches']) < 1:
                output_strings[1] += "No mutations found in %s\n" % (gene_name,)
            else:
                # Write mutations found to output file
                total_unknown_str += "\n%s\n" % (gene_name)
                all_results, total_known, total_unknown = write_output(gene, gene_name, GENES[gene]['mis_matches'], known_mutations, unknown_flag)

                # Add results to output strings
                output_strings[0] += "\n" + all_results
                output_strings[1] += total_known + "\n"

                # Add unknown mutations the total results of unknown mutations
                total_unknown_str += total_unknown + "\n"

    # Add unknown results to all results
    if unknown_flag is True:
        output_strings[1] += "\n\nUnknown Mutations \n" + total_unknown_str

    # Write output files
    for i in range(len(output_files)):
        file_fullpath = out_path + "/" + output_files[i]
        with open(file_fullpath, "w") as file_:
            file_.write(output_strings[i])
    print(output_strings[1])
    print(output_strings[0])
