#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
#                         Script Description                                  #
# name:                                                                       #
#     lib_seqrutinator.py                                                     #
#                                                                             #
# description:                                                                #
#     Automatic protein superfamily analysis                                  #
#                                                                             #
# Authors:                                                                    #
#     Nicolas Stocchi, Agustin Amalfitano, Fernando Villareal, Marcelo        #
#     Atencio, Arjen ten Have                                                 #
###############################################################################


###############################################################################
#######################   AUTO-DATAMINING LIBRARY   ###########################
### Libraries #################################################################

import glob
import os
import re
import numpy as np
import statistics
from statsmodels import robust
import scipy.stats
from scipy.stats import shapiro
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import multiprocessing as mp
from subprocess import check_output
from datetime import datetime
startTime = datetime.now()

path = os.path.dirname(os.path.abspath(__file__))

###############################################################################
######################   Functions with general uses   ########################
###############################################################################

# Input files (complex function)

numbers = re.compile(r'(\d+)')


def numericalSort(value):

    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])

    return parts


def input_files(ftype):

    list_of_files = sorted((glob.glob(ftype)), key=numericalSort)

    return list_of_files

# Single List Reader


def list_reader(filename):
    single_list = []

    f = open(filename, 'r')
    for line in f:
        single_list.append(line.split()[0])
    f.close()
    return single_list

# Rename files


def rename_files(list_of_files, dst1):

    for i in list_of_files:
        os.system("cp " + i + " " + str(dst1))

    l = range(1, len(list_of_files) + 1)

    news = []

    for i, j in zip(list_of_files, l):

        new_name = "SF_" + str(j) + "_v1.fa"
        news.append(new_name)
        os.rename(i, new_name)

    # Generate file map
    f0 = open("file_map", "w")
    f0.write("Original name\tNew name\n")
    for i, j in zip(list_of_files, news):
        f0.write(str(i) + '\t' + str(j) + '\n')
    f0.close()

    os.system("mv file_map " + str(dst1))

    return()

# Rename a Mafft file


def rename_mafft(mafft_file):

    temp_mafft = str((mafft_file.split(".")[:-1])[0])

    new_temp = str(temp_mafft) + ".faa"
    os.system("mv Mafft_" + str(mafft_file) + " " + str(new_temp))

    return new_temp

# Arguments maker


def arg_maker(list_of_files, cores):

    list_of_arguments = []
    la_temp = []

    for group in list_of_files:

        if len(la_temp) == cores:

            list_of_arguments.append(la_temp)
            la_temp = []
            la_temp.append(group)

        else:
            la_temp.append(group)

        if group == list_of_files[-1] and la_temp not in list_of_arguments:
            list_of_arguments.append(la_temp)

    return list_of_arguments

# Fetching


def fetching(hits, target, output):

    with open(output, "w") as f:

        for seq_record in SeqIO.parse(target, "fasta"):

            if seq_record.id in hits:
                f.write(">" + str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")

    f.close()

    return()

# Anti-fetching


def anti_fetching(hits, target, output):

    with open(output, "w") as f:

        for seq_record in SeqIO.parse(target, "fasta"):

            if seq_record.id not in hits:
                f.write(">" + str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")

    f.close()

    return()

# Anti-fetching faa to fsa


def output_fsa_generator(bad_hits, input_file, good_output, bad_output):
    with open(good_output, 'w') as good_f, open(bad_output, 'w') as bad_f:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            if seq_record.id not in bad_hits:
                good_f.write(">" + str(seq_record.id) + "\n")
                string = ""
                for elem in seq_record.seq:
                    if elem != "-":
                        string += elem
                good_f.write(string + "\n")
            else:
                bad_f.write(">" + str(seq_record.id) + "\n")
                string = ""
                for elem in seq_record.seq:
                    if elem != "-":
                        string += elem
                bad_f.write(string + "\n")
    return()

# Fetching Parallel


def fetching_parallel(list_of_arguments):

    check_len = []

    # LOG #####################################################################
    print(list_of_arguments[0])
    print(list_of_arguments[1])
    print(list_of_arguments[2])
    ###########################################################################

    with open(str(path) + "/" + str(list_of_arguments[2]), "w") as f:

        for seq_record in SeqIO.parse(list_of_arguments[1], "fasta"):

            if seq_record.id in list_of_arguments[0]:
                f.write(">" + str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")
                check_len.append(seq_record.id)

                if len(check_len) == len(list_of_arguments[0]):
                    break

    f.close()

    # LOG #####################################################################
    print(str(list_of_arguments[2]) + " generated: " +
    str((datetime.now() - startTime)))
    ###########################################################################

    return()

# Multi-Fetching


def multi_fetc(list_of_hits, list_of_fastas, n_target):

    if isinstance(n_target, str):  # Check if it is a string or a list
        dataset = n_target

        for lg, lf in zip(list_of_hits, list_of_fastas):
            list_of_arguments = []

            for lgg, lfg in zip(lg, lf):
                com = [lgg, dataset, lfg]
                list_of_arguments.append(com)

            # Multiprocessing
            pool = mp.Pool(processes=len(list_of_arguments))
            pool.map(fetching_parallel, list_of_arguments)
            pool.close()  # IMPORTANT: Close the pool
            pool.join()

    else:

        list_of_targets = list(n_target)

        for lg, lt, lf in zip(list_of_hits, list_of_targets, list_of_fastas):
            list_of_arguments = []

            for lgg, ltg, lfg in zip(lg, lt, lf):
                com = [lgg, ltg, lfg]
                list_of_arguments.append(com)

            # Multiprocessing
            pool = mp.Pool(processes=len(list_of_arguments))
            pool.map(fetching_parallel, list_of_arguments)
            pool.close()  # IMPORTANT: Close the pool
            pool.join()

    return()

def number_seq(fasta_file):

    out = check_output("grep '>' " + str(fasta_file) + " | wc -l", shell=True)

    val = int(out[:-1])

    return val


# Sequences Extractor (REALLY IMPORTANT)


def seqs_extractor(fasta_file):

    names = []
    seqs = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        names.append(seq_record.id)
        seqs.append(seq_record.seq)

    return names, seqs

def seqs_extractor_desc(fasta_file):

    names = []
    seqs = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        names.append(seq_record.description)
        seqs.append(seq_record.seq)

    return names, seqs


# Scores Extractor


def scores_extractor(hmms):

    names = []
    scores = []

    fh = open(hmms, "r")
    linesfh = fh.readlines()

    for xh in linesfh:

        try:
            names.append(xh.split()[0])
            scores.append(float(xh.split()[5]))

        except IndexError:
            names.append(xh.split()[0])
            scores.append(float(xh.split()[1]))

    fh.close()

    return names, scores

# Scores Extractor by hmms


def scores_extractor_by_hmms(fasta, hmms):

    f_names, f_seqs = seqs_extractor(fasta)

    names = []
    scores = []

    fh = open(hmms, "r")
    linesfh = fh.readlines()

    for xh in linesfh:
        names.append(xh.split()[0])
        scores.append(float(xh.split()[5]))
    fh.close()

    real_n = []
    real_s = []

    for n, s in zip(names, scores):
        if n in f_names:
            real_n.append(n)
            real_s.append(s)

    return real_n, real_s

# Remove gaps


def remove_gaps(fasta_file):

    names, seqs = seqs_extractor(fasta_file)

    seqs_wo_gaps = []

    for seq in seqs:
        str_seq = ""

        for aa in seq:
            if aa != "-":
                str_seq += aa

        seqs_wo_gaps.append(str_seq)

    filename = str((fasta_file.split(".")[:-1])[0])
    output_file = str(filename) + ".fsa"

    f1 = open(output_file, "w")
    for name, seq in zip(names, seqs_wo_gaps):
        f1.write(">" + str(name) + "\n")
        f1.write(str(seq) + "\n")
    f1.close()

    return output_file

# Split big fasta


def split_big_fasta(fasta_file, cores):

    val = number_seq(fasta_file)
    print(val)

    res = val / cores

    m = int(res) + 1

    names = []
    seqs = []
    n = 1

    for seq_record in SeqIO.parse(fasta_file, "fasta"):

        if len(names) < m:

            #print(len(names))
            names.append(str(seq_record.id))
            seqs.append(str(seq_record.seq))

        if len(names) == m:

            # LOG #############################################################
            #print(len(names))
            print(("Split " + str(n) + " is finished " +
            str((datetime.now() - startTime))))
            ###################################################################

            f = open("split_" + str(n) + "_" + str(fasta_file), "w")
            for k, l in zip(names, seqs):
                f.write(">" + str(k) + "\n")
                f.write(str(l) + "\n")
            f.close()
            names = []
            seqs = []
            n = n + 1

    # Last group
    if names != []:
        f = open("split_" + str(n) + "_" + str(fasta_file), "w")
        for k, l in zip(names, seqs):
            f.write(">" + str(k) + "\n")
            f.write(str(l) + "\n")
        f.close()

        # LOG #################################################################
        print(("Split " + str(n) + " is finished " +
        str((datetime.now() - startTime))))
        #######################################################################

    list_of_splits = input_files("split_*")

    return list_of_splits

# Split a list based on the amount of ips


def split_ips_list(list_of_files, ips, num_ips):

    res = len(list_of_files) / len(ips)

    m = int(res) + 1

    list_for_ips = []
    lip = []

    for f_file in list_of_files:

        if len(lip) < m:
            lip.append(f_file)

        if len(lip) == m:
            list_for_ips.append(lip)
            lip = []

    # Last group
    if lip != []:
        list_for_ips.append(lip)

    ips_names_list = []

    for ip, nip in zip(list_for_ips, range(1, len(list_for_ips) + 1)):
        ips_names_list.append(str(nip) + "_ipfile_" + str(num_ips))

        fip = open(str(nip) + "_ipfile_" + str(num_ips), "w")
        for i in ip:
            fip.write(str(i) + "\n")
        fip.close()

    os.system("cat " + str(path) + "/*_ipfile_" + str(num_ips) + " > " + str(path) + "/total_ipfiles_" + str(num_ips))

    return list_for_ips, ips_names_list

# Sorting filelist (from high to low)


def sorting_filelist_HL(len_list, list_to_sort):

    sorted_list = [a for b, a in sorted(zip(len_list, list_to_sort), reverse=True)]

    return sorted_list

# Sorting filelist (from low to high)


def sorting_filelist_LH(len_list, list_to_sort):

    sorted_list = [a for b, a in sorted(zip(len_list, list_to_sort), reverse=False)]

    return sorted_list

# Old sorting filelist


def sorting_filelist(len_list, list_to_sort):

    sorted_list = [a for b, a in sorted(zip(len_list, list_to_sort), reverse=True)]

    return sorted_list

# Truncate


def truncate(flot):
    return float('%.1f'%(flot))

# Bad characters seq remover


def bad_char_detector(fasta_file, bad_char_file, clean_file):
    m = list('0123456789BJOUXZ?')

    names, seqs = seqs_extractor(fasta_file)
    bad_seqs = []
    for name, seq in zip(names, seqs):
        file_chars = list(seq)
        for fc in file_chars:
            if fc in m:
                print (("Invalid character detected in " + str(name)))
                bad_seqs.append(name)
                break
    fetching(bad_seqs, fasta_file, bad_char_file)
    anti_fetching(bad_seqs, fasta_file, clean_file)
    return()

###############################################################################
###############################   Processes   #################################
###############################################################################
# MAFFT common (auto-MAFFT)


def MAFFT(fasta_file):

    os.system("mafft " + str(fasta_file) +
    " > Mafft_" + str(fasta_file))
    return()

# MAFFT GINSI


def MAFFT_ginsi(fasta_file):

    os.system("mafft --reorder --maxiterate 1000 --retree 1 --globalpair " +
    str(fasta_file) + " > Mafft_" + str(fasta_file))

    return()

# MAFFT global


def MAFFT_global(fasta_file):

    os.system("mafft --globalpair " + str(fasta_file) + " > Mafft_" + str(fasta_file))

    return()

# MAFFT for SEQrutinator


def MAFFT_seqrutinator(fasta_file):

    names, seqs = seqs_extractor(fasta_file)

    if len(names) <= 500:
        os.system("mafft --reorder --maxiterate 1000 --retree 1 --globalpair " +
        str(fasta_file) + " > Mafft_" + str(fasta_file))

    else:
        MAFFT_global(fasta_file)

    return()

# MAFFT parallel


def MAFFT_parallel(list_of_arguments):

    os.system("mafft " + str(path) + "/" +
    str(list_of_arguments[0]) + " > " + str(path) + "/Mafft_" +
    str(list_of_arguments[0]))

    return()

# MAFFT-add


def MAFFT_add(to_add, pre_msa, out_msa, cores):

    ta_names, ta_seqs = seqs_extractor(to_add)

    if len(ta_names) > 1000:
        os.system("mafft --retree 1 --add " + str(to_add) + " " + str(pre_msa) +
        " > Mafft_" + str(out_msa))

    else:
        os.system("mafft --retree 2 --add " + str(to_add) + " " + str(pre_msa) +
        " > Mafft_" + str(out_msa))

    return()

# MAFFT-add parallel


def parallel_MAFFT_add(list_of_arguments):

    val = number_seq(list_of_arguments[0])

    if val > 1000:
        os.system("mafft " + " --retree 1 --add " +
        str(list_of_arguments[0]) + " Mafft_" + str(list_of_arguments[1]) +
        " > Mafft_" + str(list_of_arguments[2]))

    else:
        os.system("mafft " + " --retree 2 --add " +
        str(list_of_arguments[0]) + " Mafft_" + str(list_of_arguments[1]) +
        " > Mafft_" + str(list_of_arguments[2]))

    return()

# Multi-MAFFT-add


def multi_MAFFT_add(list_of_groups, list_to_add, list_nn_names):

    for lg, la, ln in zip(list_of_groups, list_to_add, list_nn_names):

        list_of_arguments = []

        for lgg, laa, lnn in zip(lg, la, ln):
            com = [laa, lgg, lnn]
            list_of_arguments.append(com)

        # Multiprocessing
        pool = mp.Pool(processes=len(list_of_arguments))
        pool.map(parallel_MAFFT_add, list_of_arguments)
        pool.close()  # IMPORTANT: Close the pool
        pool.join()

    return()

# FAMSA


def FAMSA(fasta_file):
    
    os.system("./famsa " + str(fasta_file) + " Mafft_" + str(fasta_file))
    # Mafft temporary flag is used to avoid to change almost every format file
    
    return()


# Check BMGE version
def check_bmge_version(bv_arg):
    # Execute the command to get BMGE version information
    output = "bmge_version.txt"
    os.system(f'java -jar BMGE.jar -h > {output}')
    
    bvers = open(output, "r")
    lines = bvers.readlines()
    
    if '\n' in lines:
        detected_version = 2 
    elif '   incorrect options\n' in lines:
        detected_version = 1
    else:
        print("Error: BMGE.jar file not found. Make sure this file is in the same directory as Seqrutinator")
        return False
    
    bvers.close()
    
    os.remove(output)
        # Check if the detected version matches the specified bv argument
    if bv_arg != detected_version:
        print(f"Error: BMGE version {detected_version} detected, but -bv argument is {bv_arg}.")
        return False

        # Version and bv argument match
    return True

# run BMGE
def BMGE(version, fasta, hh, BMGE_out):

    if version == 1: 
        os.system(f'java -jar BMGE.jar -i {fasta} + -t AA -h {hh}  > {BMGE_out}')

    if version == 2: 
        os.system(f'java -jar BMGE.jar -i {fasta} -t AA -o output_int_bmge -e {hh} > {BMGE_out}') 
        os.remove('output_int_bmge')

    return()

# BMGE reader
def BMGE_reader(version,BMGE_output):

    real_lines = []

    if version == 1:
        
        fb = open(BMGE_output, "r")
        lines = fb.readlines()

        for line, n_line in zip(lines, range(1, len(lines))):

            #print(n_line, line)

            if line not in real_lines:
                real_lines.append(str(line))

        fb.close()

        in_list = []
        out_list = []

        for rl in real_lines:
            if "before" in rl:
                in_seq = int(rl.split(" ")[5])
                out_char = int(rl.split(" ")[8])
                in_list.append(in_seq)
                out_list.append(out_char)

            if "after" in rl:
                list_rl = str(rl.split("\r")[0])
                in_seq = int(list_rl.split(" ")[6])
                out_char = int(list_rl.split(" ")[9])
                in_list.append(in_seq)
                out_list.append(out_char)

        alert = "NO"

        if in_list[0] != out_list[0]:
            alert = "YES"

        else:
            alert = "NO"

    if version == 2:
        
        fb = open(BMGE_output, "r")
        lines = fb.readlines()
        fb.close()

        in_list = []
        out_list = []

        for rl in lines:
            if "input" in rl:
                in_seq = int(rl.split(" ")[4])
                in_char = int(rl.split(" ")[6])
                in_list.append(in_seq)
                in_list.append(in_char)

            if "output" in rl:
                out_seq = int(rl.split(" ")[3])
                out_char = int(rl.split(" ")[5])
                out_list.append(out_seq)
                out_list.append(out_char)

        alert = "n.a."

    return alert, in_list[0], out_list[0], in_list[1], out_list[1]


# Hmmbuild


def hmmbuild(fasta_file):

    os.system("hmmbuild --informat afa --amino hmm_" + str(fasta_file) +
    " " + str(fasta_file))

    return()

# Hmmbuild parallel


def hmmbuild_parallel(list_of_arguments):

    os.system("hmmbuild --informat afa --amino " + str(path) +
    "/hmm_" + str(list_of_arguments[0]) + " " + str(path) + "/Mafft_" +
    str(list_of_arguments[0]))

    return()

# Hmmsearch_common


def hmmsearch_common(input_file, target_file):
    real_input = input_file.split(".")[0]
    output_file = "output_" + str(real_input)

    os.system("hmmsearch " + str(input_file) + " " + str(target_file) +
    " > " + str(output_file))

    return output_file

# Hmmsearch-tab (everything above HMMER inclusion threshold)


def hmmsearch_tab_hmm(input_file, output_file, target_file):

    os.system("hmmsearch --noali --tblout hmms_" + str(output_file) +
    " " + str(input_file) + " " + str(target_file))

    f = open("hmms_" + str(output_file), 'r')
    lines = f.readlines()
    f.close()

    w = open("hmms_" + str(output_file), 'w')
    for line in lines:
        if not '#' in line:
            w.write(line)
    w.close()

    return()

# Hmmsearch-tab (everything above HMMER inclusion threshold)


def hmmsearch_tab(input_file, output_file, target_file):

    os.system("hmmsearch --noali --tblout hmms_" + str(output_file) +
    " hmm_" + str(input_file) + " " + str(target_file))

    f = open("hmms_" + str(output_file), 'r')
    lines = f.readlines()
    f.close()

    w = open("hmms_" + str(output_file), 'w')
    for line in lines:
        if not '#' in line:
            w.write(line)
    w.close()

    return()

# Hmmsearch parallel


def hmmsearch_parallel(list_of_arguments):

    os.system("hmmsearch --noali --tblout " + str(path) + "/hmms_" +
    str(list_of_arguments[1]) + " " + str(path) + "/hmm_" +
    str(list_of_arguments[0]) + " " + str(list_of_arguments[2]))

    f = open(str(path) + "/hmms_" + str(list_of_arguments[1]), 'r')
    lines = f.readlines()
    f.close()

    w = open(str(path) + "/hmms_" + str(list_of_arguments[1]), 'w')
    for line in lines:
        if not '#' in line:
            w.write(line)
    w.close()

    return()

# Hmmsearch parallel without uncomment function to avoid problems with
# multi-hmmsearch


def hmmsearch_parallel_2(list_of_arguments):

    os.system("hmmsearch --noali --tblout " + str(path) + "/hmms_" +
    str(list_of_arguments[1]) + " " + str(path) + "/hmm_" +
    str(list_of_arguments[0]) + " " + str(list_of_arguments[2]))

    return()

# Multi-Hmmsearch


def multi_hmmsearch(list_of_splits, group, cores):

    list_of_arguments = []

    for ls, ln in zip(list_of_splits, range(1, len(list_of_splits) + 1)):
        com = [group, "tesp_" + str(ln) + "_" + str(group), ls]
        # tesp means TEmporary SPlit
        list_of_arguments.append(com)

    # Multiprocessing
    pool = mp.Pool(processes=cores)
    pool.map(hmmsearch_parallel_2, list_of_arguments)
    pool.close()  # IMPORTANT: Close the pool
    pool.join()

    # Concatenate hmmsearches
    new_hmms = "hmms_" + str(group)

    # WATCH OUT with cat hmms_tesp_*, it is necessary to specify the group
    os.system("cat " + str(path) + "/hmms_tesp_*_" + str(group) + " > " + str(path) + "/" + str(new_hmms))

    n_seqs = []
    scores = []

    fh = open(str(path) + "/" + str(new_hmms), "r")
    linesfh = fh.readlines()
    for xh in linesfh:
        if '#' not in xh:
            n_seqs.append(xh.split()[0])
            # SOLVE THIS
            scores.append(float(xh.split()[5]))
    fh.close()

    s = sorted(zip(scores, n_seqs), reverse=True)
    scores, n_seqs = map(list, zip(*s))

    fn = open(str(path) + "/" + str(new_hmms), "w")  # WATCH OUT
    for a, b in zip(n_seqs, scores):
        fn.write(str(a) + "\t" + str(b) + "\n")
    fn.close()

    return new_hmms


# Multi-processing


def multi_pro(list_of_arguments, func, dataset):

    if func == 'M':
        function = MAFFT_parallel

    if func == 'HB':
        function = hmmbuild_parallel

    if func == 'HS':
        function = hmmsearch_parallel

        if isinstance(dataset, str):  # Check if it is a string or a list

            for arg in list_of_arguments:
                new_arg = []

                for group in arg:
                    com = [group, group, dataset]
                    new_arg.append(com)

                pool = mp.Pool(processes=len(new_arg))
                pool.map(function, new_arg)
                pool.close()  # IMPORTANT: Close the pool
                pool.join()

        else:

            for arg, dat in zip(list_of_arguments, dataset):
                new_arg = []

                for group, n_data in zip(arg, dat):
                    com = [group, group, n_data]
                    new_arg.append(com)

                pool = mp.Pool(processes=len(new_arg))
                pool.map(function, new_arg)
                pool.close()  # IMPORTANT: Close the pool
                pool.join()

    if func == 'M' or func == 'HB':

        for arg in list_of_arguments:
            arg2 = []

            for a2 in arg:
                arg2.append([a2])

            # Multiprocessing
            pool = mp.Pool(processes=len(arg2))
            pool.map(function, arg2)
            pool.close()  # IMPORTANT: Close the pool
            pool.join()

    return()

def hits_generator(fasta, hits_list):

    hits_file = "hits_" + str(fasta)

    fo = open(hits_file, "w")
    for h in hits_list:
        fo.write(str(h) + "\n")
    fo.close()

    return hits_file

# Remove gappy sequences

def remove_gappy_seqs(target, filename, gap_len, output_file):

    gap_substring_len =  "-" * gap_len

    names, seqs = seqs_extractor(filename)

    selected_seqs = []
    gappy_seqs = []

    for n, s in zip(names, seqs):

        if gap_substring_len not in s:
            selected_seqs.append(n)

        else:
            gappy_seqs.append(n)

    fetching(selected_seqs, target, output_file)

    return gappy_seqs

# Remove big sequences


def big_seqs_remover(target, filename, seq_len, output_file):
    gap_list = []

    names, seqs = seqs_extractor(filename)

    for name, seq in zip(names, seqs):
        gap_i_list = [x for x in seq if x == "-"]
        len_gap_i = len(gap_i_list)
        gap_list.append(len_gap_i)

    average_gap = np.mean(gap_list)
    std = np.std(gap_list)

    threshold = average_gap - seq_len * std

    selected_seqs = []
    big_seqs = []

    for name, seq, g_list in zip(names, seqs, gap_list):

        if g_list > threshold:
            selected_seqs.append(name)

        else:
            big_seqs.append(name)

    fetching(selected_seqs, target, output_file)

    return big_seqs

# Shapiro-Wilk Test


def shapiro_test(data):
    # normality test
    stat, p = shapiro(data)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
        print('Sample looks Gaussian')
        print ('(fail to reject H0)')
        normal = True
    else:
        print('Sample does not look Gaussian (reject H0)')
        normal = False
    return(normal)

# z_scores


def f_zscores(list1):
    z_scores2 = list(scipy.stats.zscore(list1))
    # sort the list (not necessary I guess)
    z_scores2.sort(reverse = True)
    return z_scores2

###############################################################################
############################   MAIN MODULES AD   ##############################
###############################################################################
# Short Sequences Remover (SSR) ###############################################


def SSR(fasta_file, real_fasta, ref_seq, ref_len, ref_perc, ali):

    new_name = "new_name"

    # Extract names and sequences from MSA
    names, seqs = seqs_extractor(fasta_file)

    # Identify if it is a faa or fsa
    wog = 0

    for seq in seqs:

        if "-" in seq:
            wog = wog + 1
            break

    if wog == 0:

        # Count the amount of amino for each sequence
        amino_seqs = []

        for seq in seqs:
            amino_seq = [x for x in seq if x != "-"]
            amino_seqs.append(len(amino_seq))

        # Threshold calculation
        threshold = ref_len * (ref_perc / 100)

        # Check for short sequences
        remove_names = []
        remove_seqs = []
        remove_len = []

        good_names = []
        good_seqs = []

        for name, seq, len_i in zip(names, seqs, amino_seqs):

            if len_i < threshold:
                remove_names.append(name)
                remove_seqs.append(seq)
                remove_len.append(len_i)

            else:
                if name != ref_seq:  # Extract reference sequence
                    good_names.append(name)
                    good_seqs.append(seq)
                else:
                    remove_names.append(name)
                    remove_seqs.append(seq)
                    remove_len.append(len_i)

        if remove_names != []:

            new_name = "SSR_" + str(real_fasta) + ".fsa"

            f1 = open(new_name, "w")
            for name_g, seq_g in zip(good_names, good_seqs):
                f1.write(">" + str(name_g) + "\n")
                f1.write(str(seq_g) + "\n")
            f1.close()

            # Sort short sequences
            s = sorted(zip(remove_len, remove_names))
            remove_len, remove_names = map(list, zip(*s))

            f2 = open("SSR_" + str(real_fasta) + "_data.txt", "w")
            f2.write("Threshold: " + str(threshold) + "\n")
            for seq_i, len_i in zip(remove_names, remove_len):
                f2.write(str(seq_i) + "\t" + str(len_i) + "\n")
            f2.close()

            f3 = open("SSR_" + str(real_fasta) + "_removed.fsa", "w")
            for name_i, seq_i in zip(remove_names, remove_seqs):
                f3.write(">" + str(name_i) + "\n")
                f3.write(str(seq_i) + "\n")
            f3.close()

            f4 = open("SSR_" + str(real_fasta) + "_removed.txt", "w")
            for seq_i in (remove_names):
                f4.write(str(seq_i) + "\n")
            f4.close()

            if ali == 1:
                MAFFT_ginsi(new_name)
            if ali == 2:
                MAFFT_seqrutinator(new_name)
            if ali == 3:
                FAMSA(new_name)

    if wog == 1:

        wog_fasta = remove_gaps(fasta_file)
        names, seqs = seqs_extractor(wog_fasta)

        # Count the amount of amino for each sequence
        amino_seqs = []
        for seq in seqs:
            amino_seq = [x for x in seq if x != "-"]
            amino_seqs.append(len(amino_seq))

        # Threshold calculation
        threshold = ref_len * (ref_perc / 100)

        # Check for short sequences
        remove_names = []
        remove_seqs = []
        remove_len = []

        good_names = []
        good_seqs = []

        for name, seq, len_i in zip(names, seqs, amino_seqs):
            if len_i < threshold:
                remove_names.append(name)
                remove_seqs.append(seq)
                remove_len.append(len_i)
            else:
                if name != ref_seq:  # Extract reference sequence
                    good_names.append(name)
                    good_seqs.append(seq)
                else:
                    remove_names.append(name)
                    remove_seqs.append(seq)
                    remove_len.append(len_i)

        if remove_names != []:
            new_name = "SSR_" + str(real_fasta) + ".fsa"
            f1 = open(new_name, "w")
            for name_g, seq_g in zip(good_names, good_seqs):
                f1.write(">" + str(name_g) + "\n")
                f1.write(str(seq_g) + "\n")
            f1.close()

            # Sort short sequences
            s = sorted(zip(remove_len, remove_names))
            remove_len, remove_names = map(list, zip(*s))

            f2 = open("SSR_" + str(real_fasta) + "_data.txt", "w")
            f2.write("Threshold: " + str(threshold) + "\n")
            for seq_i, len_i in zip(remove_names, remove_len):
                f2.write(str(seq_i) + "\t" + str(len_i) + "\n")
            f2.close()

            f3 = open("SSR_" + str(real_fasta) + "_removed.fsa", "w")
            for name_i, seq_i in zip(remove_names, remove_seqs):
                f3.write(">" + str(name_i) + "\n")
                f3.write(str(seq_i) + "\n")
            f3.close()

            f4 = open("SSR_" + str(real_fasta) + "_removed.txt", "w")
            for seq_i in remove_names:
                f4.write(str(seq_i) + "\n")
            f4.close()

            if ali == 1:
                MAFFT_ginsi(new_name)
            if ali == 2:
                MAFFT_seqrutinator(new_name)
            if ali == 3:
                FAMSA(new_name)

    # Check if I need the fsa or Mafft_fsa
    return new_name, remove_names


###############################################################################
# Non-homologues hits remover (NHHR)


def NHHR(fasta_file, real_fasta, alpha, thr, ref_len, ref_perc, ali):

    names, seqs = seqs_extractor(fasta_file)

    new_fsa = "NHHR_" + str(real_fasta) + ".fsa"
    removed_seqs = []

    if ref_perc != 0:
        # HMMER ###############################################################

        # Remove gaps
        wog_fasta = remove_gaps(fasta_file)
        os.system("cp " + str(wog_fasta) + " NHHR_" + str(real_fasta) + ".fsa")

        # Hmmbuild
        hmmbuild(fasta_file)

        # Hmmsearch
        hmmsearch_tab(fasta_file, fasta_file, wog_fasta)
        os.system("mv hmm_" + str(fasta_file) + " NHHR_" + str(real_fasta) + ".hmm")

        H_names, H_scores = scores_extractor_by_hmms(fasta_file, "hmms_" + str(fasta_file))
        os.system("mv hmms_" + str(fasta_file) + " NHHR_" + str(real_fasta) + ".txt")

        # Statistics
        # Std
        std = np.std(H_scores)

        # Average
        average = np.mean(H_scores)

        # Threshold
        threshold = average - alpha * std

        # Median
        median = np.median(H_scores)

        # Quartiles
        Q3 = np.percentile(H_scores, 75)
        Q1 = np.percentile(H_scores, 25)
        IQR = Q3 - Q1

        # Set the statistical threshold
        if thr == 1:
            threshold = average - alpha * std
            method = "mean - " + str(alpha) + " * SD"

        if thr == 2:
            threshold = Q1 - 1.5 * IQR
            method = "IQR"

        # LOG #################################################################
        print ("\nStatistical Analysis")
        print ("Average: "), average
        print ("Std: "), std
        print ("Threshold: "), threshold
        print ("Median: "), median
        print ("IQR: "), IQR
        #######################################################################

        removed_seqs = []
        removed_score2 = []
        removed_len2 = []

        really_bad_seqs = [x for x in names if x not in H_names]

        if really_bad_seqs:
            for elem in really_bad_seqs:
                removed_seqs.append(elem)

        # Remove sequences with low score and medium or high length
        # Count the amount of amino for each sequence

        amino_seqs = []

        minimum_ref_len = ref_len * (ref_perc / 100)

        for seq in seqs:
            amino_seq = [x for x in seq if x != "-"]
            amino_seqs.append(len(amino_seq))

        # New threshold pondering by length

        for name, score in zip(H_names, H_scores):
            # Check low score

            if score < threshold:

                # Check length
                for am_seq, name2 in zip(amino_seqs, names):

                    if name == name2:

                        # If medium/high length
                        if am_seq >= minimum_ref_len:
                            removed_seqs.append(name)
                            removed_score2.append(score)
                            removed_len2.append(am_seq)
                            break

                        else:
                            break

        # Function of output files
        new_fsa = "NHHR_" + str(real_fasta) + ".fsa"
        seq_re = "NHHR_" + str(real_fasta) + "_seqs_removed.txt"
        bad_output = "NHHR_" + str(real_fasta) + "_removed.fsa"
        output_fsa_generator(removed_seqs, fasta_file, new_fsa, bad_output)
        seq_lengths = "NHHR_" + str(real_fasta) + "_seqs_lengths.txt"
        bad_output_txt = "NHHR_" + str(real_fasta) + "_removed.txt"
        histog = "NHHR_histogram_" + str(real_fasta) + ".png"

        # Histogram ###########################################################
        plt.hist(H_scores, bins='auto')
        plt.title("HMMER scores distribution")
        plt.axvline(threshold, color='red', linestyle='dashed', linewidth=1, label='Threshold: ' + str(round(threshold,2)))
        plt.legend(loc='upper right')
        plt.xlabel("HMMER score")
        plt.savefig(histog)
        plt.close()
        #######################################################################

        # Save seqs lengths output file, sort the lists in zip style
        s = sorted(zip(amino_seqs, names), reverse=False)
        len3, names3 = map(list, zip(*s))

        with open(seq_lengths, "w") as f1:
            for am_len3, name3 in zip(len3, names3):
                perce = am_len3 / ref_len
                f1.write(str(name3) + " " + str(am_len3) + "/" + str(perce) + "\n")
        f1.close()

        with open(bad_output_txt, "w") as f2:
            for seq in removed_seqs:
                f2.write(str(seq) + "\n")
        f2.close()
        #######################################################################

        # Add a parameter to do mafft
        if removed_seqs:
            if ali == 1:
                MAFFT_ginsi(new_fsa)
            if ali == 2:
                MAFFT_seqrutinator(new_fsa)
            if ali == 3:
                FAMSA(new_fsa)

        # Remove "really bad seq" in "remove seqs" to save correctly output
        removed_seqs_2 = [x for x in removed_seqs if x not in really_bad_seqs]

        # Save sequences, score and length
        with open(seq_re, "w") as f3:
            f3.write("Threshold: " + str(threshold) + "\n")
            f3.write("Threshold Method: " + str(method) + "\n")
            f3.write("Minimum length: " + str(ref_perc / 100) + "*" + str(ref_len) + "= " + str(minimum_ref_len) + "\n\n")
            f3.write("Std: " + str(std) + "\n")
            f3.write("Average: " + str(average) + "\n")
            f3.write("Median: " + str(median) + "\n")
            f3.write("IQR: " + str(IQR) + "\n")

            if removed_seqs_2:
                f3.write("\nSequence (HMMER score/ length)\n")
                for name, score, length in zip(removed_seqs_2, removed_score2, removed_len2):
                    f3.write(str(name) + ": (" + str(score) + "/ "+ str(length) + ")\n")

            else:
                f3.write("There are no sequences TOO bad with score > 0 (low score and medium length)")

        # Save the really really bad sequences
            if really_bad_seqs:
                f3.write("\nReally Bad Sequences (hmm score = 0)\n")
                for seq in really_bad_seqs:
                    f3.write(str(seq) + "\n")
        f3.close()

    return new_fsa, removed_seqs

###############################################################################
# Gap Instigator Identifier (GIR)


def gap_instigator_identifier(fasta_file, thr, aa, gir_plt, prob_seq_file):
    # THE MOST IMPORTANT LIST:
    real_names = []

    # Extract names and sequences from MSA
    names, seqs = seqs_extractor(fasta_file)
    #print(seqs)
    #for seq in seqs:
    #    print(seq)

    # Lists
    l = len(seqs)  # Amount of sequences
    len_c = len(seqs[0])  # Amount of columns
    range_len_c = range(0, len_c)
    len_cols = list(range(1, len(seqs[0]) + 1))  # Length of columns

    cs = 0  # Cumulative score
    gs = 0  # Global score
    limg = 1  # Gap limit non-occupancy duplicator

    scores_col = []  # Score columns +1, 0
    scores_col2 = []  # Score columns +1, -1 * 2 progressive

    blocks_ps = []  # Possible problematic sequences by blocks of gaps
    blocks_sc = []  # Scores of blocks
    blocks_gs = []  # Scores of blocks 2
    all_globals = []  # All global scores detected
    all_globals_st_end = []  # # Start and end of the global block [start, end]
    blocks_st_end = []  # Start and end of the block [start, end]
    pos_seqs = []  # Possible problematic sequences

    # Calculate column score

    for lc in range_len_c:
        n = 0
        for name, seq in zip(names, seqs):
            if seq[lc] == "-":
                #print(seq[lc])
                n = n + 1
                #print(n)
        s = (n / l) * 100

        if s >= thr:
            for name, seq in zip(names, seqs):
                if seq[lc] != "-":
                    pos_seqs.append(name)
            cs = cs + 1
            gs = gs + 1
            limg = 1
            scores_col.append(cs)
            scores_col2.append(gs)
            if lc == range_len_c[-1] and cs >= aa:
                blocks_ps.append(pos_seqs)
                blocks_sc.append(cs)
                blocks_st_end.append([lc - cs + 1, lc])
                if gs >= aa:
                    blocks_gs.append(gs)
            if lc == range_len_c[-1] and cs < aa:
                if gs >= aa:
                    blocks_gs.append(gs)
                    blocks_st_end.append([lc - gs + 1, lc])
        else:
            if cs >= aa:
                blocks_ps.append(pos_seqs)
                blocks_sc.append(cs)
                blocks_st_end.append([lc - cs + 1, lc])
                if gs >= aa:
                    blocks_gs.append(gs)
            if gs >= aa:
                all_globals.append(gs)
                all_globals_st_end.append([lc - gs + 1, lc])
            pos_seqs = []  # INDENTATION
            if cs != 0:  # If gap block ended
                if gs > 0:
                    gs = gs - limg
            if cs == 0:
                if gs > 0:
                    limg = limg * 2
                    gs = gs - limg
            if gs < 0:
                gs = 0
            cs = 0
            scores_col.append(cs)
            scores_col2.append(gs)
        #print(str(lc) + " P: " + str(s) + "% ready")

    total_score = sum(blocks_sc)

    # Find possible sequences

    # Problematic sequences: names, sum scores and block size
    total_pr_names = []  # Filtered duplicates prob_names
    prob_names = []
    prob_scores = []
    prob_global_scores = []
    prob_block = []
    prob_block_scores = []

    lines_for_f1 = []

    plt.figure()
    plt.title('GIR Gap Scores (GGSs) in MSA columns')
    plt.ylabel('GGSs')
    plt.xlabel('MSA columns')
    plt.plot(len_cols, scores_col, color='blue')
    plt.plot(len_cols, scores_col2, linestyle='dotted', color='green')
    labels_leg = ['Continuous GGS', 'Combined GGS']
    plt.legend(labels_leg, bbox_to_anchor=(0, 0.925, 0.75, .075), ncol=2,mode='expand',borderaxespad=0.)
    if max(scores_col) >= max(scores_col2):
        var = int(max(scores_col) / aa) + 2
    if max(scores_col2) > max(scores_col):
        var = int(max(scores_col2) / aa) + 2
    for i in range(var):
        line_plot = i * aa
        #print(line_plot)
        plt.plot(np.linspace(0, len(len_cols), num=500), np.linspace(line_plot,line_plot, num=500), linestyle='dotted', color='black', alpha=0.5)
    #plt.show()
    plt.savefig(gir_plt)
    plt.close()    

    range_blocks = range(1, len(blocks_ps) + 1)
    sort_blocks = sorted(zip(blocks_sc, blocks_ps, blocks_gs, blocks_st_end, range_blocks), reverse=True)
    try:
        blocks_sc, blocks_ps, blocks_gs, blocks_st_end, range_blocks = map(list, zip(*sort_blocks))

        for blp, bls, blg, blc, bll in zip(blocks_ps, blocks_sc, blocks_gs, blocks_st_end, range_blocks):
            #print blp
            blp_fil = []
            occur = []
            for b in blp:
                if b not in blp_fil:
                    blp_fil.append(b)

            for bb in blp_fil:
                ocr = [1 for x in blp if x == bb]
                occur.append(sum(ocr))


            s = sorted(zip(occur, blp_fil), reverse=True)
            occur, blp_fil = map(list, zip(*s))

            for bf, oc in zip(blp_fil, occur):
                if oc >= aa:
                    if bf not in total_pr_names:
                        total_pr_names.append(bf)
                    prob_names.append(bf)

                    prob_scores.append(oc)
                    prob_block.append(bll)
                    glc = blg - bls + oc  # Global by column
                    prob_global_scores.append(glc)
                    prob_block_scores.append(bls)
        
        # Based on possible sequences, find problematic sequences

        f1 = open(prob_seq_file, "w")
        if lines_for_f1 != []:
            f1.write("Important:\n")
            for f1line in lines_for_f1:
                f1.write(f1line)

        if total_pr_names != []:
            f1.write("\nMax possible problematic combined score: " + str(total_score) + "\n")
            f1.write("Sequence\tContinuous_Score\tCombined_Score\tInstigated_blocks: Number of block(block len): Cont. block score/Comb. block score\n")

        total_scores = []
        global_scores = []
        local_scores_bl = []

        for pr_seq in total_pr_names:
            pr_score = []
            pr_glob = []
            pr_blocks = []
            pr_bl_score = []
            for p1, p2, p3, p4, p5 in zip(prob_names, prob_scores, prob_block, prob_global_scores, prob_block_scores):
                if p1 == pr_seq:
                    pr_score.append(p2)
                    pr_blocks.append(p3)
                    pr_glob.append(p4)
                    pr_bl_score.append(p5)

            ts = sum(pr_score)
            gs = sum(pr_glob)

            total_scores.append(ts)
            global_scores.append(gs)
            real_names.append(pr_seq)

            block_info = []

            for b1, b2, b3, b4 in zip(pr_blocks, pr_bl_score, pr_score, pr_glob):
                bi = str(b1) + "(" + str(b2) + "): " + str(b3) + "(Cont.)/" + str(b4) + "(Comb.)"
                block_info.append(bi)

            str_blocks = '\t'.join(str(e) for e in block_info)
            local_scores_bl.append(str_blocks)
        s = sorted(zip(total_scores, global_scores, local_scores_bl, real_names), reverse=True)
        total_scores, global_scores, local_scores_bl, real_names = map(list, zip(*s))

        for r1, r2, r3, r4 in zip(real_names, total_scores, global_scores, local_scores_bl):
            f1.write(str(r1) + "\t" + str(r2) + "\t" + str(r3) + "\t" + str(r4) + "\n")
        f1.close()

    except ValueError:
        f3 = open(prob_seq_file, "w")
        f3.write("Your alignment is good\n")
        if all_globals != [] and blocks_sc == []:
            f3.write("But just in case, check this out:\n")
            range_blocks = range(1, len(all_globals) + 1)
            sort_blocks = sorted(zip(all_globals, all_globals_st_end, range_blocks), reverse=True)
            all_globals, all_globals_st_end, range_blocks = map(list, zip(*sort_blocks))

            for blg, blc, bll in zip(all_globals, all_globals_st_end, range_blocks):
                f3.write('Gap block ' + str(bll) + " (" + str(blc[0]) + " - " + str(blc[1]) + ") has a global score: " + str(blg) + "\n")
        f3.close()
    return real_names

# Gap Instigator Remover


def GIR(met, fasta_file, real_fasta, thr, aa, ali):

    new_name = "new_name"

    prob_lines = []

    it = 0  # Iteration as variable

    GIR_plot = "GIR_it" + str(it + 1) + "_Plot_" + str(real_fasta) + ".png"
    GIR_prob_seq = "GIR_it" + str(it + 1) + "_Data_" + str(real_fasta) + ".txt"
    GIR_rem_fsa = "GIR_" + str(real_fasta) + "_removed.fsa"
    GIR_rem = "GIR_" + str(real_fasta) + "_removed.txt"

    if met == 1:  # One by one
        real_hits = []

        # First GIR
        real_names = gap_instigator_identifier(fasta_file, thr, aa, GIR_plot, GIR_prob_seq)
        #print(real_names)

        # Remove gaps from fasta file
        # wog_fasta: without gaps
        wog_fasta = remove_gaps(fasta_file)

        while real_names != []:
            # Iteration number
            it = it + 1

            GIR_plot = "GIR_it" + str(it + 1) + "_Plot_" + str(real_fasta) + ".png"
            GIR_prob_seq = "GIR_it" + str(it +1) + "_Data_" + str(real_fasta) + ".txt"

            if it == 1:
                # Anti-fetching of fasta file
                hit = real_names[0]
#                print(hit)
                new_name = "GIR_it" + str(it) + "_" + str(real_fasta) + ".fsa"
                anti_fetching(hit, wog_fasta, new_name)
                real_hits.append(real_names[0])

            else:
                # Anti-fetching of fasta file
                hit = real_names[0]
#                print(hit)
                new_name = "GIR_it" + str(it) + "_" + str(real_fasta) + ".fsa"
                anti_fetching(hit, "GIR_it" + str(it - 1) + "_" + str(real_fasta) + ".fsa", new_name)
                real_hits.append(real_names[0])

            # MAFFT of new fasta
            if ali == 1:
                MAFFT_ginsi(new_name)
            if ali == 2:
                MAFFT_seqrutinator(new_name)
            if ali == 3:
                FAMSA(new_name)

            fasta_file = rename_mafft(new_name)

            # GIR of new aligned fasta
            partial_names = gap_instigator_identifier(fasta_file, thr, aa, GIR_plot, GIR_prob_seq)
            real_names = list(partial_names)

        fetching(real_hits, wog_fasta, GIR_rem_fsa)

        f1 = open(GIR_rem, "w")
        for rh in real_hits:
            f1.write(str(rh) + "\n")
        f1.close()

        f2 = open(GIR_prob_seq, "r")
        for line in f2:
            if line != "Your alignment is good\n" and line != "But just in case, check this out:\n":
                prob_lines.append(line)
        f2.close()

    if met == 2:  # Batch

        real_hits = []

        # First GIR
        real_names = gap_instigator_identifier(fasta_file, thr, aa, GIR_plot, GIR_prob_seq)
        #print(real_names)
        # Remove gaps from fasta file
        wog_fasta = remove_gaps(fasta_file)

        while real_names != []:
            # Iteration number
            it = it + 1

            GIR_plot = "GIR_it" + str(it + 1) + "_Plot_" + str(real_fasta) + ".png"
            GIR_prob_seq = "GIR_it" + str(it +1) + "_Data_" + str(real_fasta) + ".txt"

            if it == 1:
                # Anti-fetching of fasta file
                new_name = "GIR_it" + str(it) + "_" + str(real_fasta) + ".fsa"
                anti_fetching(real_names, wog_fasta, new_name)

                for rn in real_names:
                    real_hits.append(rn)
            else:
                # Anti-fetching of fasta file
                new_name = "GIR_it" + str(it) + "_" + str(real_fasta) + ".fsa"
                anti_fetching(real_names, "GIR_it" + str(it - 1) + "_" + str(real_fasta), new_name)

                for rn in real_names:
                    real_hits.append(rn)

            # MAFFT of new fasta
            if ali == 1:
                MAFFT_ginsi(new_name)
            if ali == 2:
                MAFFT_seqrutinator(new_name)
            if ali == 3:
                FAMSA(new_name)    

            # GIR of new aligned fasta
            partial_names = gap_instigator_identifier("Mafft_" + str(new_name), thr, aa, GIR_plot, GIR_prob_seq)
            real_names = list(partial_names)

        fetching(real_hits, wog_fasta, GIR_rem_fsa)

        f1 = open(GIR_rem, "w")
        for rh in real_hits:
            f1.write(str(rh) + "\n")
        f1.close()

        f2 = open(GIR_prob_seq, "r")
        for line in f2:
            if line != "Your alignment is good\n" and line != "But just in case, check this out:\n":
                prob_lines.append(line)
        f2.close()

    return new_name, real_hits, prob_lines

###############################################################################
# Continuous Gap Sequence Remover


def CGSR(fasta_file, real_fasta, thr, minimum, ali):

    final_removed_seqs = []
    all_data = []
    removed_seqs_i = []

    new_name = "new_name"
    iteration = 1
    thr = thr / 100

    wog_fasta = remove_gaps(fasta_file)

    while True:

        if iteration > 1:

            new_name = "CGSR_It" + str(iteration) + "_" + str(real_fasta) + ".fsa"
            anti_fetching(removed_seqs_i, wog_fasta, new_name)

            wog_fasta = new_name

            if ali == 1:
                MAFFT_ginsi(new_name)
            if ali == 2:
                MAFFT_seqrutinator(new_name)
            if ali == 3:
                FAMSA(new_name)    

            fasta_file = rename_mafft(new_name)

        #print(fasta_file)
        removed_seqs_i = []

        names, seqs = seqs_extractor(fasta_file)
        #print(names)
        #print(seqs)

        l = len(seqs)  # Amount of sequences
        #print(l)
        len_c = len(seqs[0])  # Amount of columns

        gap_column = []

        for lc in range(0, len_c):
            n = 0

            for name, seq in zip(names, seqs):
                if seq[lc] == "-":
                    n = n + 1

            # s counts the amount of gaps in column i
            s = (n / l)

            if s <= thr:
                gap_column.append(0)

            else:
                gap_column.append(1)

        # We want the total amount of gaps in amino columns, total_unusual_len have
        # this numbers, and total_unusual_name have the sequence name
        total_unusual_name = []
        total_unusual_len = []
        partial_gaps = []

        for name, seq in zip(names, seqs):
            n = 0
            N = 0

            #print(name)
            #print(seq)

            unusual_rowgap_lengh = []

            for lc, gap_col in zip(range(0, len_c), gap_column):

                if gap_col == 0:

                    # Check if seq[i] is a gap
                    if seq[lc] == "-":
                        n += 1
                        N += 1

                    else:
                        if n > 0:
                            unusual_rowgap_lengh.append(n)
                            # Restart partial sum
                            n = 0

            # Check the final part of the sequence
            if n > 0:
                unusual_rowgap_lengh.append(n)

            # Sort unusual_rowgap_lengh
            s_i = sorted(unusual_rowgap_lengh, reverse=True)

            if len(s_i) > 10:
                s_i = [s_i[j] for j in range(0, 10)]

            partial_gaps.append(s_i)
            total_unusual_name.append(name)
            total_unusual_len.append(N)

        # Sort the lists in zip style
        s = sorted(zip(total_unusual_len, total_unusual_name), reverse=True)
        total_unusual_len, total_unusual_name = map(list, zip(*s))

        #print(total_unusual_len)
        #print(total_unusual_name)
        max_local_gap_subtring = []

        with open("CGSR_It" + str(iteration) + "_" + str(real_fasta) + "_Gapdata.txt", "w") as f2:

            for num, name in zip(total_unusual_len, total_unusual_name):

                f2.write(str(name) + " " + str(num) + "\n")

                for ss_i, name2 in zip(partial_gaps, names):

                    if name == name2:

                        if ss_i:
                            max_ss_i = max(ss_i)
                            max_local_gap_subtring.append(max_ss_i)

                        else:
                            max_local_gap_subtring.append(0)

                        for elem in ss_i:
                            f2.write(str(elem) + " ")
                        f2.write("\n")

                        break

        s = sorted(zip(max_local_gap_subtring, total_unusual_name), reverse=True)
        max_local_gap_subtring, total_unusual_name = map(list, zip(*s))

        # Statistical analysis
        distribution_elems = []

        #print(partial_gaps)

        for elem in partial_gaps:

            #print(elem)
            t_i = [x for x in elem if x > 10]

            #print(t_i)

            if t_i:
                for elem2 in t_i:
                    distribution_elems.append(elem2)

        #print(distribution_elems)

        if len(distribution_elems) <= 1:
            all_data.append("WARNING: YOUR ALIGNMENT IS TOO GOOD. THERE ARE NO GAPS > 10 WITH AMINO COLUMNS")
            break

        # MAD
        mad = robust.scale.mad(distribution_elems,
            c=0.6744897501960817, axis=0)
        #print "MAD", mad

        # Median
        median = np.median(distribution_elems)

        # Average gaps
        average_gap = np.mean(distribution_elems)

        # Std
        std = np.std(distribution_elems)


        #ithreshold = average_gap + 2 * std
        #mad_threshold = median + mad
        Q3 = np.percentile(distribution_elems, 75)
        Q1 = np.percentile(distribution_elems, 25)
        IQR = Q3 - Q1
        threshold = Q3 + 1.5 * IQR
        # Histogram
        plt.hist(distribution_elems, bins='auto')
        plt.title("Gap size distribution")
        plt.axvline(threshold, color='k', linestyle='dashed', linewidth=1, label= 'Threshold: ' + str(threshold))
        plt.xlabel('Continuous gap size')
        plt.legend(loc='upper right')
        plt.savefig('CGSR_It' + str(iteration) + "_" + str(real_fasta) + '_Gap_histo.png')
        plt.close()


        data_list = []

        data_list.append('Iteration ' + str(iteration) + "\n")
        data_list.append("Threshold: " + str(threshold) + "\n")
        data_list.append("Q1: " + str(Q1) + "\n")
        data_list.append("Q3: " + str(Q3) + "\n")
        data_list.append("IQR: " + str(IQR) + "\n")
        data_list.append("Mean: " + str(average_gap) + "\n")
        data_list.append("Median: " + str(median) + "\n")
        data_list.append("Mad: " + str(mad) + "\n")
        data_list.append("Std: " + str(std) + "\n")
        data_list.append("Candidate Sequences:\n")

        for local_sc, name in zip(max_local_gap_subtring, total_unusual_name):

            if local_sc >= max([threshold, minimum]):

                # Final removed list
                final_removed_seqs.append(name)

                # Removed list iteration i
                removed_seqs_i.append(name)
                data_list.append(str(name) + " " + str(local_sc) + "\n")
                #print name, local_sc

            else:
                break

        if max(max_local_gap_subtring) < threshold:

            thr2 = max([Q3, minimum])

            if max(max_local_gap_subtring) >= thr2:
                for local_sc, name in zip(max_local_gap_subtring, total_unusual_name):

                        if local_sc >= thr2:
                            removed_seqs_i.append(name)
                            final_removed_seqs.append(name)
                            data_list.append(str(name) + " " + str(local_sc) + "\n")

                        else:
                            break

        all_data.append(data_list)

        # stop while
        if not removed_seqs_i:
            break

        iteration += 1

    if final_removed_seqs:
        with open("CGSR_" + str(real_fasta) + "_removed.txt", "w") as f:
            for name in final_removed_seqs:
                f.write(str(name) + "\n")
        f.close()

    if iteration == 1:
        new_name = fasta_file

    return all_data, new_name, final_removed_seqs

###############################################################################
# Pseudogene Remover (PR)


def PR(filename, s5, auto, ali):
    new_filename = filename.split(".")

    # Read msa and extract sequences and names

    names, seqs = seqs_extractor(filename)
    

    # Remove gaps
    wog_fasta = remove_gaps(filename)

    if filename != str(new_filename[0] + ".faa"):
        if ali == 1:
            MAFFT_ginsi(filename)
        if ali == 2:
            MAFFT_seqrutinator(filename)
        if ali == 3:
            FAMSA(filename)    

        filename = rename_mafft(filename)

    # Hmmbuild
    hmmbuild(filename)
    os.system("cp hmm_" + str(filename) + " PR_" + str(new_filename[0]) + "_it1.hmm")
    # Hmmsearch
    hmmsearch_tab(filename, filename, wog_fasta)
    os.system("cp hmms_" + str(filename) + " PR_" + str(new_filename[0]) + "_it1.txt")
    # Save original scores and names
    H_names, H_scores = scores_extractor_by_hmms(wog_fasta, "hmms_" + str(filename))
    H_names_original = list(H_names)

    # Check if any sequence have 0 hmmer score and exit #######################
    stop_program = [x for x in names if x not in H_names]
    cero_score = [x for x in H_scores if x == 0]

    if stop_program != [] or cero_score != []:
        print("There are sequences with no hmmer score, please check your data")
        print(stop_program)

        anti_fetching(stop_program + cero_score, wog_fasta, filename)
        fetching(stop_program + cero_score, wog_fasta, "PR_NO_scores_sequences.fasta")

        if ali == 1:
            MAFFT_ginsi(filename)
        if ali == 2:
            MAFFT_seqrutinator(filename)
        if ali == 3:
            FAMSA(filename)    

        filename = rename_mafft(filename)

        # Hmmbuild
        hmmbuild(filename)

        # Hmmsearch
        hmmsearch_tab(filename, filename, wog_fasta)

        # Save original scores and names
        H_names, H_scores = scores_extractor_by_hmms(wog_fasta, "hmms_" + str(filename))
        H_names_original = list(H_names)

    ###########################################################################
    seqs_removed = []
    i = 0

    while True:

        print("Iteration " + str(i))
        shap_test = shapiro_test(H_scores)

        # Check normality
        if shap_test:
            print("In this iteration the distribution is normal")
            
        else:
            print("In this iteration the distribution is NOT normal")
            
        # Statistics
        std = statistics.pstdev(H_scores)
        average = statistics.mean(H_scores)
        Q3 = np.percentile(H_scores, 75)
        Q1 = np.percentile(H_scores, 25)
        IQR = Q3 - Q1

        # Compute scoredrop
        # Copy H_scores
        H_scores2 = list(H_scores[1:])
        H_scores1 = list(H_scores[:-1])

        # Remove first element
        diferential = [a - b for a, b in zip(H_scores1, H_scores2)]
        score_drop = max(diferential)

        # Save the score drop threshold
        for a, b in zip(H_scores1, H_scores2):
            if a - b == score_drop:
                score_drop_threshold = a
                break

        # Different thresholds
        m2sd = average - 2 * std
        m3sd = average - auto * std  # auto variable used
        IQR_threshold = Q1 - 1.5 * IQR

        # Count the amount of seqs below each threshold
        b_score_drop_threshold = len([x for x in H_scores if x < score_drop_threshold])

        b_m2sd = len([x for x in H_scores if x < m2sd])
        b_m3sd = len([x for x in H_scores if x < m3sd])
        b_IQR_threshold = len([x for x in H_scores if x < IQR_threshold])

        plt.scatter(range(len(H_scores)), H_scores, marker='o')
        plt.title("HMMER Score Distribution Plot")
        plt.axhline(y=score_drop_threshold, color='red', linestyle='dashed', linewidth=1,
        label='$\Delta$ max: ' + str(score_drop_threshold) + " (" + str(b_score_drop_threshold) + " seqs below)")
        plt.axhline(y=m2sd, color='blue', linestyle='dashed', linewidth=1,
             label='M-2sd: ' + str(truncate(m2sd)) + " (" + str(b_m2sd) + " seqs below)")
        plt.axhline(y=m3sd, color='green', linestyle='dashed', linewidth=1,
        label='M-' + str(auto) + 'sd: ' + str(truncate(m3sd)) + " (" + str(b_m3sd) + " seqs below)")
        diferential5 = [5 * x for x in diferential]
        plt.plot(range(len(diferential)), diferential5)
        plt.axhline(y=IQR_threshold, color='black', linestyle='dashed', linewidth=1,
        label='Lower bound: ' + str(truncate(IQR_threshold)) + " (" + str(b_IQR_threshold) + " seqs below)")

        if shap_test:
            normal_label = "This distribution is normal"
        else:
            normal_label = "this distribution is NOT normal"

        plt.plot([], [], ' ', label=normal_label)
        plt.ylabel("HMMER score")
        plt.xlabel("sequences")
        plt.legend()

        h_s_p = "PR_" + str(new_filename[0]) + "_it" + str(i +1) + ".png"
        plt.savefig(h_s_p)
        plt.close()

        
        if s5 == 1:
            PR_thresh = m3sd
        if s5 == 2:
            PR_thresh =  IQR_threshold  
        
        fetch_threshold = len([x for x in H_scores if x < PR_thresh])
        if fetch_threshold == 0:
            #print ("mean - 3sd threshold have not sequence below it, break while")
            break
        else:
            terminate = []
            H_names_reloaded = []
            for name, score in zip(H_names, H_scores):
                # check low score
                if score < PR_thresh:
                        terminate.append(name)
                else:
                    H_names_reloaded.append(name)
        #######################################################################

        new_group = list(H_names_reloaded)

        # Use temporal files for recompute MSAs, profiles, and searches
        name_it = "temporal_" + str(new_filename[0]) + ".fsa"
        seqs_removed = ["not_empty"]

        # fetching (parameters: target = "wog_fasta", list = "new_group")
        fetching(new_group, wog_fasta, name_it)

        # Mafft
        if ali == 1:
            MAFFT_ginsi(name_it)
        if ali == 2:
            MAFFT_seqrutinator(name_it)
        if ali == 3:
            FAMSA(name_it)

        new_f_name = rename_mafft(name_it)

        #Hmmbuild
        hmmbuild(new_f_name)
        os.system("cp hmm_" + str(new_f_name) + " PR_" + str(new_filename[0]) + "_it" + str(i +2) + ".hmm")
        #Hmmsearch
        hmmsearch_tab(new_f_name, new_f_name, wog_fasta)
        os.system("cp hmms_" + str(new_f_name) + " PR_" + str(new_filename[0]) + "_it" + str(i +2) + ".txt")
        
        # Rename the files go back into the loop
        # Check new worst sequence
        Group_H_names = []
        Group_H_scores = []

        with open("hmms_" + str(new_f_name), "r") as hf:
            for line in hf:
                seq_i_name = str(line.split()[0])
                if seq_i_name in new_group:
                    Group_H_names.append(seq_i_name)
                    Group_H_scores.append(float(line.split()[5]))

        # Rewrite the variables to do the next iteration
        H_scores = list(Group_H_scores)
        H_names = list(Group_H_names)
        i = i + 1

    # Outputs #################################################################
    PR_remain = "PR_" + str(new_filename[0]) + "_remain_seqs.txt"
    PR_removed_file = "PR_" + str(new_filename[0]) + "_removed.txt"
    PR_removed = []

    if seqs_removed:
        rrr = [x for x in H_names_original if x not in H_names]

        with open(PR_removed_file, "w") as f:
            for seq in rrr:
                f.write(str(seq) + "\n")
                PR_removed.append(seq)

        with open(PR_remain, "w") as f:
            for seq in H_names:
                f.write(str(seq) + "\n")
    ###########################################################################
    # Remove temporal files

    if seqs_removed:
        final_fsa = "PR_" + str(new_filename[0]) + ".fsa"
        final_faa = "PR_" + str(new_filename[0]) + ".faa"

        os.system("mv " + str(name_it) + " " + str(final_fsa))
        os.system("mv " + str(new_f_name) + " " + str(final_faa))
    else:
        final_fsa = "NOT.fsa"
        final_faa = "NOT.faa"
    ###########################################################################

    return PR_removed, PR_removed_file, PR_remain, final_fsa, final_faa
