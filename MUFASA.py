#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:07:48 2019

@author: Nicolas Stocchi
"""

import re
import os
import glob
import argparse
import logging
from Bio import SeqIO
from datetime import datetime
startTime = datetime.now()

path = os.path.dirname(os.path.abspath(__file__))

# Command line arguments ######################################################
parser = argparse.ArgumentParser()
parser.add_argument('-i', default='hmm_profile', help='HMMER profile, default = hmm_profile')
parser.add_argument('-ext', default='*.fsa', help='Target files extension. default = *.fsa')
parser.add_argument('-c', default=4, help='Cores, default = 4')

args = vars(parser.parse_args())

hmm_profile = str(args['i'])
ext = str(args['ext'])
cores = int(args['c'])

### Logging ###################################################################
logger = logging.getLogger('MUFASA')
hdlr = logging.FileHandler('MUFASA.log', 'w')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.DEBUG)

# LOG #########################################################################
logger.info('MUFASA initiated - Total time: ' +
str((datetime.now() - startTime)))
print('MUFASA initiated - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################
# FUNCTIONS ###################################################################
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

# Fetching


def fetching(hits, target, output):
    with open(output, "w") as f:
        for seq_record in SeqIO.parse(target, "fasta"):
            if seq_record.id in hits:
                f.write(">" + str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")
    f.close()
    return()

# Hmmsearch_common


def hmmsearch_common(input_file, target_file):
    os.system("hmmsearch " + str(input_file) + " " + target_file + 
    " > output.txt")

# Hmmsearch-tab (everything above HMMER inclusion threshold)


def hmmsearch_tab(input_file, real_fasta, target_file):
    os.system("hmmsearch --noali --tblout " + str(real_fasta) + "_hmms.txt" +
    " " + str(input_file) + " " + str(target_file))
    f = open(str(real_fasta) + "_hmms.txt", 'r')
    lines = f.readlines()
    f.close()
    w = open(str(real_fasta) + "_hmms.txt", 'w')
    for line in lines:
        if not '#' in line:
            w.write(line)
    w.close()
    return()

# MAFFT GINSI cores


def MAFFT_ginsi_c(fasta_file, cores):
    os.system("mafft --thread " + str(cores) + " --threadit " + 
    str(cores) + " --reorder --maxiterate 1000 --retree 1 --globalpair " + 
    str(fasta_file) + " > " + str(real_fasta) + "_hits.faa")
    return()


###############################################################################

list_of_fastas = input_files(ext)

# LOG #########################################################################
logger.info('List of fasta loaded - Total time: ' +
str((datetime.now() - startTime)))
print('List of fasta loaded - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################

list_of_hmms = []

os.system("mkdir Results_MUFASA")
dst = path + "/Results_MUFASA"
os.system("mkdir Results_MUFASA/hmmsearch")
hst = path + "/Results_MUFASA/hmmsearch"
os.system("mkdir Results_MUFASA/MAFFT")
mst = path + "/Results_MUFASA/MAFFT"
os.system("mkdir Results_MUFASA/Hits")
hht = path + "/Results_MUFASA/Hits"

f1 = open("Hits_summary.tsv", "w")
f1.write("Fasta\tTotal_Seq\tSeq_Detected\n")

for fasta in list_of_fastas:

    real_fasta = fasta.split('.')[0]

    # LOG #####################################################################
    logger.info('Hmmsearch using ' + str(hmm_profile) + ' with ' + str(fasta) + ' initiated' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    print('Hmmsearch using ' + str(hmm_profile) + ' with ' + str(fasta) + ' initiated' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

    
    hmmsearch_common(hmm_profile, fasta)

    o1 = open("output.txt", "r")
    for line in o1:
        if "Initial search space (Z):" in line:
            total_seqs = [int(s) for s in line.split() if s.isdigit()]
        if "Domain search space  (domZ):" in line:
            seqs_det = [int(s) for s in line.split() if s.isdigit()]
    o1.close()

    f1.write(str(fasta) + "\t" + str(total_seqs[0]) + "\t" + str(seqs_det[0]) + "\n")

    hmmsearch_tab(hmm_profile, real_fasta, fasta)

    hits = list_reader(str(real_fasta) + "_hmms.txt")

    fetching(hits, fasta, str(real_fasta) + "_hits.fsa")

    MAFFT_ginsi_c(str(real_fasta) + "_hits.fsa", cores)


    # LOG #####################################################################
    logger.info('Hmmsearch using ' + str(hmm_profile) + ' with ' + str(fasta) + ' finished' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    print('Hmmsearch using ' + str(hmm_profile) + ' with ' + str(fasta) + ' finished' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

f1.close()

os.system("rm output.txt")
os.system("mv Hits_summary.tsv " + str(dst))
os.system("mv *hits.fsa " + str(hht))
os.system("mv *.faa " + str(mst))
os.system("mv *hmms.txt " + str(hst))
os.system("mv *.log " + str(dst))

# LOG #########################################################################
logger.info('MUFASA finished - Total time: ' +
str((datetime.now() - startTime)))
print('MUFASA finished - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################
