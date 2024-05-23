#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:07:48 2019

@author: Nicolas Stocchi

Updated on Wed Apr 03 13:29:19 2024 

@by: Fernando Villarreal
"""

import re
import os
import glob
import argparse
import logging
from Bio import SeqIO
from datetime import datetime
import sys
startTime = datetime.now()

path = os.path.dirname(os.path.abspath(__file__))

# Command line arguments ######################################################
parser = argparse.ArgumentParser()
parser.add_argument('-i', default='hmm_profile', help='HMMER profile, default = hmm_profile')
parser.add_argument('-ext', default='*.fsa', help='Target files extension. default = *.fsa')
parser.add_argument('-a', default='y', help='Align hits with MAFFT? Use y (default) or n')
parser.add_argument('-c', default=4, help='Cores, default = 4')
parser.add_argument('-r', default='', help='Reference sequence for mafft-add and trimming (must be in fasta format)')
parser.add_argument('-t', default='5-400', help='Trimming N- and C-end based on reference sequence (residues corresponding to the provided indexes are retained)')

args = vars(parser.parse_args())

hmm_profile = str(args['i'])
ext = str(args['ext'])
align = str(args['a'])
cores = int(args['c'])
refseq = str(args['r'])
trimming = (args['t'].split('-'))
trimming = [int(x) for x in trimming]
n_end = int(trimming[0])
c_end = int(trimming[1]) 

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
def anti_fetching(hits, target, output):

    with open(output, "w") as f:

        for seq_record in SeqIO.parse(target, "fasta"):

            if seq_record.id not in hits:
                f.write(">" + str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")

    f.close()

    return()

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
def MAFFT_ginsi_c(fasta_file, real_fasta, cores):
    os.system(f"mafft --thread {cores} --threadit {cores} --reorder --maxiterate 1000 --retree 1 --globalpair {fasta_file} >  {real_fasta}_hits.faa")
    return()

def MAFFT_add(to_add, fasta_file, cores):

    real_fasta = "_".join(fasta_file.split('_')[:-1])

    os.system(f'mafft --thread {cores} --threadit {cores} --retree 2 --ep 0.0 --add {to_add} --globalpair --maxiterate 16 {fasta_file} > {real_fasta}_add.faa')

    return()

def process_sequences(fasta_file, start_pos, end_pos):
    # Extract sequences from the FASTA file
    names, seqs = seqs_extractor(fasta_file)
 
    # Get the reference sequence
    ref_seq = seqs[-1]

    # Initialize variables to track the adjusted start and end positions
    adj_start_pos = 0
    adj_end_pos = 0

    # Adjust start_pos based on characters in the reference sequence
    for i, char in enumerate(ref_seq):
        if char != "-":
            adj_start_pos += 1
        if adj_start_pos == start_pos:
            adj_start_pos = i + 1
            break

    # Adjust end_pos based on characters in the reference sequence
    for i, char in enumerate(ref_seq):
        if char != "-":
            adj_end_pos += 1
        if adj_end_pos == end_pos:
            adj_end_pos = i + 1
            break

    # Remove columns before adj_start_pos and after adj_end_pos in all sequences
    processed_seqs = []
    for seq in seqs:
        processed_seq = seq[adj_start_pos - 1:adj_end_pos]
        processed_seqs.append(processed_seq)

    return names, processed_seqs

if align != 'y' and refseq != '':
    print(f'If you want to perform MAFFT-add and reference trimming, -a must be set to "y"')
    exit()

list_of_fastas = input_files(ext)

#### Pre-chekup
#### Sequences IDs with spaces are not compatible with hmmer, so check

for fasta in list_of_fastas:
    names, seqs = seqs_extractor_desc(fasta)
    for id, seq in zip(names,seqs):
        if ' ' in str(id):
            # LOG #####################################################
            logger.critical('One or more IDs in the sequences of your fasta files have spaces. This is problematic for hmmer. Please remove spaces in the IDs or consider to rename sequences altogether (recommended)')
            print('One or more IDs in the sequences of your fasta files have spaces. This is problematic for hmmer')
            print('Please remove spaces in the IDs or consider to rename sequences altogether (recommended)')        
            sys.exit()
            ###########################################################
if refseq != '':
    refseq_seq = []
    refseq_name = []
    for seq_record in SeqIO.parse(refseq, "fasta"):
        refseq_name.append(seq_record.id)
        refseq_seq.append(seq_record.seq)
    if (len(refseq_name)) != 1:
        raise ValueError(f'Error: the file called with -r must have only one sequence. Please check')
    if c_end > len(refseq_seq[0]):
        raise IndexError(f'Error with the entered C-end trimming value of {c_end}. The reference sequence {refseq_name[0]} has {len(refseq_seq[0])} residues. Please check')

# LOG #########################################################################
logger.info('List of fasta loaded - Total time: ' +
str((datetime.now() - startTime)))
print('List of fasta loaded - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################

profile_name = hmm_profile.split(".")[0]

dst = path + "/Results_MUFASA_" + str(profile_name)
os.system("mkdir " + str(dst))
hst = dst + "/hmmsearch"
os.system("mkdir " + str(hst))
hht = dst + "/Hits"
os.system("mkdir " + str(hht))
if align == 'y':
    mst = dst + "/MAFFT"
    os.system("mkdir " + str(mst))

if refseq != '':
    ast = dst + "/MAFFT_add"
    os.system("mkdir " + str(ast))
    tst = dst + "/reference_trimmed"
    os.system("mkdir " + str(tst))

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
    if align == 'y':
        MAFFT_ginsi_c(str(real_fasta) + "_hits.fsa", real_fasta, cores)


    # LOG #####################################################################
    logger.info('Hmmsearch using ' + str(hmm_profile) + ' with ' + str(fasta) + ' finished' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    print('Hmmsearch using ' + str(hmm_profile) + ' with ' + str(fasta) + ' finished' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

f1.close()

if refseq != '':

    logger.info(f'Block of MSA trimming based on reference sequence iniciated. Total time {(datetime.now() - startTime)}' )

    refseq_name = []

    for seq_record in SeqIO.parse(refseq, "fasta"):
            refseq_name.append(seq_record.id)
        
    # Mafft add the reference to each fasta
    
    list_of_fastas_add = input_files('*_hits.faa')
    
    for fasta in list_of_fastas_add:
        MAFFT_add(refseq, fasta, cores)
        # LOG #####################################################################
        logger.info('Mafft add of sequence ' + str(refseq_name[0]) + ' to ' + str(fasta) + ' finished' + 
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('Mafft add of sequence ' + str(refseq_name[0]) + ' to ' + str(fasta) + ' finished' + 
        ' - Total time: ' + str((datetime.now() - startTime)))
        ###########################################################################

    # Trim N- and C-end based on the coordinates
    list_of_fasta_trim = input_files('*_add.faa') 
    
    for fasta in list_of_fasta_trim:

        names, processed_seqs = process_sequences(fasta, n_end, c_end)
     
        real_fasta = f"_".join(fasta.split('_')[:-1])
        output_file = str(real_fasta) + '_trimmed.faa'
        # write the output file, starting with the reference sequence
        with open(output_file, "w") as f:
            f.write(f">{names[-1]}\n")
            f.write(f"{processed_seqs[-1]}\n")
            for name, seq in zip(names[:-1], processed_seqs[:-1]):
                f.write(f">{name}\n")
                f.write(f"{seq}\n")
 
        # LOG #####################################################################
        logger.info('N- and C-end trimming based on ' + str(refseq_name[0]) + '\'s positions ' + str(trimming[0]) + '-' + str(trimming[1]) + ' for ' + str(fasta) + ' finished' + 
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('N- and C-end trimming based on ' + str(refseq_name[0]) + '\'s positions ' + str(trimming[0]) + '-' + str(trimming[1]) + ' for ' + str(fasta) + ' finished' +  
        ' - Total time: ' + str((datetime.now() - startTime)))
        ###########################################################################

os.system("rm output.txt")
os.system("mv Hits_summary.tsv " + str(dst))
os.system("mv *hits.fsa " + str(hht))
os.system("mv *hmms.txt " + str(hst))
os.system("mv *.log " + str(dst))
if refseq != '':
    os.system("mv *_add.faa " + str(ast))
    os.system("mv *_trimmed.faa " + str(tst))
if align == 'y':
    os.system("mv *.faa " + str(mst))

# LOG #########################################################################
logger.info('MUFASA finished - Total time: ' +
str((datetime.now() - startTime)))
print('MUFASA finished - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################
