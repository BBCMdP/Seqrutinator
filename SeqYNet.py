#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:14:28 2019

@authors: Nicolas Stocchi, Agustin Amalfitano, Fernando Villarreal, Marcelo
Atencio, Arjen ten Have
Final revision 03/25/2024 FV
"""

import os
import sys
import logging
import argparse
import lib_seqrutinator as lib
from datetime import datetime
startTime = datetime.now()

path = os.path.dirname(os.path.abspath(__file__))

###############################################################################
########   SeqYNet (Performs Seqrutinator on several datasets)     ############
###############################################################################
# Command Line Arguments
parser = argparse.ArgumentParser()

# General parameters
parser.add_argument('-ext', default='*.fasta', help='Extension of fasta files, default = \'*.fasta\'')
parser.add_argument('-s', default=1, help='Sort filelist from high to low (1) or from low to high (2)')
parser.add_argument('-m', default='12345', help='Pipelines')
parser.add_argument('-ali', default=1, help='Use MAFFT G-INS-i (1); Use either MAFFT G-INS-i (n<=500) or Global (n>500) (2); Use FAMSA, recommended only for really large datasets (3)')
parser.add_argument('-ref1', default='first_seq', help='First sequence in the MSA is the reference sequence (by default) or user\'s input reference sequence')
parser.add_argument('-ref2', default=0, help='User\'s input reference sequence length')
parser.add_argument('-bv', default=1, help='BMGE version 1 (either 1.0 or 1.12) or 2')
parser.add_argument('-BMGE', default=0, help='BMGE deactivated (0), BMGE > 0 is h option for BMGE')


# Module 1: SSR
parser.add_argument('-p1', default=0.65, help='Proportion of sequence length coverage for SSR (0 to  1)')

# Module 2: NHHR
parser.add_argument('-p2', default=0.65, help='Proportion of sequence length coverage for NHHR (0 to  1)')
parser.add_argument('-m2', default=1, help='HMMER Score (1) or TCS (2) method for NHHR ')
parser.add_argument('-s2', default=2, help="Mean - alphaSD (1) or Q1 - 1.5IQR (2)")
parser.add_argument('-a2', default=3, help='Alpha for NHHR')

# Module 3: GIR
parser.add_argument('-m3', default=1, help='Method one by one (1) or batch (2) for GIR')
parser.add_argument('-p3', default=0.9, help='Proportion of gaps to define a gap column for GIR (>= VALUE, from 0 to 1)')
parser.add_argument('-aa3', default=30, help='aa window of contiguos gap columns for GIR')

# Module 4: CGSR
parser.add_argument('-p4', default=0.5, help='Proportion of gaps to define a gap column for CGSR (>= VALUE, from 0 to 1)')
parser.add_argument('-aa4', default=30, help='aa window of contiguos gap columns for CGSR')

# Module 5: OR
parser.add_argument('-a5', default=3, help='Alpha for mean - alphaSD (3 is recommended as befault option and 2.35 for normal distributions)')
parser.add_argument('-s5', default=2, help='Mean - alphaSD (1) or Q1 - 1.5IQR (2)')

args = vars(parser.parse_args())

# General Parameters
ext = str(args['ext'])
sort_option = int(args['s'])
pipeline = str(args['m'])  # Terminal Condition
ali = int(args['ali'])
ref_seq = str(args['ref1'])
ref_len = int(args['ref2'])
BMGE = float(args['BMGE'])
bmge_version = int(args['bv'])

# Module 1: SSR
p1 = float(args['p1'])

# Module 2: NHHR
p2 = float(args['p2'])
m2 = int(args['m2'])
s2 = int(args['s2'])  # t2
a2 = float(args['a2'])  # alpha_1

# Module 3: GIR
m3 = int(args['m3'])
p3 = float(args['p3'])
aa3 = int(args['aa3'])

# Module 4: CGSR
p4 = float(args['p4'])
aa4 = int(args['aa4'])  # t4

# Module 5: OR
a5 = float(args['a5'])
s5 = float(args['s5'])

### LOG #######################################################################
logger = logging.getLogger('SeqYNet')
hdlr = logging.FileHandler('SeqYNet.log', 'w')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.DEBUG)

# LOG #########################################################################
logger.info('SeqYNet initiated - Total time: ' +
str((datetime.now() - startTime)))
print('SeqYNet initiated - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################

pre_fasta_list = []

filelist = lib.input_files(ext)

for ff in filelist:
    if "hmms_" not in ff and "Mafft_" not in ff:
        pre_fasta_list.append(ff)

command_parameters = []

# General Commands
if pipeline != '12345':
    command_parameters.append("-m " + str(pipeline))

if ali != 0:
    command_parameters.append("-ali " + str(ali))

if ref_seq != 'first_seq':
    command_parameters.append("-ref1 " + str(ref_seq))

if ref_len != 0:
    command_parameters.append("-ref2 " + str(ref_len))

if BMGE != 0:
    command_parameters.append("-BMGE " + str(BMGE))

if bmge_version != 1:
    command_parameters.append("-bv " + str(bmge_version))

# Module 1: SSR
if p1 != 0.65:
    command_parameters.append("-p1 " + str(p1))

# Module 2: NHHR
if p2 != 0.65:
    command_parameters.append("-p2 " + str(p2))

if m2 != 1:
    command_parameters.append("-m2 " + str(m2))

if s2 != 2:
    command_parameters.append("-s2 " + str(s2))

if a2 != 3:
    command_parameters.append("-a2 " + str(a2))

# Module 3: GIR
if m3 != 1:
    command_parameters.append("-m3 " + str(m3))

if p3 != 0.9:
    command_parameters.append("-p3 " + str(p3))

if aa3 != 30:
    command_parameters.append("-aa3 " + str(aa3))

# Module 4: CGSR
if p4 != 0.5:
    command_parameters.append("-p4 " + str(p4))

if aa4 != 30:
    command_parameters.append("-aa4 " + str(aa4))

# Module 5: PR
if a5 != 3:
    command_parameters.append("-a5 " + str(a5))
if s5 != 2:
    command_parameters.append("-s5 " + str(s5))

print(command_parameters)

parameters = []

if command_parameters:
    if len(command_parameters) > 1:
        parameters = ' '.join(str(x) for x in command_parameters)

print(parameters)


if BMGE > 0:

    if lib.check_bmge_version(bmge_version):
        # LOG #####################################################################
        logger.info('BMGE (version ' + str(bmge_version) + ') will be used. Selected h for BMGE is ' + str(BMGE) +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('BMGE (version ' + str(bmge_version) + ') will be used. Selected h for BMGE is ' + str(BMGE) +
        ' - Total time: ' + str((datetime.now() - startTime)))
        ###########################################################################
    else:
        # LOG #####################################################################
        logger.critical('There is a problem with BMGE. Check the executable is in the ' +
                        'directory and if it is, which version you are using (must match -bv).' + 
                         ' - Total time: ' +  str((datetime.now() - startTime)))
        print('There is a problem with BMGE. Check the executable is in the directory and if it is, which version you are using (must match -bv).')
        sys.exit()
        ###########################################################################


if ref_seq != 'first_seq' and ref_len != 0:
    # LOG #####################################################################
    logger.critical('You cannot apply two rules for the sequence length variable,' +
    'please apply by reference sequence name (ref1) OR by length (ref2), but NOT both.' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print("You put the wrong inputs using ref1 and ref2, you can only use one. please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

len_pre_fasta = []

for pf in pre_fasta_list:
    pf_names, pf_seqs = lib.seqs_extractor(pf)
    len_pre_fasta.append(len(pf_names))

    # LOG #####################################################################
    logger.info(str(pf) + ' loaded' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print(str(pf) + ' loaded' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

if sort_option == 1:
    list_of_fasta = lib.sorting_filelist_HL(len_pre_fasta, pre_fasta_list)

    # LOG #####################################################################
    logger.info('List of fasta sorted from high to low' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    print('List of fasta sorted from high to low' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

if sort_option == 2:
    list_of_fasta = lib.sorting_filelist_LH(len_pre_fasta, pre_fasta_list)

    # LOG #####################################################################
    logger.info('List of fasta sorted from low to high' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    print('List of fasta sorted from low to high' + 
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

BMGE_lines = []

for fasta in list_of_fasta:
    
    real_fasta = fasta.split('.')[0]
    
    # LOG #####################################################################
    logger.info('Seqrutinator initiated for ' + str(fasta) + ' - Total time: ' +
    str((datetime.now() - startTime)))
    print('Seqrutinator initiated for ' + str(fasta) + ' - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################
    
    if not command_parameters:
        os.system("python3 seqrutinator.py -f " + str(fasta))
    
    else:
        if len(command_parameters) == 1:
            print("python3 seqrutinator.py -f " + str(fasta) + " " + str(command_parameters[0]))
            os.system("python3 seqrutinator.py -f " + str(fasta) + " " + str(command_parameters[0]))

        if len(command_parameters) > 1:
            print("python3 seqrutinator.py -f " + str(fasta) + " " + str(parameters))
            os.system("python3 seqrutinator.py -f " + str(fasta) + " " + str(parameters))

    hmms_list = lib.input_files("hmms_*_" + str(real_fasta) + ".*")
    
    for hmms in hmms_list:
        os.system("rm " + str(hmms))
    
    if BMGE > 0:
        bf = open(f'{real_fasta}_BMGE_table.csv')
        bf_lines = bf.readlines()
        
        BMGE_lines.append(str(real_fasta) + "\n")
        
        for bf_l in bf_lines:
            BMGE_lines.append(bf_l)
    
        bf.close()

        os.system("mv " + str(real_fasta) + "_BMGE_table.csv " "Results_" + str(real_fasta))
    
    os.system("cp Results_" + str(real_fasta) + "/" + str(fasta) + " " + str(path))

    # LOG #####################################################################
    logger.info('Seqrutinator finished for ' + str(fasta) + ' - Total time: ' +
    str((datetime.now() - startTime)))
    print('Seqrutinator finished for ' + str(fasta) + ' - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

BMGE_file = open("BMGE_Summary.tsv", "w")
for bl in BMGE_lines:
    BMGE_file.write(bl)
BMGE_file.close()

# LOG #########################################################################
logger.info('SeqYNet finished - Total time: ' +
str((datetime.now() - startTime)))
print('SeqYNet finished - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################

