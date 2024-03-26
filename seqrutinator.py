#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
#                         Script Description                                  #
# name:                                                                       #
#     seqrutinator_v02gamma.py                                                #
#                                                                             #
# description:                                                                #
#     Automatic protein superfamily analysis                                  #
#                                                                             #
# Authors:                                                                    #
#     Nicolas Stocchi, Agustin Amalfitano, Fernando Villareal, Marcelo        #
#     Atencio, Arjen ten Have                                                 #
#                                                                             #
# Version release:                                                            #
#     version 02gamma 2022-07-05                                              #
#     version 02delta 2022-09-14 (name issues corrected)                      #
#     version 02epsilon 2022-12-14 (output file standardized                  #
###############################################################################

# libraries

import os
import sys
import logging
import argparse
import lib_seqrutinator as lib
from datetime import datetime
startTime = datetime.now()

path = os.path.dirname(os.path.abspath(__file__))

###############################################################################
#############################   SEQRUTINATOR   ################################
###############################################################################
# Command Line Arguments
parser = argparse.ArgumentParser()

# General parameters
parser.add_argument('-m', default='12345', help='Pipelines')
parser.add_argument('-f', help='Fasta file (NOTE THAT THIS IS A REQUIRED PARAMETER)')
parser.add_argument('-ali', default=1, help='Use MAFFT G-INS-i (1); Use either MAFFT G-INS-i (n<=500) or Global (n>500) (2); Use FAMSA, recommended only for really large datasets (3)')
parser.add_argument('-ref1', default='first_seq', help='Use to change input reference sequence')
parser.add_argument('-ref2', default=0, help='Use to directly input reference sequence length')
parser.add_argument('-bv', default=1, help='BMGE version 1 (either 1.0 or 1.12) or 2')
parser.add_argument('-BMGE', default=0, help='BMGE deactivated (0), BMGE > 0 is h option for BMGE')

# Module 1: SSR
parser.add_argument('-p1', default=0.65, help='Proportion (0 to 1) of sequence length coverage for SSR')

# Module 2: NHHR
parser.add_argument('-p2', default=0.65, help='Proportion (0 to 1) of sequence length coverage for NHHR')
parser.add_argument('-s2', default=1, help="Mean - alphaSD (1) or Q1 - 1.5IQR (2) for NHHR")
parser.add_argument('-a2', default=3, help='Alpha for NHHR')

# Module 3: GIR
parser.add_argument('-m3', default=1, help='Method one by one (1) or batch (2) for GIR')
parser.add_argument('-p3', default=0.9, help='Proportion (0 to 1) of gaps to define a gap column for GIR (>= VALUE)')
parser.add_argument('-aa3', default=30, help='aa window of contiguos gap columns for GIR')

# Module 4: CGSR
parser.add_argument('-p4', default=0.5, help='Proportion (0 to 1) of gaps to define a gap column for CGSR (>= VALUE)')
parser.add_argument('-aa4', default=30, help='aa window of contiguos gap columns for CGSR')

# Module 5: PR
parser.add_argument('-a5', default=3, help='Alpha for mean - alphaSD (3 is recommended as default option and 2.35 for normal distributions)')
parser.add_argument('-s5', default=1, help='Mean - alphaSD (1) or Q1 - 1.5IQR (2) for PR')

args = vars(parser.parse_args())

# General Parameters
pipeline = str(args['m'])  # Terminal Condition
fasta_file = str(args['f'])  # TC
ali = int(args['ali'])
BMGE = float(args['BMGE'])
bmge_version = int(args['bv'])

# Module 1: SSR
ref_seq = str(args['ref1'])  # TC
ref_len = int(args['ref2'])
p1 = float(args['p1']) * 100

# Module 2: NHHR
p2 = float(args['p2']) * 100
s2 = int(args['s2'])  # t2
a2 = float(args['a2'])  # alpha_1

# Module 3: GIR
m3 = int(args['m3'])
p3 = float(args['p3']) * 100
aa3 = int(args['aa3'])

# Module 4: CGSR
p4 = float(args['p4']) * 100
aa4 = int(args['aa4'])  # t4

# Module 5: PR
a5 = float(args['a5'])
s5 = float(args['s5'])

real_fasta = str((fasta_file.split(".")[:-1])[0])
print(real_fasta)

### LOG #######################################################################
logger = logging.getLogger('Seqrutinator')
hdlr = logging.FileHandler('Seqrutinator_' + str(real_fasta) + '.log', 'w')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.DEBUG)

# LOG #########################################################################
logger.info('Seqrutinator initiated - Total time: ' +
str((datetime.now() - startTime)))
print('Seqrutinator initiated - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################
# Terminal Conditions and parameters logged

# General parameters
# -m parameter
# Split the pipeline in modules

default_pipeline = ['1', '2', '3', '4', '5']

splitted_pipeline = list(pipeline)

rep_pipeline = set(list(splitted_pipeline))

if len(rep_pipeline) < len(splitted_pipeline):
    # LOG #####################################################################
    print("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
    logger.info("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
    sys.exit()
    ###########################################################################

if len(splitted_pipeline) > 5:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -m parameter pipelines,' +
    ' you have more than 5 steps in your pipeline' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print("You put the wrong input on -m parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

if len(splitted_pipeline) == 0:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -m parameter as pipeline,' +
    ' you do not have enough steps in your pipeline - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -m parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('Input pipeline: ' + str(pipeline) +
' - Total time: ' + str((datetime.now() - startTime)))
print('Input pipeline: ' + str(pipeline) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

nsp = 0

for sp in splitted_pipeline:

    nsp = nsp + 1

    # LOG #####################################################################

    if sp == '1':
        logger.info('Step ' + str(nsp) + ': SSR' +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('Step ' + str(nsp) + ': SSR' +
        ' - Total time: ' + str((datetime.now() - startTime)))

    if sp == '2':
        logger.info('Step ' + str(nsp) + ': NHHR' +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('Step ' + str(nsp) + ': NHHR' +
        ' - Total time: ' + str((datetime.now() - startTime)))

    if sp == '3':
        logger.info('Step ' + str(nsp) + ': GIR' +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('Step ' + str(nsp) + ': GIR' +
        ' - Total time: ' + str((datetime.now() - startTime)))

    if sp == '4':
        logger.info('Step ' + str(nsp) + ': CGSR' +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('Step ' + str(nsp) + ': CGSR' +
        ' - Total time: ' + str((datetime.now() - startTime)))

    if sp == '5':
        logger.info('Step ' + str(nsp) + ': PR' +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print('Step ' + str(nsp) + ': PR' +
        ' - Total time: ' + str((datetime.now() - startTime)))

    ###########################################################################

    if sp not in default_pipeline:

        # LOG #################################################################
        logger.critical('You put the wrong input in -m parameter,' + str(sp) +
        ' is not an existing module. - Total time: ' +
        str((datetime.now() - startTime)))
        print("You put the wrong input on -m parameter, please check the log and try again" +
        ' - Total time: ' + str((datetime.now() - startTime)))
        sys.exit()
        #######################################################################

# -f parameter
if fasta_file == 'fasta_file':

    # LOG #####################################################################
    logger.critical('You put the wrong input in -f parameter as fasta file,' +
    ' you did not give an input fasta file - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -f parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

else:

    # LOG #####################################################################
    logger.info('Fasta file: ' + str(fasta_file) +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print('Fasta file: ' + str(fasta_file) +
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

    # Extract names and sequences from MSA
    names, seqs =lib.seqs_extractor_desc(fasta_file)
 
    IUPAC = ["-", "a", "A", "c", "C", "d", "D", "e", "E", "f", "F", "g", "G",
    "h", "H", "i", "I", "k", "K", "l", "L", "m", "M", "n", "N", "p", "P", "q",
    "Q", "r", "R", "s", "S", "t", "T", "v", "V", "w", "W", "y", "Y"]
    for id, seq in zip(names,seqs):
        for char in seq:
            if char not in IUPAC:
                    # LOG #####################################################
                    logger.critical('Fasta file has sequences with non IUPAC symbols, for example ' + str(char) + ' in sequence ' + str(id) + 
                    '\nPlease check your data and remove all the sequences with non IUPAC characters')
                    print('Fasta file has sequences with non IUPAC symbols, for example ' + str(char) + ' in sequence ' + str(id) + 
                    '\nPlease check your data and remove all the sequences with non IUPAC characters')
                    sys.exit()
                    ###########################################################
        if ' ' in str(id):
            # LOG #####################################################
            logger.critical('One or more IDs in the sequences of your fasta file have spaces. This is problematic for Seqrutinator. Please remove spaces in the IDs or consider to rename sequences altogether (recommended)')
            print('One or more IDs in the sequences of your fasta file have spaces. This is problematic for Seqrutinator')
            print('Please remove spaces in the IDs or consider to rename sequences altogether (recommended)')        
            sys.exit()
            ###########################################################

    names, seqs =lib.seqs_extractor(fasta_file)

    wog_con = 0

    for seq in seqs:

        if "-" in seq:
            wog_con = wog_con + 1
            break

    if wog_con == 0:  # If fasta input is not aligned
        if ali == 1:
            lib.MAFFT_ginsi(fasta_file)
        if ali == 2:
            lib.MAFFT_seqrutinator(fasta_file)
        if ali == 3:
            lib.FAMSA(fasta_file)

        mafft_name = lib.rename_mafft(fasta_file)
        fasta_file = mafft_name

# ali parameter

if ali < 1 or ali > 3:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -ali parameter,' +
    ' you did not specify MAFFT options - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -ali parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

if ali == 1:

    # LOG #####################################################################
    logger.info('MAFFT GINSI will be used in all modules - Total time: ' +
    str((datetime.now() - startTime)))
    print('MAFFT GINSI will be used in all modules - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

if ali == 2:

    # LOG #####################################################################
    logger.info('MAFFT GINSI (n<= 500) or Global (>500) will be used in all modules - Total time: ' +
    str((datetime.now() - startTime)))
    print('MAFFT GINSI (n<= 500) or Global (>500) will be used in all modules - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

if ali == 3:

    # LOG #####################################################################
    logger.info('FAMSA2 will be used in all modules - Total time: ' +
    str((datetime.now() - startTime)))
    print('FAMSA2 will be used in all modules - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

# BMGE parameter

if BMGE < 0 or BMGE > 1:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -BMGE parameter,' +
    ' you did not specify if you want to do BMGE - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -BMGE parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

if BMGE == 0:

    # LOG #####################################################################
    logger.info('BMGE will not be done - Total time: ' +
    str((datetime.now() - startTime)))
    print('BMGE will not be done - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

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

# -ref1 or ref2 parameters
if ref_seq != 'first_seq' and ref_len != 0:
    # LOG #####################################################################
    logger.critical('You cannot apply two rules for the sequence length variable,' +
    'please apply by reference sequence name (ref1) OR by length (ref2), but NOT both.' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print("You put the wrong inputs using ref1 and ref2, you can only use one. please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

list_rs = []

if ref_seq != 'first_seq':

    wog_fasta = lib.remove_gaps(fasta_file)
    ref_names, ref_seqs = lib.seqs_extractor(wog_fasta)

    for rf1, rf2 in zip(ref_names, ref_seqs):
        if rf1 == ref_seq:
            ref_size = len(list(rf2))
            list_rs.append(ref_size)
            break

    # LOG #####################################################################
    logger.info('You chose ' + str(ref_seq) + ' as your reference sequence' +
    ', so Seqrutinator will use ' + str(ref_seq) + ' rather than the first sequence of your MSA. - Total time: ' + str((datetime.now() - startTime)))
    print('You chose ' + str(ref_seq) + ' as your reference sequence' +
    ', so Seqrutinator will use ' + str(ref_seq) + ' rather than the first sequence of your MSA. - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################
    
elif ref_len != 0:
    # LOG #####################################################################
    logger.info('You chose that the SSR cut-off is determined by the lenght you fixed at ' + str(ref_len) + ' by which the concept of reference sequence does no apply for SSR. - Total time: ' + str((datetime.now() - startTime)))
    print('You chose that the SSR cut-off is determined by the lenght you fixed at ' + str(ref_len) + ' by which the concept of reference sequence does no apply for SSR. - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################
    
else:
    wog_fasta = lib.remove_gaps(fasta_file)
    ref_names, ref_seqs = lib.seqs_extractor(wog_fasta)

    # LOG #####################################################################
    logger.info('You used -ref1 by default, so Seqrutinator will use the first sequence of your MSA (' + str(ref_names[0]) +
    ') as reference. If you want you use another sequence as reference, indicate the name of the reference sequence using the -ref1 parameter.')
    logger.info('Total time: ' + str((datetime.now() - startTime)))
    print('You chose the first sequence of your MSA by default: ' + str(ref_names[0]) +
    ' - Total time: ' + str((datetime.now() - startTime)))
    ###########################################################################

    ref_size = len(list(ref_seqs[0]))
    list_rs.append(ref_size)

if ref_len != 0:
    list_rs.append(ref_len)

    # LOG #####################################################################
    logger.info('You chose ' + str(ref_len) +
    ' as the reference sequence length - Total time: ' +
    str((datetime.now() - startTime)))
    print('You chose ' + str(ref_len) +
    ' as the reference sequence length - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

ref_size = int(list_rs[0])

# LOG #########################################################################
logger.info('Reference sequence length: ' + str(ref_size) +
' - Total time: ' + str((datetime.now() - startTime)))
print('Reference sequence length: ' + str(ref_size) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

###############################################################################
# Module 1 parameters: SSR

# -p1 parameter
if p1 < 0 or p1 > 100:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -p1 parameter,' +
    ' you chose a wrong percentage for SSR' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print("You put the wrong input on -p1 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('Coverage of reference sequence length for SSR: ' + str(p1) +
'% - Total time: ' + str((datetime.now() - startTime)))
print('Coverage of reference sequence length for SSR: ' + str(p1) +
'% - Total time: ' + str((datetime.now() - startTime)))
###############################################################################
# Module 2 parameters: NHHR

# -p2 parameter
if p2 < 0 or p2 > 100:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -p2 parameter,' +
    ' you chose a wrong percentage for NHHR' +
    ' - Total time: ' + str((datetime.now() - startTime)))
    print("You put the wrong input on -p2 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('Coverage of reference sequence length for NHHR: ' + str(p2) +
'% - Total time: ' + str((datetime.now() - startTime)))
print('Coverage of reference sequence length for NHHR: ' + str(p2) +
'% - Total time: ' + str((datetime.now() - startTime)))
###############################################################################
# -s2 parameter

if s2 < 1 or s2 > 2:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -s2 parameter,' +
    ' you did not specify the threshold: Mean - 3SD (1) or Q1 - 1.5IQR (2) for NHHR - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -s2 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

if s2 == 1:

    # LOG #####################################################################
    logger.info('Mean - ' + str(a2) + 'SD was chosen as a threshold for NHHR - Total time: ' +
    str((datetime.now() - startTime)))
    print('Mean - ' + str(a2) + 'SD was chosen as a threshold for NHHR - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

if s2 == 2:

    # LOG #####################################################################
    logger.info('Q1 - 1.5IQR was chosen as a threshold for NHHR - Total time: ' +
    str((datetime.now() - startTime)))
    print('Q1 - 1.5IQR was chosen as a threshold for NHHR - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

# -a2 parameter

# LOG #########################################################################
logger.info('Selected alpha for NHHR: ' + str(a2) +
' - Total time: ' + str((datetime.now() - startTime)))
print('Selected alpha for NHHR: ' + str(a2) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

###############################################################################
# Module 3 parameters: GIR
# -m3 parameter

if m3 < 1 or m3 > 2:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -m3 parameter,' +
    ' you did not specify the method: (1) one by one or (2) batch for GIR - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -m3 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

if m3 == 1:

    # LOG #####################################################################
    logger.info('One by one was chosen for GIR - Total time: ' +
    str((datetime.now() - startTime)))
    print('One by one was chosen for GIR - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

if m3 == 2:

    # LOG #####################################################################
    logger.info('Batch was chosen for GIR - Total time: ' +
    str((datetime.now() - startTime)))
    print('Batch was chosen for GIR - Total time: ' +
    str((datetime.now() - startTime)))
    ###########################################################################

# -p3 parameter

if p3 < 0 or p3 > 100:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -p3 parameter,' +
    ' you chose a wrong percentage for considering gap columns for GIR - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -p3 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('Percentage of gap column for GIR: ' + str(p3) +
'% - Total time: ' + str((datetime.now() - startTime)))
print('Percentage of gap column for GIR: ' + str(p3) +
'% - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

# -aa3 parameter

if aa3 < 0:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -aa3 parameter,' +
    ' you chose a wrong aa window length for GIR- Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -aa3 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('aa window of contiguos gap columns for GIR: ' + str(aa3) +
' - Total time: ' + str((datetime.now() - startTime)))
print('aa window of contiguos gap columns for GIR: ' + str(aa3) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

# -p4 parameter

if p4 < 0 or p4 > 100:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -p4 parameter,' +
    ' you chose a wrong percentage for CGSR - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -p4 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('Minimum admisible percentage of gaps for a gap column for CGSR: ' +
str(p4) +  '% - Total time: ' + str((datetime.now() - startTime)))
print('Minimum admisible percentage of gaps for a gap column for CGSR: ' +
str(p4) +  '% - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

# -aa4 parameter

if aa4 < 0:

    # LOG #####################################################################
    logger.critical('You put the wrong input in -aa4 parameter,' +
    ' you chose a wrong min admisible threshold for CGSR - Total time: ' +
    str((datetime.now() - startTime)))
    print("You put the wrong input on -aa4 parameter, please check the log and try again" +
    ' - Total time: ' + str((datetime.now() - startTime)))
    sys.exit()
    ###########################################################################

# LOG #########################################################################
logger.info('aa window of contiguos gap columns for CGSR: ' + str(aa4) +
' - Total time: ' + str((datetime.now() - startTime)))
print('aa window of contiguos gap columns for CGSR: ' + str(aa4) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

# PR parameters
# -s5 parameter

# LOG #########################################################################
logger.info('Selected alpha for PR: ' + str(a5) +
' - Total time: ' + str((datetime.now() - startTime)))
print('Selected alpha for PR: ' + str(a5) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################

# -a5 parameter
# LOG #########################################################################
logger.info('Selected alpha for PR: ' + str(a5) +
' - Total time: ' + str((datetime.now() - startTime)))
print('Selected alpha for PR: ' + str(a5) +
' - Total time: ' + str((datetime.now() - startTime)))
###############################################################################



original_fasta = fasta_file
real_fasta = str((fasta_file.split(".")[:-1])[0])
#print(real_fasta)

temporary_fasta = "temporary_fasta"

os.system("mkdir Results_" + str(real_fasta))

for sp in splitted_pipeline:
    if sp == '1':
        os.system("mkdir Results_" + str(real_fasta) + "/1_SSR")
    if sp == '2':
        os.system("mkdir Results_" + str(real_fasta) + "/2_NHHR")
    if sp == '3':
        os.system("mkdir Results_" + str(real_fasta) + "/3_GIR")
    if sp == '4':
        os.system("mkdir Results_" + str(real_fasta) + "/4_CGSR")
    if sp == '5':
        os.system("mkdir Results_" + str(real_fasta) + "/5_PR")
dst = path + "/Results_" + str(real_fasta)
sst = path + "/Results_" + str(real_fasta) + "/1_SSR"
nst = path + "/Results_" + str(real_fasta) + "/2_NHHR"
gst = path + "/Results_" + str(real_fasta) + "/3_GIR"
cst = path + "/Results_" + str(real_fasta) + "/4_CGSR"
pst = path + "/Results_" + str(real_fasta) + "/5_PR"

if BMGE > 0:
    BMGE_lines_split = []

used_modules = []
num_seqs = []

for sp in splitted_pipeline:

    if sp != splitted_pipeline[0]:  # First module
        fasta_file = temporary_fasta

    if sp == '1':  # SSR

        spm = 'SSR'

        if spm not in used_modules:
            used_modules.append(spm)

        else:
            # LOG #############################################################
            print("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            logger.info("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            sys.exit()
            ###################################################################

        # LOG #################################################################
        partialTime = datetime.now()
        logger.info('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        #######################################################################

        SSR_output, SSR_removed = lib.SSR(fasta_file, real_fasta, ref_seq, ref_size, p1, ali)

        # LOG #################################################################
        logger.info("Sequence removed: " + str(ref_seq) +
        ' - Total time: ' + str((datetime.now() - startTime)))
        print("Sequence removed: " + str(ref_seq) +
        ' - Total time: ' + str((datetime.now() - startTime)))
        #######################################################################

        if SSR_removed == []:

            temporary_fasta = fasta_file
            num_seqs.append(0)

            # LOG #############################################################
            logger.info('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            print('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            ###################################################################

        else:

            temp_SSR = lib.rename_mafft(SSR_output)
            temporary_fasta = temp_SSR

            if BMGE > 0:

                # LOG #########################################################
                logger.info('BMGE initiated for ' + str(temp_SSR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('BMGE initiated for ' + str(temp_SSR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

                try:
                    lib.BMGE(bmge_version,temp_SSR, BMGE, "BMGE_output.txt")
                    BMGE_alert, seq_bef, seq_aft, BMGE_bef,BMGE_aft = lib.BMGE_reader(bmge_version,"BMGE_output.txt")

                    BMGE_lines_split.append(str(spm) + "\t" + str(seq_bef) + "\t" + str(seq_aft) + "\t" + str(BMGE_bef) + "\t" + str(BMGE_aft) + "\t" + str(BMGE_alert) + "\n")

                    # LOG #####################################################
                    logger.info('BMGE finished for ' + str(temp_SSR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE finished for ' + str(temp_SSR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

                except:
                    BMGE_lines_split.append(str(spm) + "\tSomething went wrong with BMGE.\n")

                    # LOG #####################################################
                    logger.info('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################


            r_lines = []

            r1 = open(str(spm) + '_' + str(real_fasta) + '_removed.txt', 'r')
            for line in r1:
                r_lines.append(line)
            r1.close()

            num_seqs.append(int(len(r_lines)))

            if len(r_lines) == 1:

                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequence was removed with ' + str(spm) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequence was removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

            else:

                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequences were removed with ' + str(spm) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

        # LOG #################################################################
        logger.info('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        #######################################################################

    if sp == '2':
        spm = 'NHHR'

        if spm not in used_modules:
            used_modules.append(spm)

        else:
            # LOG #############################################################
            print("WARNING: YOU CANNOT REPEAT MODULES, FILENAMES ISSUES")
            logger.info("WARNING: YOU CANNOT REPEAT MODULES, FILENAMES ISSUES")
            sys.exit()
            ###################################################################

        # LOG #################################################################
        partialTime = datetime.now()
        logger.info('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        #######################################################################

        NHHR_output, NHHR_removed = lib.NHHR(fasta_file, real_fasta, a2, s2, ref_size, p2, ali)

        if NHHR_removed == []:

            temporary_fasta = fasta_file
            num_seqs.append(0)

            # LOG #############################################################
            logger.info('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            print('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            ###################################################################

        else:
            temp_NHHR = lib.rename_mafft(NHHR_output)
            temporary_fasta = temp_NHHR

            if BMGE > 0:

                # LOG #########################################################
                logger.info('BMGE initiated for ' + str(temp_NHHR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('BMGE initiated for ' + str(temp_NHHR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

                try:
                    lib.BMGE(bmge_version,temp_NHHR, BMGE, "BMGE_output.txt")
                    BMGE_alert, seq_bef, seq_aft, BMGE_bef, BMGE_aft = lib.BMGE_reader(bmge_version,"BMGE_output.txt")

                    BMGE_lines_split.append(str(spm) + "\t" + str(seq_bef) + "\t" + str(seq_aft) + "\t" + str(BMGE_bef) + "\t" + str(BMGE_aft) + "\t" + str(BMGE_alert) + "\n")

                    # LOG #####################################################
                    logger.info('BMGE finished for ' + str(temp_NHHR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE finished for ' + str(temp_NHHR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

                except:
                    BMGE_lines_split.append(str(spm) + "\tSomething went wrong with BMGE.\n")

                    # LOG #####################################################
                    logger.info('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

            r_lines = []
            r1 = open(str(spm) + '_' + str(real_fasta) + '_removed.txt', 'r')
            for line in r1:
                r_lines.append(line)
            r1.close()

            num_seqs.append(int(len(r_lines)))

            # LOG #############################################################
            if len(r_lines) == 1:
                logger.info(str(len(r_lines)) + ' sequence was removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequence was removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
            else:
                logger.info(str(len(r_lines)) + ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
            ###################################################################

        # LOG #################################################################
        logger.info('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        #######################################################################

    if sp == '3':  # GIR

        spm = 'GIR'

        if spm not in used_modules:
            used_modules.append(spm)

        else:
            # LOG #############################################################
            print("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            logger.info("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            sys.exit()
            ###################################################################

        # LOG #################################################################
        partialTime = datetime.now()
        logger.info('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        #######################################################################

        GIR_output, GIR_removed, GIR_prob = lib.GIR(m3, fasta_file, real_fasta, p3, aa3, ali)

        if GIR_prob:

            if len(GIR_prob) > 1:
                # plural
                # LOG #########################################################
                logger.info('WARNING: ' + str(len(GIR_prob)) +
                ' gap blocks were detected with a global score > ' + str(aa3) +
                ', please check the last GIR_Problematic_Sequences file for more info.' +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('WARNING: ' + str(len(GIR_prob)) +
                ' gap blocks were detected with a global score > ' + str(aa3) +
                ', please check the last GIR_Problematic_Sequences file for more info.' +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

            if len(GIR_prob) == 1:
                # singular
                # LOG #########################################################
                logger.info('WARNING: ' + str(len(GIR_prob)) +
                ' gap block was detected with a global score > ' + str(aa3) +
                ', please check the last GIR_Problematic_Sequences file for more info.' +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('WARNING: ' + str(len(GIR_prob)) +
                ' gap block was detected with a global score > ' + str(aa3) +
                ', please check the last GIR_Problematic_Sequences file for more info.' +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

        if GIR_removed == []:

            temporary_fasta = fasta_file
            num_seqs.append(0)

            # LOG #############################################################
            logger.info('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            print('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            ###################################################################

        else:
            temp_GIR = lib.rename_mafft(GIR_output)
            temporary_fasta = temp_GIR

            if BMGE > 0:

                # LOG #########################################################
                logger.info('BMGE initiated for ' + str(temp_GIR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('BMGE initiated for ' + str(temp_GIR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

                try:
                    lib.BMGE(bmge_version,temp_GIR, BMGE, "BMGE_output.txt")
                    BMGE_alert, seq_bef, seq_aft, BMGE_bef, BMGE_aft = lib.BMGE_reader(bmge_version,"BMGE_output.txt")

                    BMGE_lines_split.append(str(spm) + "\t" + str(seq_bef) + "\t" + str(seq_aft) + "\t" + str(BMGE_bef) + "\t" + str(BMGE_aft) + "\t" + str(BMGE_alert) + "\n")

                    # LOG #####################################################
                    logger.info('BMGE finished for ' + str(temp_GIR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE finished for ' + str(temp_GIR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

                except:
                    BMGE_lines_split.append(str(spm) + "\tSomething went wrong with BMGE.\n")

                    # LOG #####################################################
                    logger.info('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

            r_lines = []
            r1 = open(str(spm) + '_' + str(real_fasta) + '_removed.txt', 'r')
            
            for line in r1:
                r_lines.append(line)
            r1.close()

            num_seqs.append(int(len(r_lines)))

            if len(r_lines) == 1:
                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequence was removed with ' + str(spm) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequence was removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

            else:
                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

        # LOG #################################################################
        logger.info('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        #######################################################################

    if sp == '4':  # CGSR

        spm = 'CGSR'

        if spm not in used_modules:
            used_modules.append(spm)

        else:
            # LOG #############################################################
            print("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            logger.info("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            sys.exit()
            ###################################################################

        # LOG #################################################################
        partialTime = datetime.now()
        logger.info('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        #######################################################################

        CGSR_data, CGSR_output, CGSR_removed = lib.CGSR(fasta_file, real_fasta, p4, aa4, ali)

        for dat, l_dat in zip(CGSR_data, range(1, len(CGSR_data) + 1)):
            f = open(str(spm) + "_It" + str(l_dat) + "_" + str(real_fasta) + "_data.txt", "w")
            for d in dat:
                f.write(str(d))
            f.close()

        if CGSR_removed == []:

            temporary_fasta = fasta_file
            num_seqs.append(0)

            # LOG #############################################################
            logger.info('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            print('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            ###################################################################

        else:
            lib.fetching(CGSR_removed, original_fasta, "CGSR_" + str(real_fasta) + "_removed.fsa")
            fsa_CGSR = lib.remove_gaps("CGSR_" + str(real_fasta) + "_removed.fsa")

            temp_CGSR = lib.rename_mafft(CGSR_output)
            temporary_fasta = temp_CGSR

            if BMGE > 0:

                # LOG #########################################################
                logger.info('BMGE initiated for ' + str(temp_CGSR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('BMGE initiated for ' + str(temp_CGSR) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

                try:
                    lib.BMGE(bmge_version,temp_CGSR, BMGE, "BMGE_output.txt")
                    BMGE_alert, seq_bef, seq_aft, BMGE_bef, BMGE_aft = lib.BMGE_reader(bmge_version,"BMGE_output.txt")

                    BMGE_lines_split.append(str(spm) + "\t" + str(seq_bef) + "\t" + str(seq_aft) + "\t" + str(BMGE_bef) + "\t" + str(BMGE_aft) + "\t" + str(BMGE_alert) + "\n")

                    # LOG #####################################################
                    logger.info('BMGE finished for ' + str(temp_CGSR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE finished for ' + str(temp_CGSR) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

                except:
                    BMGE_lines_split.append(str(spm) + "\tSomething went wrong with BMGE.\n")

                    # LOG #####################################################
                    logger.info('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

            r_lines = []

            r1 = open(str(spm) + '_' + str(real_fasta) + '_removed.txt', 'r')
            for line in r1:
                r_lines.append(line)
            r1.close()

            num_seqs.append(int(len(r_lines)))

            if len(r_lines) == 1:
                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequence was removed with ' + str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequence was removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

            else:
                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequences were removed with ' + str(spm) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

        # LOG #################################################################
        logger.info('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        #######################################################################

    if sp == '5':  # PR

        spm = 'PR'

        if spm not in used_modules:
            used_modules.append(spm)

        else:
            # LOG #############################################################
            print("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            logger.info("WARNING: YOU CAN NOT REPEAT MODULES, FILENAMES ISSUES")
            sys.exit()
            ###################################################################

        # LOG #################################################################
        partialTime = datetime.now()
        logger.info('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) +
        ') initiated - Total time: ' + str((datetime.now() - startTime)))
        #######################################################################

        PR_removed, PR_removed_file, PR_remain, PR_fsa, PR_faa = lib.PR(fasta_file, s5, a5, ali)

        if PR_removed == []:

            temporary_fasta = fasta_file
            num_seqs.append(0)

            # LOG #############################################################
            logger.info('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            print('No sequences were removed with ' + str(spm) +
            ' - Total time: ' + str((datetime.now() - startTime)))
            ###################################################################

        else:
            temporary_fasta = PR_faa

            if BMGE > 0:

                # LOG #########################################################
                logger.info('BMGE initiated for ' + str(PR_faa) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print('BMGE initiated for ' + str(PR_faa) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                ###############################################################

                try:
                    lib.BMGE(bmge_version,PR_faa, BMGE, "BMGE_output.txt")
                    BMGE_alert, seq_bef, seq_aft, BMGE_bef, BMGE_aft = lib.BMGE_reader(bmge_version,"BMGE_output.txt")

                    BMGE_lines_split.append(str(spm) + "\t" + str(seq_bef) + "\t" + str(seq_aft) + "\t" + str(BMGE_bef) + "\t" + str(BMGE_aft) + "\t" + str(BMGE_alert) + "\n")

                    # LOG #####################################################
                    logger.info('BMGE finished for ' + str(PR_faa) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE finished for ' + str(PR_faa) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

                except:
                    BMGE_lines_split.append(str(spm) + "\tSomething went wrong with BMGE.\n")

                    # LOG #####################################################
                    logger.info('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    print('BMGE did not finished for ' + str(spm) +
                    ' - Total time: ' + str((datetime.now() - startTime)))
                    ###########################################################

            r1 = open(str(spm) + '_' + str(real_fasta) + '_removed.txt', 'w')
            for tr in PR_removed:
                r1.write(str(tr) + "\n")
            r1.close()

            r_lines = PR_removed
            num_seqs.append(int(len(r_lines)))

            if len(r_lines) == 1:
                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequence was removed with ' + str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequence was removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

            else:
                # LOG #########################################################
                logger.info(str(len(r_lines)) +
                ' sequences were removed with ' + str(spm) +
                ' - Total time: ' + str((datetime.now() - startTime)))
                print(str(len(r_lines)) + ' sequences were removed with ' +
                str(spm) + ' - Total time: ' +
                str((datetime.now() - startTime)))
                ###############################################################

        # LOG #################################################################
        logger.info('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        print('Module ' + str(sp) + ' (' + str(spm) + ') finished' +
        " - Partial time: " + str((datetime.now() - partialTime)) +
        " - Total time: " + str((datetime.now() - startTime)))
        #######################################################################

if BMGE > 0:
    bfo = open(str(real_fasta) + "_BMGE_table.csv", "w")
    bfo.write("Modules\tSeqs Before\tSeqs After\tCols Before\tCols After\tALERTS\n")

    for bls in BMGE_lines_split:
        bfo.write(bls)

    bfo.close()

    os.system("rm BMGE_output.txt")

# Move files to Results directory

if temporary_fasta != "temporary_fasta":
    final_res = "RESULT_" + str(real_fasta) + ".faa"

    os.system("cp " + str(temporary_fasta) + " " + str(final_res))

    res_names, res_seqs = lib.seqs_extractor(final_res)

    os.system("mv " + str(final_res) + " " + str(dst))

pr_names, pr_seqs = lib.seqs_extractor(original_fasta)

# LOG #########################################################################
logger.info('Total Sequences\t' + str(len(pr_names)))
logger.info('Modules\tRemoved Sequences')
###############################################################################

for usm, nse in zip(used_modules, num_seqs):

    # LOG #####################################################################
    logger.info(str(usm) + '\t' + str(nse))
    ###########################################################################

    usm_list = lib.input_files(str(usm) + "_*")
    maf_list = lib.input_files("GIR*")

    for usi in usm_list:
        if real_fasta in usi:
            os.system("mv " + str(usi) + " " + str(dst))

    for mfi in maf_list:
        if real_fasta in mfi:
            os.system("mv " + str(mfi) + " " + str(dst))

if temporary_fasta != "temporary_fasta":
    # LOG #####################################################################
    logger.info('RESULT\t' + str(len(res_names)))
    ###########################################################################

# LOG #########################################################################
logger.info('Seqrutinator finished - Total time: ' +
str((datetime.now() - startTime)))
print('Seqrutinator finished - Total time: ' +
str((datetime.now() - startTime)))
###############################################################################

move_files = lib.input_files(str(real_fasta) + ".*")
for mf in move_files:
    if mf != original_fasta:
        os.system("mv " + str(mf) + " " + str(dst))
    else:
        os.system("cp " + str(mf) + " " + str(dst))

os.system("mv Seqrutinator_" + str(real_fasta) + ".log " + str(dst))



for sp in splitted_pipeline:
    if sp == '1':
        os.system("mv Results_" + str(real_fasta) + "/SSR* " + str(sst))
    if sp == '2':
        os.system("mv Results_" + str(real_fasta) + "/NHHR* " + str(nst))
    if sp == '3':
        os.system("mv Results_" + str(real_fasta) + "/GIR* " + str(gst))
    if sp == '4':
        os.system("mv Results_" + str(real_fasta) + "/CGSR* " + str(cst))
    if sp == '5':
        os.system("mv Results_" + str(real_fasta) + "/PR* " + str(pst))
        os.system("rm hmm_*")
        os.system("rm hmms_*")
if BMGE > 0:
    os.system("cp " + str(real_fasta) + "_BMGE_table.csv " + str(dst))