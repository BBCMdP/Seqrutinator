# Seqrutinator
Public Repo for the Seqrutinator tool and its python dependencias

Seqrutinator is an objective, flexible pipeline for the scrutiny of sequence sets from complex, eukaryotic protein superfamilies. It removes sequences from pseudogenes, incorrect gene models or with sequencing errors. Testing Seqrutinator on major superfamilies BAHD, CYP and UGT correctly removed 1.94% of SwissProt entries, 14% of entries from the model plant Arabidopsis thaliana but 80% of entries from Pinus taeda’s recent complete proteome. Most of the removed sequences were partial. The scrutiny of BAHDomes, CYPomes and UGTomes obtained from 16 plant proteomes show consistent numbers suggesting good performance. MSAs, phylogenies and particularly functional clustering improved drastically upon Seqrutinator application.

Seqrutinator is (about to be) under review but can be accessed as preprint at https://doi.org/10.1101/2022.03.22.485366
Below you will find the contents of Supplemental document 1 which is basically a detailed descripion of Seqrutinator and accompaying scripts. My apologies for the fact that the output tables are no really tables here

           Introduction
Seqrutinator is a modular pipeline written in Python 3 to identify and remove protein superfamily sequences that likely either correspond to a non-functional homologue (NFH) or are incorrect, the latter due to either sequencing or gene modelling errors. Seqrutinator is made for complex single-domain protein superfamilies. It requires intermediate or large datasets (~>50 input sequences). This document describes in detail how the modules function, how the pipeline is made and which outputs are generated. Seqrutinator requires a single input file with the sequence set in fasta format, either aligned or not. Dependencies are described/listed in Materials and Methods, which includes shell commands that are used to run software packages such as MAFFT.
	Each module is accessed by a number, according to the default pipeline 12345: 1) Short Sequence Remover (SSR) → 2) Non-Homologous Hit Remover (NHHR) → 3) Gap Instigator Remover (GIR) → 4) Continuous Gap Sequence Remover → 5) Pseudogene Remover (PR). Modules 2 to 5 need an MSA as input but when applied as first module, the script automatically aligns the sequences if they are not aligned. We recommend the user trims the MSA such that it represents the mature protein and lacks non-homologous terminal sequences. Seqrutinator is therefore accompanied by MUFASA (MUltiple FASta Aligner) which is used to obtain the initial crude but aligned sequence set. For trimming and in order to include a general reference we recommend aligning a single sequence for which a structure has been resolved to the MSA by means of MAFFT --add. Seqrutinator can be applied in batch (i.e. using multiple input files) using an additional script named SeqYNet. All three scripts are described below. Descriptions concern the default settings but we describe many options to overrule default settings. 

           MUFASA (MUltiple FASta Aligner)
MUFASA performs hmmsearch and aligns positively identified sequences. It only comes with options regarding input files as demonstrated by:
 
$ python3 MUFASA.py -h
usage: MUFASA.py [-h] [-i I] [-ext EXT] [-c C]

options:
  -h, --help  show this help message and exit
  -i I        HMMER profile, default = hmm_profile
  -ext EXT    Target files extension. default = *.fsa
  -c C        Cores, default = 4

MUFASA takes as input a single HMMER profile and performs separate hmmsearch of all fasta files available in the same folder as where MUFASA is executed. By default, all files with extension fsa will be used as input sequence set and hmm_profile is the default HMMER profile. For each input sequence set, it fetches all identified sequences and saves them in a single fasta file, which represent a superfamily proteome. Next, all obtained superfamily proteomes are separately aligned using MAFFT G-INS-i. 

Table 1: Input and Output of MUFASA
Input
Output
bug1.fsa
bug1_hits.txt
bug1_hits.fsa; .faa

Hits_summary.tsv

MUFASA.log
bug2.fsa
bug2_hits.txt
bug2_hits.fsa; .faa



           Seqrutinator
Seqrutinator is a fully flexible, modular pipeline. It consists of a main python script with the actual modules and some recurring algorithms in a python module (lib_seqrutinator.py). Seqrutinator comes with a large number of options. Parameters have numbers according to the number of the module, Parameters that are not numbered are general. Parameters are as demonstrated by:

python3 seqrutinator.py -h
usage: seqrutinator.py [-h] [-m M] [-f F] [-ali ALI] [-ref1 REF1]
                         [-ref2 REF2] [-BMGE BMGE] [-p1 P1] [-p2 P2] [-s2 S2]
                         [-a2 A2] [-m3 M3] [-p3 P3] [-aa3 AA3] [-p4 P4]
                         [-aa4 AA4] [-a5 A5]

optional arguments:
  -h, --help  show this help message and exit
  -m M        Pipelines
  -f F        Fasta file NOTE THAT THIS IS A REQUIRED PARAMETER
  -ali ALI    Use MAFFT G-INS-i (1); Use either MAFFT G-INS-i (n<=500) or
              Global (n>500) (2); Use FAMSA, recommended only for really large
              datasets (3)
  -ref1 REF1  Use to change input reference sequence
  -ref2 REF2  Use to directly input reference sequence length
  -BMGE BMGE  BMGE deactivated (0), BMGE > 0 is h option for BMGE
  -p1 P1      Proportion (0 to 1) of sequence length coverage for SSR
  -p2 P2      Proportion (0 to 1) of sequence length coverage for NHHR
  -s2 S2      Mean - alphaSD (1) or Q1 - 1.5IQR (2)
  -a2 A2      Alpha for NHHR
  -m3 M3      Method one by one (1) or batch (2) for GIR
  -p3 P3      Proportion (0 to 1) of gaps to define a gap column for GIR (>=
              VALUE)
  -aa3 AA3    aa window of contiguos gap columns for GIR

1 The Short Sequence Remover Module
SSR opens the offered fasta file and determines the length l of the first sequence,  the default reference sequence. This can be overruled by ref1 or ref2. Next, it determines the length of each sequence and separately saves good and bad sequences using the default threshold of 0.65*l, which can be overruled by p1. Good sequences are aligned by MAFFT G-INS-i. The additional output consists of a file that contains the accession codes of the removed sequences and a file that shows the inclusion threshold and the lengths of the removed sequences. 
2 The Non-Homologous Hit Remover Module
NHHR removes outliers identified by low HMMER score. It constructs a HMMER profile of the input MSA and uses it to run a hmmsearch against the unaligned sequences. It directly removes all sequences that score 0 and then uses the total scores of the sequences to perform a distribution analysis. In order to prevent removing short sequences (in case SSR does not precede NHHR in the applied pipe), sequences shorter than 65% (overrule by p2) of the reference sequence are excluded from this analysis. It determines the mean and the standard deviation of all the data in the array.
	Finally, it determines which sequences score below and which above the threshold, set at 3σ, and saves removed and accepted sequences separately, accepted sequences are aligned. NHHR is not iterated since its objective is to remove sequences that are not homologous and will be removed in a single screening, under the assumption that they score low. The additional output consists a file that contains the accession codes of the removed sequences, the hmmsearch ouputs and a plot that shows a distribution of the calculated scores with the applied cut-off threshold.
	Threshold settings can be overruled. Any value alpha can be used for the sigma rule (a2). Alternatively, the cut-off can be based on interquartiles (s2). The Inter Quartile Range or IQR is defined as the difference between the middle of the first half of the sequence scores and the middle of the second half of the sequence scores. This provides another measure that is often used to define outliers, such as the lower 1.5 * Inter Quartile Range (IQR) whisker (< Q1 – 1.5 x IQR). 
3 The Gap Instigator Remover
GIR removes sequences that instigate large (>29, overrule using aa3) gaps as detected in the offered MSA. This comes with two uncertainties of which only the firts one is tackled automatically. A “sequence-specific insert” does in general not instigate the same gap in all other sequences. Since the subsequence of the insert contains information that is used by MAFFT in the MSA reconstruction, residues and or subsequences from other sequences will be aligned to the sequence-specific insert. As such, we define the majority gap column, where majority means > 90% (overule by p3) of the sequences. GIR determines which columns are majority gap columns, subsequently scans the MSA and computes the number of continuous regions of majority gap columns. Secondly, in complex cases one or more conserved residues may align somewhere in the instigated gap, thereby splitting it in two smaller gaps. GIR makes two counts. The continuous gap score resembles the local alignment score: it sums up as the majority gap continues, but drops to zero when a normal column is encountered. The combined gap score resembles the global alignment score: it sums and simply distracts each normal column that is encountered. The continuous gap score is used for the automated identification of GIR-NFH. The combined gap score is merely reported as part of the graphical output that is generated for each iteration and as is described in the main text. It can be used for the manual removal of additional sequences. GIR removes the sequence that has the gap with the highest continuous gap score above the threshold and saves and aligns the remaining sequences in order to iterate the procedure. By default GIR is a one-by-one (overrule m3) sequence remover with iteration.
4 The Continuous Gap Sequence Remover
CGSR removes sequences that have one or more large continuous gaps detected in the offered MSA. Since this gap can be the result of a family-specific deletion, gaps in high gap columns (>50%, overrule using p4) are not considered. Given the high complexity of MSAs with sequences that have these large continuous gaps, these MSAs tend to be highly irregular with many errors. In order to prevent removal of  sequences that are improperly aligned rather than erroneous, CGSR determines all gap sizes and selects sequences with gaps larger than 30 residues (overrule by parameter aa4) which size is in the upper 1.5 * IQR whisker (> Q3 + 1.5 IQR) of the gap-size distribution. This second rule does not apply to a last, additional iteration. CGSR separately saves removed and accepted sequences. The latter are aligned for iteration. Hence, CGSR is a controlled batch-sequence remover with iteration. For each iteration, a histogram with the distribution of the high gap columns indexes and the corresponding upper whisker is generated.
	As a final remark, similar to what can happen to a sequence-specific insert, sequence-specific deletions can result in the misalignment of one or more residues that flank the deletion. As a result, such deletions can appear as split which impedes their detection. This problem is not tackled since the solution similar to that was used for GIR (combined gap score) tends to identify many loop regions and leads to many false positives.
5 The Pseudogene Remover
The pseudogene remover is the same outlier remover as NHHR except that is is iterated.  PR comes with optional settings for the determination of the cut-off and provides a graph that, besides the score plot, contains a plot for the score-drops between two consecutive hits in the HMMER search (delta-score) as well as four putative cut-offs. Besides the default σ3, it presents the σ2 threshold, the upper 1.5*IQR whisker (> Q3 + 1.5 IQR) and the major score-drop Δmax as putative cut-off threshold. Application of the major score drop corresponds to the idea that non-outliers (read functional homologues under functional constraint) evolve with similar evolutionary rates which results in a continuous HMMER score distribution whereas NFHs lack the constraint and will evolve faster and show much lower scores. In certain cases, the major score drop can sometimes detect this. Note however that the dataset supposedly contains sequences from various subfamilies that have somewhat different constraints and that this can also result in large score-drops. The more complex a superfamily is, the less reliable it is to use the largest score-drop as cut-off. Hence, the largest score drop should only be applied when it occurs at or near one of the other three thresholds.
Example: 
python3 seqrutinator.py -f sample.fsa -m 1324 -p1 85 -m2 2 s2 2 -a2 5 -m3 2 -p3 80 -aa3 35 -p4 85
Runs seqrutinator in the order SSR-GIR-NHHR-CGSR, hence it will remove sequences <85% of the length of the first sequence in sample.fsa, it will then iteratively remove sequences that instigate gap regions of 35 (where a gap column has more than 80% gaps) in batch mode, next it removes sequences with gap regions of > 30 where 85% of the sequences are allowed to also have a gap at the same column. Do note that the latter is actually not a good idea since it can result in the exclusion of a sequence that has a large gap in the MSA where many sequences actually have gaps.

Table 2: Input and output generated by default Seqrutinator pipeline. GIR, CGSR and PR are principly iterated, shown are the results with two iterations. * If no sequences are removed by the last applied module, the MSA of the former module is used. ** Only the final dataset is saved. *** Only the last HMMER profile of PR iterations is saved. **** Not made when no sequences are removed. ***** Identical to the last MSA. X stands for the number of the last iteration of the former module.
Module & Iterate

Input*
Output



Identified as


Test or datafile
Content
Accepted/Removed
SSR
bug1.faa
SSR_bug1_data.txt
threshold and lengths of removed
SSR_bug1.fsa; .faa
SSR_bug1_removed.txt; .fsa
NHHR
SSR_bug1.faa; fsa
NHHR_bug1.hmm; NHHR_hmms_bug1.txt:
NHHR_histogram_bug1.png;
NHHR_bug1_seqs_lengths.txt
NHHR_bug1_seqs_removed.txt
hmmer profile
hmmsearch output
histogram of hmmsearch ouput
lengths of sequences
summary of removal
NHHR_bug1.fsa; faa
NHHR_bug1_removed.txt; .fsa
GIR I1
NHHR_bug1.faa*
GIR_it1_Data_bug1.txt; GIR_it1_Plot_bug1.png
gap scores iteration 1
plot of gap scores iteration 1
GIR_1_bug1_hits.fsa; faa
GIR I2
GIR_1_bug1.faa
GIR_it2_Data_bug1.txt; GIR_it2_Plot_bug1.png
gap scores iteration 2
plot of gap scores iteration 2
GIR_2_bug1_hits.fsa; faa
GIR_removed_bug1.txt; fsa**
CGSR I1
GIR_X-bug1.faa
CGSR_It1_bug1_gapdata.txt
CGSR_It1_bug1_gap_histo.png
CGSR_It1_bug1_data.txt
gaplenght data iteration 1
histogram of gaplengths iteration 1
summary of removal iteration 1
CGSR_It1_bug1.fsa; faa
CGSR I2
CGSR_1_bug1.faa
CGSR_It2_bug1_gapdata.txt
CGSR_It2_bug1_gap_histo.png
CGSR_It2_bug1_data.txt
gaplenght data iteration 2
histogram of gaplengths iteration 2
summary of removal iteration 2
CGSR_It2_bug1.fsa; .faa
CGSR_bug1_removed.txt; fsa
PR I1
CGSR_X_bug1.faa
PR_bug1_It1.hmm***
PR_bug1_data_It1.txt
PR_bug1_plot_It1.png
hmmer profile
hmmsearch output
Elaborate score plot
PR_It1_bug1.fsa; faa****
PR I2
PR_1_bug1.faa
PR_bug1_data_2.txt

PR_It2_bug1.fsa; faa
PR_bug1_removed.txt; .fsa




Result_big1.faa*****

At the end of each module seqrutinator performs a BMGE test that determines the number of columns before and after trimming and saves the numbers in a tsv file. The rationale of this is explained in the main text of the manuscript. 

           SeqYNet.py
The additional script SeqYNet.py can be used to run Seqrutinator with multiple input files at once. Both scripts, the library and the set of fasta files, each prepared as explained, are to be placed in the same folder. All input files must have the extension fasta for running: The SeqYNet script comes with the same options and results for each fasta file will be saved in different subfolders, named on the base of the input fasta file.

Materials and Methods
Alignments. MUFASA and Seqrutinator use by default MAFFT G-INS-i 
mafft --reorder --maxiterate 1000 --retree 1 --globalpair <input.fsa> <output.faa>
MAFFT [52]⁠ can be installed to most Linux systems using Synaptic or is otherwise available at the MAFFT website [68]⁠.
Note that MAFFT G-INS-i comes with a heavy computational cost. For very large datasets (e.g. > 500 sequences of ~400 aa) we recommend to use FAMSA [25]⁠, which is an available option. FAMSA should be made executable in the same folder. FAMSA can be obtained from its Github repository [69]⁠.
HMMER The outlier removers use HMMER [34]⁠ with the following commands:
	hmmbuild --informat afa --amino <output.hmm> <input.faa>		
	hmmsearch --noali –tblout <output.txt> <input.fsa>		
Block Mapping and Gathering of Entropy [32]⁠. The pipeline and modules execute the command 
	java -jar BMGE.jar -i <msa.faa> -h 0.8						
Note that the command uses standard settings except for the entropy that is set to 0.8. Note that the actual trimming is not performed. We recommend BMGE version 1.0 that can be obtained from the Pasteur Instititute [70]⁠. The jar-file should be in the same directory as where Seqrutinator is executed. 
Dependencies of seqrutinator are indicated in the script.
