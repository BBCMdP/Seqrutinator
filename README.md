# Seqrutinator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10980626.svg)](https://doi.org/10.5281/zenodo.10980626)

- [Seqrutinator](#seqrutinator)
  - [Introduction](#introduction)
  - [MUFASA (MUltiple FASta Aligner)](#mufasa-multiple-fasta-aligner)
    - [Usage](#usage)
    - [Options](#options)
  - [Seqrutinator](#seqrutinator-1)
    - [Usage](#usage-1)
    - [Options](#options-1)
    - [Modules](#modules)
      - [Module 1: The Short Sequence Remover Module](#module-1-the-short-sequence-remover-module)
      - [Module 2: The Non-Homologous Hit Remover Module](#module-2-the-non-homologous-hit-remover-module)
      - [Module 3: The Gap Instigator Remover](#module-3-the-gap-instigator-remover)
      - [Module 4: The Continuous Gap Sequence Remover](#module-4-the-continuous-gap-sequence-remover)
      - [Module 5: The Pseudogene Remover](#module-5-the-pseudogene-remover)
  - [SeqYNet](#seqynet)
    - [Materials, Methods and Requirements](#materials-methods-and-requirements)
  - [Auxiliary scripts](#auxiliary-scripts)
    - [Sequence renamer (seq\_renamer)](#sequence-renamer-seq_renamer)
      - [Usage](#usage-2)
      - [Options](#options-2)
      - [Requirements](#requirements)
    - [Sequence Non-IUPAC Purge script (SNIP)](#sequence-non-iupac-purge-script-snip)
      - [Introduction](#introduction-1)
      - [Requirements](#requirements-1)
      - [Single fasta file](#single-fasta-file)
      - [Multiple fasta files](#multiple-fasta-files)
  - [References](#references)

Seqrutinator is an objective, flexible pipeline for the scrutiny of sequence sets from complex, eukaryotic protein superfamilies. It removes sequences from pseudogenes, incorrect gene models or with sequencing errors. Testing Seqrutinator on major superfamilies BAHD, CYP and UGT correctly removed 1.94% of SwissProt entries, 14% of entries from the model plant _Arabidopsis thaliana_ but 80% of entries from _Pinus taeda_’s recent complete proteome. Most of the removed sequences were partial. The scrutiny of BAHDomes, CYPomes and UGTomes obtained from 16 plant proteomes show consistent numbers suggesting good performance. MSAs, phylogenies and particularly functional clustering improved drastically upon Seqrutinator application.

Seqrutinator is under review but can be accessed as preprint at https://doi.org/10.1101/2022.03.22.485366

Below you will find the contents of Supplemental document 1, which is basically a detailed descripion of Seqrutinator and accompaying scripts. My apologies for the fact that the output tables are no really tables here

## Introduction
Seqrutinator is a modular pipeline written in Python 3 to identify and remove protein superfamily sequences that likely either correspond to a non-functional homologue (NFH) or are incorrect, the latter due to either sequencing or gene modelling errors. Seqrutinator is made for complex single-domain protein superfamilies. It requires intermediate or large datasets (~>50 input sequences). This document describes in detail how the modules function, how the pipeline is made and which outputs are generated. Seqrutinator requires a single input file with the sequence set in fasta format, either aligned or not. Dependencies are described/listed in Materials and Methods, which includes shell commands that are used to run software packages such as MAFFT.

Each module is accessed by a number, according to the default pipeline 12345: 1) Short Sequence Remover (SSR) → 2) Non-Homologous Hit Remover (NHHR) → 3) Gap Instigator Remover (GIR) → 4) Continuous Gap Sequence Remover → 5) Pseudogene Remover (PR). Modules 2 to 5 need an MSA as input but when applied as first module, the script automatically aligns the sequences if they are not aligned. We recommend the user trims the MSA such that it represents the mature protein and lacks non-homologous terminal sequences.
Seqrutinator is therefore accompanied by MUFASA (MUltiple FASta Aligner) which is used to obtain the initial crude but aligned sequence set. For trimming and in order to include a general reference we recommend aligning a single sequence for which a structure has been resolved to the MSA by means of `mafft --add`. Seqrutinator can be applied in batch (i.e. using multiple input files) using an additional script named SeqYNet. All three scripts are described below. Descriptions concern the default settings but we describe many options to overrule default settings.

## MUFASA (MUltiple FASta Aligner)

MUFASA performs hmmsearch and aligns positively identified sequences. It only comes with options regarding input files as demonstrated by:
```bash
$ python3 MUFASA.py -h
```

### Usage
```bash
MUFASA.py [-h] [-i I] [-ext EXT] [-c C] [-r R] [-t T]
```

### Options
```bash
  -h, --help  show this help message and exit
  -i I        HMMER profile, default = hmm_profile
  -ext EXT    Target files extension. default = *.fsa
  -c C        Cores, default = 4
  -r R        Reference sequence for mafft-add and trimming (must be in fasta format)
  -t T        Trimming N- and C-end based on reference's positions i-j (residues i and j will be conserved)
```

MUFASA takes as input a single HMMER profile and performs separate hmmsearch of all fasta files available in the same folder as where MUFASA is executed. By default, all files with extension `fsa` will be used as input sequence set and `hmm_profile` is the default HMMER profile. For each input sequence set, it fetches all identified sequences and saves them in a single fasta file, which represent a superfamily proteome. Next, all obtained superfamily proteomes are separately aligned using MAFFT G-INS-i.
Example:
```bash
$ python3 MUFASA.py -i PF02458.hmm -ext '*.fsa'
```
It will use the hmmer profile PF02458.hmm to run hmmsearch using on each file with extension `.fsa`. The outputs provided by this module are shown in Table 1. Note that each type of output will be stored in a dedicated folder.

In addition, since typically N- and C-end regions tend to be more variable (e.g., signal peptides), we recommend trimming these regions before to run Seqrutinator.
It is highly recommended to use the same reference sequence, which should be well characterized biochemically and physiologically, and which structure has been resolved (highly recommended).
MUFASA includes a module that can do this. This module is activated by incorporate one reference sequence in fasta format (with a different extension as the proteomes explored in the first module) using the argument `-r`. This reference sequence will be added to the MSA with all the homologues detected by MUFASA in each proteome (using MAFFT-add).
Finally, following the coordinates given with argument `-t i-j`, the script will trim the MSA based on residues i and j from the reference sequence (columns upstream i and downstream j will be removed). The results after MAFFT-add and after trimming are stored (see Table 1).

Example:
```bash
$ python3 MUFASA.py -i PF02458.hmm -ext '*.fsa' -r 4G0B.fa -t 4-432
```
The module MUFASA will operate as described before. Since `-r` is called, the sequence in fasta file `4G0B` will be added to each `*.faa` file (generated by module MUFASA) using MAFFT-add, yielding `*_add.faa` files. The latter will be used to trim N- and C-ends using positions 4 and 432 from the reference sequence, resulting in `*._trim.faa` files.

![image](https://github.com/BBCMdP/Seqrutinator/assets/45858786/680b7a3b-962b-49cf-8b8c-36c397eca569)
>Table 1: Input and Output of MUFASA

## Seqrutinator
Seqrutinator is a fully flexible, modular pipeline. It consists of a main python script with the actual modules and some recurring algorithms in a python module (lib_seqrutinator.py). Seqrutinator comes with a large number of options. Parameters have numbers according to the number of the module, Parameters that are not numbered are general. Parameters are as demonstrated by:

```bash
$ python3 seqrutinator.py -h
```
### Usage
```bash
seqrutinator.py [-h] [-m M] [-f F] [-ali ALI] [-ref1 REF1] [-ref2 REF2] [-bv BV] [-BMGE BMGE] [-p1 P1] [-p2 P2]
                       [-s2 S2] [-a2 A2] [-m3 M3] [-p3 P3] [-aa3 AA3] [-p4 P4] [-aa4 AA4] [-a5 A5] [-s5 S5]
```

### Options
```bash
  -h, --help  show this help message and exit
  -m M        Pipelines
  -f F        Fasta file (NOTE THAT THIS IS A REQUIRED PARAMETER)
  -ali ALI    Use MAFFT G-INS-i (1); Use either MAFFT G-INS-i (n<=500) or Global (n>500) (2); Use FAMSA, recommended only for really large datasets (3)
  -ref1 REF1  Use to change input reference sequence
  -ref2 REF2  Use to directly input reference sequence length
  -bv BV      BMGE version 1 (either 1.0 or 1.12) or 2
  -BMGE BMGE  BMGE deactivated (0), BMGE > 0 is h option for BMGE
  -p1 P1      Proportion (0 to 1) of sequence length coverage for SSR
  -p2 P2      Proportion (0 to 1) of sequence length coverage for NHHR
  -s2 S2      Mean - alphaSD (1) or Q1 - 1.5IQR (2) for NHHR
  -a2 A2      Alpha for NHHR
  -m3 M3      Method one by one (1) or batch (2) for GIR
  -p3 P3      Proportion (0 to 1) of gaps to define a gap column for GIR (>= VALUE)
  -aa3 AA3    aa window of contiguos gap columns for GIR
  -p4 P4      Proportion (0 to 1) of gaps to define a gap column for CGSR (>= VALUE)
  -aa4 AA4    aa window of contiguos gap columns for CGSR
  -a5 A5      Alpha for mean - alphaSD (3 is recommended as default option and 2.35 for normal distributions)
  -s5 S5      Mean - alphaSD (1) or Q1 - 1.5IQR (2) for PR
```

### Modules
#### Module 1: The Short Sequence Remover Module
SSR opens the offered fasta file and determines the length l of the first sequence,  the default reference sequence. This can be overruled by `ref1` or `ref2`. Next, it determines the length of each sequence and separately saves good and bad sequences using the default threshold of `0.65*l`, which can be overruled by `p1`. Good sequences are aligned by MAFFT G-INS-i. The additional output consists of a file that contains the accession codes of the removed sequences and a file that shows the inclusion threshold and the lengths of the removed sequences.

#### Module 2: The Non-Homologous Hit Remover Module
NHHR removes outliers identified by low HMMER score. It constructs a HMMER profile of the input MSA and uses it to run a hmmsearch against the unaligned sequences. It directly removes all sequences with score=0 and then uses the total scores of the sequences to perform a distribution analysis. In order to prevent removing short sequences (in case SSR does not precede NHHR in the applied pipe), sequences shorter than 65% (overrule by `p2`) of the reference sequence are excluded from this analysis. It determines the mean and the standard deviation of all the data in the array.
Finally, it determines which sequences score below and which above the threshold and saves removed and accepted sequences separately (the latter are subsequently aligned). NHHR is not iterated since its objective is to remove sequences that are not homologous and will be removed in a single screening, under the assumption that they score low. The additional output consists a file that contains the accession codes of the removed sequences, the hmmsearch ouputs and a plot that shows a distribution of the calculated scores with the applied cut-off threshold.
Threshold settings can be overruled. By default, it uses 1.5 * Inter Quartile Range (IQR) whisker (< `Q1 – 1.5 x IQR`). The Inter Quartile Range or IQR is defined as the difference between the middle of the first half of the sequence scores and the middle of the second half of the sequence scores. Using argument `-s2`, the threshold can be modified to use the average score - alpha * st. dev. The alpha value is set to 3 by default (σ3), but can be modified using argument `-a2`. 

#### Module 3: The Gap Instigator Remover
GIR removes sequences that instigate large (>29, overrule using `aa3`) gaps as detected in the offered MSA. This comes with two uncertainties of which only the firts one is tackled automatically. A “sequence-specific insert” does in general not instigate the same gap in all other sequences. Since the subsequence of the insert contains information that is used by MAFFT in the MSA reconstruction, residues and or subsequences from other sequences will be aligned to the sequence-specific insert. As such, we define the majority gap column, where majority means > 90% (overule by `p3`) of the sequences. GIR determines which columns are majority gap columns, subsequently scans the MSA and computes the number of continuous regions of majority gap columns. Secondly, in complex cases one or more conserved residues may align somewhere in the instigated gap, thereby splitting it in two smaller gaps. GIR makes two counts. The continuous gap score resembles the local alignment score: it sums up as the majority gap continues, but drops to zero when a normal column is encountered. The combined gap score resembles the global alignment score: it sums and simply distracts each normal column that is encountered. The continuous gap score is used for the automated identification of GIR-NFH. The combined gap score is merely reported as part of the graphical output that is generated for each iteration and as is described in the main text. It can be used for the manual removal of additional sequences. GIR removes the sequence that has the gap with the highest continuous gap score above the threshold and saves and aligns the remaining sequences in order to iterate the procedure. By default GIR is a one-by-one (overrule `m3`) sequence remover with iteration.

#### Module 4: The Continuous Gap Sequence Remover
CGSR removes sequences that have one or more large continuous gaps detected in the offered MSA. Since this gap can be the result of a family-specific deletion, gaps in high gap columns (>50%, overrule using `p4`) are not considered. Given the high complexity of MSAs with sequences that have these large continuous gaps, these MSAs tend to be highly irregular with many errors. In order to prevent removal of  sequences that are improperly aligned rather than erroneous, CGSR determines all gap sizes and selects sequences with gaps larger than 30 residues (overrule by parameter `aa4`) which size is in the upper 1.5 * IQR whisker (> `Q3 + 1.5 IQR`) of the gap-size distribution. This second rule does not apply to a last, additional iteration. CGSR separately saves removed and accepted sequences. The latter are aligned for iteration. Hence, CGSR is a controlled batch-sequence remover with iteration. For each iteration, a histogram with the distribution of the high gap columns indexes and the corresponding upper whisker is generated.
As a final remark, similar to what can happen to a sequence-specific insert, sequence-specific deletions can result in the misalignment of one or more residues that flank the deletion. As a result, such deletions can appear as split which impedes their detection. This problem is not tackled since the solution similar to that was used for GIR (combined gap score) tends to identify many loop regions and leads to many false positives.

#### Module 5: The Pseudogene Remover
The pseudogene remover is the same outlier remover as NHHR except that it is iterated. PR comes with optional settings for the determination of the cut-off and provides a graph that, besides the score plot, contains a plot for the score-drops between two consecutive hits in the HMMER search (delta-score) as well as four putative cut-offs. Besides the default 1.5 * Inter Quartile Range (IQR) whisker (< `Q1 – 1.5 x IQR`), it can also use average - alpha * st. dev. (using `-s5`, σ3 by default) with a customizable alpha value (`-a5`)
**Example:**
```bash
python3 seqrutinator.py -f sample.fsa -m 1324 -p1 0.85 -m2 2 s2 2 -a2 5 -m3 2 -p3 0.8 -aa3 35 -p4 0.85
```
>Runs seqrutinator in the order SSR-GIR-NHHR-CGSR, hence it will remove sequences <85% of the length of the first sequence in sample.fsa, it will then iteratively remove sequences that instigate gap regions of 35 (where a gap column has more than 80% gaps) in batch mode. Next it removes sequences with gap regions of > 30 where 85% of the sequences are allowed to also have a gap at the same column. Do note that the latter is actually not a good idea since it can result in the exclusion of a sequence that has a large gap in the MSA where many sequences actually have gaps.


| Module & Iterate | Input*            | Output                                                                                                                                                                                                                                                                                |||
|------------------|-------------------|-----------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|-------------------------------------------------------|
|                  |                   |                                                                                                                       | Identified as                                                                                                                                                  ||
|                  |                   | Test or datafile                                                                                                      | Content                                                                                                 | Accepted/Removed                                      |
| SSR              | bug1.faa          | SSR_bug1_data.txt                                                                                                     | threshold and lengths of removed                                                                        | SSR_bug1.fsa; .faa  SSR_bug1_removed.txt; .fsa        |
| NHHR             | SSR_bug1.faa; fsa | NHHR_bug1.hmm;  NHHR_hmms_bug1.txt:  NHHR_histogram_bug1.png;  NHHR_bug1_seqs_lengths.txt  NHHR_bug1_seqs_removed.txt | hmmer profile  hmmsearch output  histogram of hmmsearch ouput  lengths of sequences  summary of removal | NHHR_bug1.fsa; faa  NHHR_bug1_removed.txt; .fsa       |
| GIR I1           | NHHR_bug1.faa*    | GIR_it1_Data_bug1.txt;  GIR_it1_Plot_bug1.png                                                                         | gap scores iteration 1  plot of gap scores iteration 1                                                  | GIR_1_bug1_hits.fsa; faa                              |
| GIR I2           | GIR_1_bug1.faa    | GIR_it2_Data_bug1.txt;   GIR_it2_Plot_bug1.png                                                                        | gap scores iteration 2  plot of gap scores iteration 2                                                  | GIR_2_bug1_hits.fsa; faa  GIR_removed_bug1.txt; fsa** |
| CGSR I1          | GIR_X-bug1.faa    | CGSR_It1_bug1_gapdata.txt  CGSR_It1_bug1_gap_histo.png  CGSR_It1_bug1_data.txt                                        | gaplenght data iteration 1  histogram of gaplengths iteration 1  summary of removal iteration 1         | CGSR_It1_bug1.fsa; faa                                |
| CGSR I2          | CGSR_1_bug1.faa   | CGSR_It2_bug1_gapdata.txt  CGSR_It2_bug1_gap_histo.png  CGSR_It2_bug1_data.txt                                        | gaplenght data iteration 2  histogram of gaplengths iteration 2  summary of removal iteration 2         | CGSR_It2_bug1.fsa; .faa  CGSR_bug1_removed.txt; fsa   |
| PR I1            | CGSR_X_bug1.faa   | PR_bug1_It1.hmm***  PR_bug1_data_It1.txt  PR_bug1_plot_It1.png                                                        | hmmer profile  hmmsearch output  Elaborate score plot                                                   | PR_It1_bug1.fsa; faa****                              |
| PR I2            | PR_1_bug1.faa     | PR_bug1_data_2.txt                                                                                                    |                                                                                                         | PR_It2_bug1.fsa; faa  PR_bug1_removed.txt; .fsa       |
|                  |                   |                                                                                                                       |                                                                                                         | Result_big1.faa*****                                  |
>Table 2: Input and output generated by default Seqrutinator pipeline. GIR, CGSR and PR are principly iterated, shown are the results with two iterations.
*If no sequences are removed by the last applied module, the MSA of the former module is used.
**Only the final dataset is saved
***Only the last HMMER profile of PR iterations is saved
****Not made when no sequences are removed
*****Identical to the last MSA. X stands for the number of the last iteration of the former module.


At the end of each module seqrutinator performs a BMGE test that determines the number of columns before and after trimming and saves the numbers in a tsv file. The rationale of this is explained in the main text of the manuscript.

## SeqYNet

The additional script SeqYNet.py can be used to run Seqrutinator with multiple input files at once. Both scripts, the library and the set of fasta files, each prepared as explained, are to be placed in the same folder. All input files must have the extension fasta for running: The SeqYNet script comes with the same options and results for each fasta file will be saved in different subfolders, named on the base of the input fasta file.

---
### Materials, Methods and Requirements

**Alignments**. MUFASA and Seqrutinator use by default MAFFT G-INS-i:
`mafft --reorder --maxiterate 1000 --retree 1 --globalpair <input.fsa> <output.faa>`
MAFFT[^1]⁠ can be installed to most Linux systems using Synaptic or is otherwise available at the MAFFT website[^2]⁠.
Note that MAFFT G-INS-i comes with a heavy computational cost. For very large datasets (e.g. > 500 sequences of ~400 aa) we recommend to use FAMSA[^3], which is an available option. FAMSA should be made executable in the same folder. FAMSA can be obtained from its Github repository[^4]⁠.

**Block Mapping and Gathering of Entropy**[^5]⁠. The pipeline and modules execute the command
`java -jar BMGE.jar -i <msa.faa> -h 0.8`
Make sure the executable file BMGE.jar. By default, Seqrutinator assumes BMGE version used is 1 (either v 1.0 or 1.12).
Note that the command uses standard settings except for the entropy that is set to 0.8. Note that the actual trimming is not performed. We recommend BMGE version 1 that can be obtained from the Pasteur Institute  [^6]⁠. The jar-file should be in the same directory as where Seqrutinator is execute.
Note, however, that version 2.0 [^7] can also be used (set the argument BMGE version `-bv 2` if this is the case).

**HMMER** [^8]. Installation via Synaptic is available. The outlier removers use HMMER modules⁠ with the following commands:
`hmmbuild --informat afa --amino <output.hmm> <input.faa>`
`hmmsearch --noali –tblout <output.txt> <input.fsa>`

**Python packages**. The `requirements.txt` file can be used to install the required python3 packages. Run for example:

```bash
$ pip install -r requirements.txt
```


## Auxiliary scripts

### Sequence renamer (seq_renamer)

Available at https://github.com/BBCMdP/Sequence-renamer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10949335.svg)](https://doi.org/10.5281/zenodo.10949335)


Sequences in fasta file may often times contain characters that are incompatible not only with Seqrutinator's dependencies (as whitespaces for hmmer), but also for downstream applications (e.g., dots are replaced by underscores by PhyML). In addition, if names are long can be chopped to a short string by other programs (such as cd-hit). This name format hell makes analysis of sequences data generated by different softwares quite challenging.
To overcome this, here we provide the script `seq_renamer.py` which allows the renaming of all sequences in a fasta file to a name with 10 characters. The user can add a text string that enables the quick identification of the sequences by the user. We recommend a string of 3 characters, which will be extended to the final 10 characters by an increasing integer. The outputs are a fasta file with the renamed sequences, and a file with an index matching the new names with the original names.

```bash
$ python3 seq_renamer.py -h
```
#### Usage
```bash
seq_renamer.py [-h] [-i I] [-id ID]
```

#### Options
```bash
  -h, --help  show this help message and exit
  -i I        Fasta file
  -id ID      ID Key sequence name
```

```bash
$ python3 seq_renamer.py -i Ecoli.fasta -id eco
Execution Successful: 0:00:00.033954
```
yields two files: a fasta file with the renamed sequences `Ecoli_renamed.fsa` and a text file `Ecoli_file_map`, which contains an index of new and original names.
An extract of the results:
```bash
$ head Ecoli_renamed.fsa -n 8
>eco0000001
MTHIVRFIGLLLLNASSLRGRRVSGIQH
>eco0000002
MFLDYFALGVLIFVFLVIFYGIIILHDIPYLIAKKRNHPHADAIHVAGWVSLFTLHVIWPFLWIWATLYRPERGWGMQSHDSSVMQLQQRIAGLEKQLADIKSSSAE
>eco0000003
MHAYLHCLSHSPLVGYVDPAQEVLDEVNGVIASARERIAAFSPELVVLFAPDHYNGFFYDVMPPFCLGVGATAIGDFGSAAGELPVPVELAEACAHAVMKSGIDLAVSYCMQVDHGFAQPLEFLLGGLDKVPVLPVFINGVATPLPGFQRTRMLGEAIGRFTSTLNKRVLFLGSGGLSHQPPVPELAKADAHMRDRLLGSGKDLPASERELRQQRVISAAEKFVEDQRTLHPLNPIWDNQFMTLLEQGRIQELDAVSNEELSAIAGKSTHEIKTWVAAFAAISAFGNWRSEGRYYRPIPEWIAGFGSLSARTEN
>eco0000004
MYYLKNTNFWMFGLFFFFYFFIMGAYFPFFPIWLHDINHISKSDTGIIFAAISLFSLLFQPLFGLLSDKLGLRKYLLWIITGMLVMFAPFFIFIFGPLLQYNILVGSIVGGIYLGFCFNAGAPAVEAFIEKVSRRSNFEFGRARMFGCVGWALCASIVGIMFTINNQFVFWLGSGCALILAVLLFFAKTDAPSSATVANAVGANHSAFSLKLALELFRQPKLWFLSLYVIGVSCTYDVFDQQFANFFTSFFATGEQGTRVFGYVTTMGELLNASIMFFAPLIINRIGGKNALLLAGTIMSVRIIGSSFATSALEVVILKTLHMFEVPFLLVGCFKYITSQFEVRFSATIYLVCFCFFKQLAMIFMSVLAGNMYESIGFQGAYLVLGLVALGFTLISVFTLSGPGPLSLLRRQVNEVA
```

```bash
$ head Ecoli_file_map -n 4
eco0000001 	 WP_001300467.1 MULTISPECIES: leu operon leader peptide [Enterobacteriaceae]
eco0000002 	 WP_000478195.1 MULTISPECIES: DUF3302 domain-containing protein [Enterobacteriaceae]
eco0000003 	 WP_000543457.1 MULTISPECIES: 2,3-dihydroxyphenylpropionate/2,3-dihydroxicinnamic acid 1,2-dioxygenase [Bacteria]
eco0000004 	 WP_000291549.1 MULTISPECIES: lactose permease [Bacteria]
```
We recommend running the script before running the sensitive search with `MUFASA`. 
#### Requirements
`seq_renamer.py` requires Biopython.

### Sequence Non-IUPAC Purge script (SNIP)

[![DOI](https://zenodo.org/badge/788514656.svg)](https://zenodo.org/doi/10.5281/zenodo.10994866)

Available at https://github.com/BBCMdP/SNIP

#### Introduction 

A common issue derived from poorly annotated genomes is that resulting proteomes may come with characters that do not correspond to actual amino acids commonly found in natural proteins (what we define here as non-IUPAC characters). 
The most common scenario is the low definition of bases in the genome, typically annotated as nucleotide "N", for which translation of the predicted coding sequence may lead to an undefined residue, often times represented with symbol "X" in the protein sequence. 

Some downstream applications (even MSA generation, profile to sequence comparisons, and even phylogeny reconstruction) can't deal with for non-IUPAC characters, and results are compromised. Moreover, some proteomes show a high "contamination" with non-IUPAC annotations, making this a problem that can easily become a nuissance. In addition, stop codon character * is typically found in complete proteome annotations, and can also interfere with several downstream applications.  

In order to purge the protein datasets from sequences with non-IUPAC annotations, we have developed the Sequence Non-IUPAC Purge script, SNIP. SNIP is a python script that can deal with single or multiple fasta files, to automatically detect and remove sequences with non-IUPAC characters. SNIP do accept gap characters (such as "-" or "."), so either aligned or unaligned fasta files are accepted. The script will also process terminal stop codon character "*", and specifically remove it from the sequence. Note, however, that if a non terminal stop codon character is found in a sequence, it will be considered as non-IUPAC character, resulting in the removal of the sequence.

>**We recommend applying this script before to run MuFasA and Seqrutinator.** 

#### Requirements

SNIP requires Biopython, and accepts by argument either single (`-s`) or multiple (`-m`) fasta files. For multiple file, it is required that all of them have the same extension. 

#### Single fasta file
The output will depend on the results. If running a single fasta file:

`python3 SNIP.py -s input.fasta` 

and sequences with non-IUPAC characters are found, two output fasta files will be generated: `input_accepted.fasta` and `input_removed.fasta`. Of course, each file includes the sequences without and with non-IUPAC characters, respectively. The terminal will print out the amount of sequences (total, accepted and removed). It will also show which is the non-IUPAC character found for each removed sequence.

```
Results for file input.fasta
Total seqs:  65809
Accepted seqs:  65368
Removed seqs:  441
```
Note, in addition, that sequences in both files will not have terminal stop codon character "*".

Now, if no sequences with non-IUPAC characters are found, but terminal stop codons are detected, the resulting sequences, without these characters, will be written in a new file with extension `_nsc.fasta` (from no stop codons).

Finally, if neither non-IUPAC or terminal stop codon characters are found, only a printout in terminal is shown expliciting it (no further files are written). 

#### Multiple fasta files
SNIP can be applied to multiple fasta files. The only requirements are that all files should (i) be in the same folder, and (ii) have the same extension. For example, can be run like this: 

`python3 SNIP.py -m '*.fa`  

in a folder with five fasta files with the provided extension:

```
SNIP.py
sfa.fa
crh.fa
smu.fa
tpl.fa
cri.fa
```

Will produce two folders, `/Seqs_Accepted` and `/Seqs_Removed`, plus a summary file `Report_multiple_files.tsv`. The files generated are the same and created conditionally as described before. The difference is that all files `*_removed.fa` will be moved to the `/Seqs_Removed` folder (if any). The `*_accepted.fa` and `*_ncs.fa` (if any) are moved to `/Seqs_Accepted`, as well as all files that were not modified by the script. This is to simplify the output recovery for the user. 

The `Report_multiple_files.tsv` summarizes the results, for example:

| File           | Total_Seqs | Seqs_accepted | Seqs_removed | term_* |
|----------------|------------|---------------|--------------|--------|
| crh.fa         | 19527      | 19514         | 13           | yes    |
| cri.fa | 75253      | 75253         | 0            | yes    |
| sfa.fa | 45611      | 45611         | 0            | no     |
| smu.fa | 27137      | 24567         | 2570         | no     |
| tpl.fa | 65809      | 65368         | 441          | yes    |
> The column term_* indicates if terminal stop codon characters were found

The user can identify here that for input file `crh.fa` has 19527 seqs, of which 13 have non-IUPAC characters (and are removed), and terminal stop codon characters were found. Thus, `crh_accepted.fa` and `crh_removed.fa` files are created and saved in the proper folder. The same goes for `smu.fa` (note that around 10% of the sequences have non-IUPAC characters) and `tpl.fa`. 
The case of `cri.fa` will result in a file `cri_nsc.fa`, since there were no sequences with non-IUPAC characters, but terminal stop codons were found. 
Finally, a copy of `sfa.fa` will be found in the `/Seqs_Accepted` folder, since neither non-IUPAC or terminal stop codons were found. 


  
## References
[^1]: K. Katoh and D. M. Standley, “MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability,” Mol. Biol. Evol., vol. 30, no. 4, pp. 772–780, Apr. 2013, doi: 10.1093/molbev/mst010.
[^2]: https://mafft.cbrc.jp/alignment/server/.
[^3]: S. Deorowicz, A. Debudaj-Grabysz, and A. Gudys, “FAMSA: Fast and accurate multiple sequence alignment of huge protein families,” Sci. Rep., vol. 6, no. 1, pp. 1–13, Sep. 2016, doi: 10.1038/srep33964.
[^4]: https://github.com/refresh-bio/FAMSA/.
[^5]: A. Criscuolo and S. Gribaldo, “BMGE (Block Mapping and Gathering with Entropy): A new software for selection of phylogenetic informative regions from multiple sequence alignments,” BMC Evol. Biol., vol. 10, no. 1, pp. 1–21, Jul. 2010, doi: 10.1186/1471-2148-10-210/FIGURES/9.
[^6]: ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/
[^7]: https://gitlab.pasteur.fr/GIPhy/BMGE
[^8]: S. R. Eddy, “Accelerated Profile HMM Searches,” PLOS Comput. Biol., vol. 7, no. 10, p. e1002195, 2011, doi: 10.1371/JOURNAL.PCBI.1002195.

