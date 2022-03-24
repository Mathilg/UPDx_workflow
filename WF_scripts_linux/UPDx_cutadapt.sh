#!/bin/bash

########################################################################
########################################################################
## PRIMER TRIMMING of the raw reads with cutadapt


# system time and date at start
date

###################################################
##### Inpu/output
##################################################

# input file:
    # samplename_L001_R1_001.fastq.gz;
    # samplename_L001_R2_001.fastq.gz
# output:
    # sample.txt
    # samplename_L001_R2_001_trim1.fastq;
    # samplename_L001_R1_001_trim1.fastq;
    # cutadapt_primer_trimming_stats1.txt

###################################################
##### Package/program Requirement
##################################################

# CUTADAPT
# advice for installation: https://cutadapt.readthedocs.io/en/stable/installation.html#quick-installation
# Here, installation of cutadapt via conda -->  "conda install -c bioconda cutadapt"

###################################################
##### Removing primers with cutadapt
##################################################

# nUPDx = paired reads 2*251bp : reads length superior to target amplicon length (180bp) == primer sequences present at both extremities of the reads (5' and 3')
  # TRIMMING IN 5' and 3'
# nUPDx R1 and R2 = Both R1 and R2 have a mix of forward and reverse reads == trimming of the primers in both direction for both files
  # TRIMMING IN 5' and 3' in both R1 and R2 files

# nUPDx primers:
  #Primer_F="CCGGAGAGGGAGCCTGAGA"
  #Primer_Frc="TCTCAGGCTCCCTCTCCGG"
  #Primer_R="GAGCTGGAATTACCGCGG"
  #Primer_Rrc="CCGCGGTAATTCCAGCTC"

# Retrieving sample name
# ex file name: 1-Pf-full_S1_L001_R2_001.fastq and 1-Pf-full_S1_L001_R1_001.fastq--> here sample name is "1-Pf-full_S1"
ls *_L001_R2_001.fastq* | cut -f1,2 -d "_" > samples.txt

#Create output file
mkdir 2_Cutadapt_Output

#run cutadapt
for sample in $(cat samples.txt)

do

echo "On sample: $sample"
cutadapt -a ^CCGGAGAGGGAGCCTGAGA...CCGCGGTAATTCCAGCTC \
-a ^GAGCTGGAATTACCGCGG...TCTCAGGCTCCCTCTCCGG \
-A ^GAGCTGGAATTACCGCGG...TCTCAGGCTCCCTCTCCGG \
-A ^CCGGAGAGGGAGCCTGAGA...CCGCGGTAATTCCAGCTC \
-m 145 -M 200 --discard-untrimmed \
-o ./2_Cutadapt_Output/${sample}_L001_R1_001_trim1.fastq -p ./2_Cutadapt_Output/${sample}_L001_R2_001_trim1.fastq \
${sample}_L001_R1_001.fastq* ${sample}_L001_R2_001.fastq* \
>> ./2_Cutadapt_Output/cutadapt_primer_trimming_stats1.txt 2>&1

done

##file summary according to the sample how many reads were trimmed
#Here’s a little one-liner to look at what fraction of reads were retained in each sample (column 2) and what fraction of bps were retained in each sample (column 3):
paste samples.txt <(grep "passing" ./2_Cutadapt_Output/cutadapt_primer_trimming_stats1.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" ./2_Cutadapt_Output/cutadapt_primer_trimming_stats1.txt | cut -f3 -d "(" | tr -d ")")


# system time and date at end
date


## Cutadapt parameters explanation
# -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC : trim primer in 5' and 3' from file R1 (primer F and primer Rrc)
# -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC : trim primer in 5' and 3' from file R2 (primer R and primer Frc)
# -m 145: remove reads with length < 200
# -M 200 : remove reads with length > 200
# --discard-untrimmed : removed trimmed reads
# -o : R1 output
# -p : R2 output
# input: R1 R2
# last line allow to print a summary of the trimming : >> cutadapt_primer_trimming_stats.txt 2>&1
