#!/bin/bash


########################################################################
########################################################################
## CONTROL CHECKING fastQC after trimming


# system time and date at start
date

###################################################
##### SETTINGS
##################################################


# input file:
    # samplename_L001_R1_001_trim1.fastq
    # samplename_L001_R2_001_trim1.fastq
# output:
    # samplename_L001_R1_001_trim1_fastqc.htlm;
    # samplename_L001_R1_001_trim1_fastqc.zip;
    # samplename_L001_R2_001_trim1_fastqc.htlm;
    # samplename_L001_R2_001_trim1_fastqc.zip;

###################################################
##### Package/program Requirement
##################################################

# fastQC : FastQC A quality control application for high throughput sequence data
# advice for installation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# if conda is installed: "conda install -c bioconda fastqc"

###################################################
##### Run FastQC on each fastq file on the command line :
##################################################

fastqc *_trim1.fastq
mkdir ../3_QC_after_trim
mv *_fastqc.* ../3_QC_after_trim

## view htlm file with safari:
# "open -a "Safari" PathToFile/samplename_L001_R1_001_trim1_fastqc.html"


# system time and date at end
date
