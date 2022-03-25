#!/bin/bash


########################################################################
########################################################################
## CONTROL CHECKING fastQC before trimming


# system time and date at start
date

###################################################
##### input/output
##################################################

# input file:
    # samplename_L001_R1_001.fastq.gz;
    # samplename_L001_R2_001.fastq.gz
# output:
    # samplename_L001_R1_001_fastqc.htlm;
    # samplename_L001_R1_001_fastqc.zip;
    # samplename_L001_R2_001_fastqc.htlm;
    # samplename_L001_R2_001_fastqc.zip;

###################################################
##### Package/program Requirement
##################################################

# fastQC : FastQC A quality control application for high throughput sequence data
# advice for installation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# if conda is installed: "conda install -c bioconda fastqc"


###################################################
##### Run FastQC on each fastq file on the command line :
##################################################

fastqc *_L001_R1_001.fastq* *_L001_R2_001.fastq*
mkdir 1_QC_before_trim
mv *_fastqc.* 1_QC_before_trim

## view htlm file with safari:
# "open -a "Safari" PathToFile/samplename_L001_R1_001_fastqc.html"


# system time and date at end
date
