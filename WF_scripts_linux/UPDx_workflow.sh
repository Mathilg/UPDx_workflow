#!/bin/bash

## Software required
## Cutadapt
## Fastqc

# UPDx WORKFLOW -- Linux part

## OVERVIEW WF

## FASTQC before trimming
## PRIMERS TRIMMING with CUTADAPT
## FASTQC after trimming


# system time and date at start
date

## Folder sontaining the scripts
Path_script="/home/mathilde/UPDx_workflow_Linux/WF_scripts_linux"

## Initial working directory
## folder containing raw reads file R1 and R2
Path_ini="/home/mathilde/A_nUPDx_EXP2/N"


###################################################
##### FASTQC BEFORE TRIMMING  Bash script
##################################################
cd $Path_ini

chmod 0755 $Path_script/UPDx_fastQC_before_trim.sh
echo "FASTQC BEFORE TRIMMING STEP"
sh $Path_script/UPDx_fastQC_before_trim.sh
echo "FASTQC BEFORE TRIMMING STEP DONE


"


###################################################
##### CUTADAPT TRIMMING  Bash script
##################################################
cd $Path_ini

chmod 0755 $Path_script/UPDx_cutadapt.sh
echo "CUTADAPT TRIMMING STEP"
sh $Path_script/UPDx_cutadapt.sh
echo "CUTADAPT TRIMMING STEP DONE


"

###################################################
##### FASTQC AFTER TRIMMING  Bash script
##################################################
cd $Path_ini/2_Cutadapt_output

chmod 0755 $Path_script/UPDx_fastQC_after_trim.sh
echo "FASTQC AFTER TRIMMING STEP"
sh $Path_script/UPDx_fastQC_after_trim.sh
echo "FASTQC AFTER TRIMMING STEP DONE


"

# system time and date at end
date
