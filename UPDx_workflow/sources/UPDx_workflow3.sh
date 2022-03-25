#!/bin/bash

##Software required
## unix/linux terminal
## Cutadapt
## FastQC
## R version 3.6.3 and corresponding packages
## Blast + program

# UPDx WORKFLOW OVERVIEW 
## FASTQC before trimming
## TRIMMING with CUTADAPT
## FASTQC after trimming
## Dada2 analysis
    # Quality Trimming
    # Error model
    # Dereplication
    # Inferring ASVs/Denoising
    # Merging
    # Generationg a count table / OTU Matrice
    # Chimera identification
    # Overview of counts throughout
    # Assigning taxonomy
## Blastn Sequences identification
## Phylogenetic analysis of the ASV sequences using custom and local database
## Report analysis generation


# system time and date at start
date

###################################################
##### SETTINGS Paths
##################################################

## Define path to main folder "UPDx_workflow" from code location
Path_to_WorkflowFolder="../"
cd ${Path_to_WorkflowFolder}
Path_ini=$(pwd)

## Path to scripts 
Path_script="${Path_ini}/sources"

## Path to input files (raw reads file R1 and R2)
Path_data="${Path_ini}/data"

## Path to custom parasitic database for phylogeny analysis
Path_database="${Path_ini}/sources/r/UPDx_references_tree.fasta"

###################################################
##### FASTQC BEFORE TRIMMING  Bash script
##################################################
cd $Path_data

chmod 0755 $Path_script/sh/UPDx_fastQC_before_trim.sh
echo "FASTQC BEFORE TRIMMING STEP"
sh $Path_script/sh/UPDx_fastQC_before_trim.sh
echo "FASTQC BEFORE TRIMMING STEP DONE


"
###################################################
##### CUTADAPT TRIMMING  Bash script
##################################################
cd $Path_data

chmod 0755 $Path_script/sh/UPDx_cutadapt.sh
echo "CUTADAPT TRIMMING STEP"
sh $Path_script/sh/UPDx_cutadapt.sh
echo "CUTADAPT TRIMMING STEP DONE


"
###################################################
##### FASTQC AFTER TRIMMING  Bash script
##################################################
cd $Path_data/2_Cutadapt_Output

chmod 0755 $Path_script/sh/UPDx_fastQC_after_trim.sh
echo "FASTQC AFTER TRIMMING STEP"
sh $Path_script/sh/UPDx_fastQC_after_trim.sh
echo "FASTQC AFTER TRIMMING STEP DONE


"

###################################################
##### Sequence analysis (dada2, blastn, threshold)
##################################################
cd $Path_data/2_Cutadapt_Output

chmod 0755 $Path_script/r/UPDx_analysis_taxonomy.R
echo "SEQUENCE ANALYSIS STEP"
#Rscript $Path_script/r/UPDx_analysis_taxonomy.R  $path_blastn --save
Rscript $Path_script/r/UPDx_analysis_taxonomy.R 
echo "SEQUENCE ANALYSIS STEP DONE


"
###################################################
##### Sequence phylogeny 
##################################################
cd $Path_data/5_Blastn_Output

chmod 0755 $Path_script/r/UPDx_analysis_phylogeny.R
echo "SEQUENCE PHYLOGENY STEP"
Rscript $Path_script/r/UPDx_analysis_phylogeny.R $Path_database --save
echo "SEQUENCE PHYLOGENY STEP DONE


"

###################################################
##### Report pdf
##################################################
cd $Path_data/6_Phylo_Output

chmod 0755 $Path_script/r/UPDx_analysis_report2.Rmd
chmod 0755 $Path_script/r/UPDx_analysis_report1.R
echo "REPORT STEP"
Rscript $Path_script/r/UPDx_analysis_report1.R $Path_script --save
echo "REPORT STEP DONE


"

# system time and date at end
date
