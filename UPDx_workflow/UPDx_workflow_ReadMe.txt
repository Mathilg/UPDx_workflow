######## READ ME "UPDx sequences analysis"
 
### UPDx PROJECT DESCRIPTION = 
# UPDx: Universal Parasite Detection assay. 
# Targeted amplicon deep sequencing (TADS) of the 18S rRNA gene for the detection and characterization of parasitic communities in biological samples. 
# UPDx protocol employs a single pair of universal primers to capture all blood-borne parasites while reducing host 18S rRNA template using restriction enzymes to digest the 18S rRNA gene at cut sites present only in the host sequence prior to PCR amplification, thus enhancing the amplification of parasite 18S rRNA for TADS. 
# Amplicon Sequencing is achieved with MiSeq Illumina, 251bp Paired end read, amplicon length expected around 180 bp. 

### UPDx WORKFLOW FOLDER CONTENT:
# UPDx_workflow_ReadMe.txt: readme file explaining workflow setup et organization.
# “data” folder: input files location (raw read files produced by paired-end sequencing)
# “data-example” folder: dataset of 6 input files corresponding to the UPDx sequencing of two parasite-free (negative controls) and or four parasite-infected blood samples (including Leishmania infantum ‘li’, Plasmodium falciparum ‘Pf’, Plasmodium vivax ‘Pv’, Babesia microti ‘Bmic’. These input files can be used as workflow implementation control. 
# “sources” folder: workflow scripts location
		# 3 Main workflow scripts: you can choose between three different versions of the workflow. 
			# script 1 "UPDx_workflow1.sh": sequences analysis + taxonomic identification (Workflow steps 1 to 5)
			# script 2 "UPDx_workflow2.sh": sequences analysis + taxonomic identification + phylogeny analysis (Workflow steps 1 to 6)
			# script 3 "UPDx_workflow3.sh": sequences analysis + taxonomic identification + phylogeny analysis + pdf report (Workflow steps 1 to 7)
		# Folder "sh": Bash scripts
		# Folder "r": R scripts
		# Folder "r_package_install": include the R script allowing the installation of the R packages required.

### UPDx WORKFLOW STEPS OVERVIEW
# This workflow is adapted from the "full example workflow for amplicon data" published by Michael D. Lee, Ph.D and available at https://astrobiomike.github.io/amplicon/
	
	# Workflow INPUT: Fastq files (R1 and R2 illumina fastq files)
	# Step 1: Sequence quality checking with FastQC before primers trimming with Cutadapt.
	# Step 2: Primers removal (at 5' and 3' sequences extremities) and length sequences filtering (length between 145-200 bp) using Cutadapt
	# Step 3: Quality checking with FastQC after primers trimming with Cutadapt.
	# Step 4: Sequences analysis using dada2 R package (available at https://benjjneb.github.io/dada2/)
		# Quality trimming/filtering
			# Visualisation of reads Quality before and after filtering (plots)
			# Filtering: Removal of PhiX sequences, filtering based on Quality Scores (QS) and sequence length. 
				# Truncating reads at the first instance of quality score less than 2
				# After truncating, sequences with length < 145 b or  > 200 b are discarded.
				# After truncating, Reads with overall QS < 15 are discarded
				# After truncating, Reads with more than 2 erroneous base calls are discarded
			# Nb: More on Quality Score on https://emea.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html
		# Error model generation of the dataset
		# Dereplication: Sequences sharing 100% of sequence identity are pooled as one sequences. It also generates a new quality-score profile of each unique sequence based on the average quality scores of each base of all of the sequences that were replicates of it.
		# Inferring true biological sequences ASVs (Amplicon Sequence Variants (ASVs)). 
			# ASV definition by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence is more likely to be of biological origin or more likely to be spurious (using error model).
		# Merging of the forward and reverse reads: Minimum of 150 bp overlap, overhang discarded.
		# Chimera identification: Chimeras are sequences formed from two or more biological sequences joined together. 
			# DADA2 identifies likely chimeras by aligning each sequence with those that were recovered in greater abundance and then seeing if there are any lower-abundance sequences that can be made exactly by mixing left and right portions of two of the more-abundant ones. These are then removed.
		# Outputs:
			# fasta files of the ASVs sequences: "AsvSequences.fa"
			# Summary of the trimming and filtering steps: "SummaryTrimming.csv" 
			# Summary of the ASV sequence metrics: "AsvSequencesMetrics.csv" 
	# Step 5: Taxonomic assignement of the ASVs using NCBI Blastn algorithm 
		# Taxonomic assignment of the ASVs using remote Blastn program 
		# The cutoff value for positivity are added to each ASV. 
 			# cutoff min (20reads, see Flaherty et al., 2018)
			# cutoff max or "Sliding cutoff" : run-dependant, based on negative control contamination 
				# In each run, at least 3 Blood negative samples: sliding Cutoff = [ ( ucontamin_run ) + 4 ( S.D. ) ] * Sreads = Cutoffmax
  					# with ucontamin_run : mean proportion of contaminating reads in negative blood samples in the particular run.
  					#S.D. : standard deviations
  					#Sreads : the number of reads sequenced for the sample.
			# If cutoff max > cutoff min, cutoof max will be applied, if cutoff max < cutoff min, cutoff min will be applied
		# ASVs assigned to “parasitic” taxa are extracted as fasta files for phylogenetic analysis
			# parasitic sequences: selection of ASVs not identified as "Chordata","Fungus", "Fungi", "synthetic", "Bacteria", "Virus", "Plantae", "Arthropoda" (that are not Ixodida*), "Viridiplantae", "uncultured_eukaryote"
			# *some parasitic sequence, such as Babesia spp. sequences, are currently badly assigned to Ixodes sequences in NCBI database so we are excluding only arthropoda sequences that are not Ixodida (this is not pertinent if you are analyzing ticks samples!)
		#Outputs: 
			# Raw outputs from Blastn: "Blastn_RawOutput.csv" 
			# Summary of the all ASV sequence metrics and taxonomy: "AllAsvSequenceTaxonomyMetrics.csv" 
			# Summary of the parasitic ASV sequence metrics and taxonomy: "ParasiticAsvSequenceTaxonomyMetrics.csv" 
			# Fasta file with all parasitic ASVs sequences considered as "positive": "AllParasiticAsvSequences.fa"
			# For each sample, fasta file including parasitic ASVs sequences considered as "positive": "Samples_ParasiticAsvSequences.fa"
	# Step 6: Phylogenetic analysis pf the parasitic ASVs sequences using a custom and local database
		#Outputs: 
			# PDF: Global trees and subtrees for each sample and for the whole parasitic ASVs sequences identified
	# Step 7: Generation of an analysis report (PDF)

### UPDx WORKFLOW SETTINGS: Tools required (version highly recommended)
# Ubuntu 18.04 LTS (Bionic Beaver)
# Cutadapt v1.15 (available at https://cutadapt.readthedocs.io/en/stable/installation.html)
# FastQC v0.11.5 (available at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
# BLAST+ standalone: ncbi-blast-2.13.0+-x64-linux (available at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
	#Blastn program has to be on the $PATH
# R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# dada2 R package (https://benjjneb.github.io/dada2/dada-installation.html)
# Other R package: can be installed using R scripts available for packages installation "sources/r_packages_install" folder
	# R scripts available for packages installation  "sources/r_packages_install" folder.
	#  Workflow1 : requires "Rcpp", "readr", "tidyr", "tidyverse", "rlang", "taxize",  "filesstrings", "seqinr", "NCmisc"
		# can be installed using "R_packages_workflow1.R" script
	#  Workflow2 : requires above packages + "ggplot2", "RColorBrewer", "reshape", "dbscan", "gridExtra", "stringr", "cluster", "phylogram",  "ggdendro", "colorspace", "ape", "Biostrings", "dplyr", "msa", "DECIPHER", "treeio", "ggtree"
		# can be installed using "R_packages_workflow1.R" and "R_packages_workflow2.R" scripts
	#  Workflow3 : requires above packages + "tinytex", "knitr", "stringr", "dplyr", "readr", "tibble", "kableExtra"
		# can be installed using "R_packages_workflow1.R", "R_packages_workflow2.R", "R_packages_workflow3.R"and scripts
# if you choose to run workflow3 you will also need: Pandoc version 1.12.3 
# if you choose to run workflow3 you will also need: libfontconfig1-dev

### UPDx WORKFLOW EXECUTION
# Install the required tools (according to the version of the workflow you choose to run)
# Download the whole “UPDx_workflow” folder.
# Data to analyze:
	# Place your fastq files (inputs) in the "data" folder
		# Data naming requirement: please follow this rule of naming:
			# no underscore in {sample} in {sample}_L001_R1_001.fastq or {sample}_L001_R2_001.fastq 
			# negative control name hast to containe one of these expressions:  "Negative", "NEG", "neg", "uRBC"
# UPDx_workflow execution:
	# 3 Main workflow scripts: you can choose between three different versions of the workflow. 
			# script 1 "UPDx_workflow1.sh": sequences analysis + taxonomic identification (Workflow steps 1 to 5)
			# script 2 "UPDx_workflow2.sh": sequences analysis + taxonomic identification + phylogeny analysis (Workflow steps 1 to 6)
			# script 3 "UPDx_workflow3.sh": sequences analysis + taxonomic identification + phylogeny analysis + pdf report (Workflow steps 1 to 7)
	# Execute the chosen script: 
		# example: workflow1:
		cd {LocationOfTheFolderOnYourComputer}/UPDx_workflow/sources
		sudo chmod 0777 UPDx_workflow1.sh
		./UPDx_workflow1.sh


