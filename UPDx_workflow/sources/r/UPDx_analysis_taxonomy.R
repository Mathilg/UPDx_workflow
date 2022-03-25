#!/usr/bin/Rscript

########################################################################
####UPDx READS ANALYSIS 
########################################################################
############ WF OVERVIEW

############ dada2 analysis
### Quality trimming/filtering
### Generating error model of the data
### De-replication
### INFERRING ASVs (Amplicon Sequence Variants): infer true biological sequences
### Chimera identification
### Adding minimum cutoff for UPDx analysis (see "Readme.txt" file)
### OUTPUTS : 
#### fasta files of the ASVs sequences: "AsvSequences.fa"
#### Summary of the trimming and filtering steps: "SummaryTrimming.csv" 
#### Summary of the ASV sequence metrics: "AsvSequencesMetrics.csv" 

############ Taxonomic assignments BLASTn
### Appending taxonomy (full and short taxonomy info)
### Merging taxonomy info and metrics of ASVs
### Adding max cutoff for UPDx analysis (see "Readme.txt" file)
### OUTPUTS
#### Raw outputs from Blastn: "Blastn_RawOutput.csv" 
#### Summary of the all ASV sequence metrics and taxonomy: "AllAsvSequenceTaxonomyMetrics.csv" 
#### Summary of the parasitic ASV sequence metrics and taxonomy: "ParasiticAsvSequenceTaxonomyMetrics.csv" 
#### Fasta file with all parasitic ASVs sequences considered as "positive": "AllParasiticAsvSequences.fa"
#### For each sample, fasta file including parasitic ASVs sequences considered as "positive": "Samples_ParasiticAsvSequences.fa"

############ Package/program Requirement
# Package installed with the specific script provided with the WF
library(dada2)
library(filesstrings)
library(readr)
library(taxize)
library(dplyr)
library(tidyr)
library(purrr)
library(seqinr)
library(NCmisc)


########################################################################
####UPDx READS ANALYSIS with Dada2
########################################################################

#####################################
############ SAMPLE NAME
# List of sample names by scanning the "samples.txt" file produced during cutadapt step
samples <- scan("../samples.txt", what="character")
samples

#####################################
############ Quality trimming/filtering
# INPUT variables: file names for the forward and reverse trimmed reads
forward_reads <- paste0(samples, "_L001_R1_001_trim1.fastq")
forward_reads 
reverse_reads <- paste0(samples, "_L001_R2_001_trim1.fastq")
reverse_reads 
# OUTPUT variables: file names for the forward and reverse filtered reads
filtered_forward_reads <- paste0(samples, "_L001_R1_001_filtered.fastq")
filtered_forward_reads 
filtered_reverse_reads <- paste0(samples, "_L001_R2_001_filtered.fastq")
filtered_reverse_reads 
# QUALITY CHECKING before trimming (bases in x-axis and quality score on the y-axis)
plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)
# FILTERING
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,                 
                              rev=reverse_reads, filt.rev=filtered_reverse_reads, 
                              rm.phix=c(TRUE,TRUE),       #remove PhiX sequences from R1 and R2
                              truncQ = c(2,2),            #Truncate reads at the first instance of a quality score less than or equal to truncQ in R1 and R2. Quality score of 2 --> The lowest value usually found in practice is Q=2 (P=0.63), which means the base call is more likely to be wrong than correct.
                              minLen= c(145,145),         #min read length for R1 and R2 (145bp)
                              maxLen= c(200,200),         #max read length for R1 and R2 (200bp)
                              maxEE=c(2,2),               #After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
                              minQ=c(15,15))              #reads with QS inf to 15 will be discarded
filtered_out                           
# Removal of empty file (all reads filtered and discarded)
# List of the filtered reads files
filtered_forward_reads2 <- list.files(pattern = "*_R1_001_filtered*") 
filtered_reverse_reads2 <- list.files(pattern = "*_R2_001_filtered*") 
# files removed by filtering
A <- unlist(filtered_forward_reads)
B <- unlist(filtered_forward_reads2)
empty_filtered_forward_reads <- setdiff(A,B)
C <- unlist(filtered_reverse_reads)
D <- unlist(filtered_reverse_reads2)
empty_filtered_reverse_reads <- setdiff(C,D)
# QUALITY CHECKING after trimming
plotQualityProfile(filtered_forward_reads2)
plotQualityProfile(filtered_reverse_reads2)

#####################################
############ Generating error model of the data
# Learning the specific error-signature of our dataset
# Error rates are learned by alternating between sample inference and error rate estimation until convergence.
err_forward_reads <- learnErrors(filtered_forward_reads2, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads2, multithread=TRUE)
# plots
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)

#####################################
############ De-replication
# Refresh Samples list
samples1 <- list.files(all.files = TRUE, pattern = "_filtered")
newName <- sub("_L001_R1_001_filtered.fastq", "", samples1)
newName2 <- sub("_L001_R2_001_filtered.fastq", "", newName)
samples1 <- unique(newName2)
# DEREPLICATION : pool sequences sharing 100% of sequence identity 
# When DADA2 de-replicates sequences, it also generates a new quality-score profile of each unique sequence based on the average quality scores of each base of all of the sequences that were replicates of it.
derep_forward <- derepFastq(filtered_forward_reads2, verbose=TRUE)
derep_forward
names(derep_forward) <- samples1 
derep_reverse <- derepFastq(filtered_reverse_reads2, verbose=TRUE)
names(derep_reverse) <- samples1

#####################################
############ INFERRING ASVs (Amplicon Sequence Variants): infer true biological sequences. 
# incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence is more likely to be of biological origin or more likely to be spurious.
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE)
dada_forward
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)
dada_reverse

#####################################
############ READS MERGING
# Merging conditions: 
# Forward and reverse reads overlap by at least 150 bases, 100% identity in the overlap region.
# trimOverhang option to TRUE : if reads are longer than overlap, region is trimmed 
merged_amplicons <- mergePairs(dada_forward, derep_forward, 
                               dada_reverse, derep_reverse, 
                               trimOverhang=TRUE, minOverlap=150)
merged_amplicons

#####################################
############ Chimera identification
# chimera : Chimeras are sequences formed from two or more biological sequences joined together (incomplete extension steps).Chimeras are common in amplicon sequencing when closely related sequences are amplified. 
# DADA2 identifies likely chimeras by aligning each sequence with those that were recovered in greater abundance and then seeing if there are any lower-abundance sequences that can be made exactly by mixing left and right portions of two of the more-abundant ones. These are then removed.
seqtab <- makeSequenceTable(merged_amplicons)
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) 
# check abundance: 
sum(seqtab.nochim)/sum(seqtab) 

#####################################
############ OUTPUTS 
# dada2 output folder
dir.create("../4_Dada2_Output")
setwd("../4_Dada2_Output")
#Export the table of the Trimming and filtering metrics (number of reads processed among the WF): "SummaryTrimming.csv"
getN <- function(x) sum(getUniques(x))
summary_tab1 <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                           filtered=filtered_out[,2])
summary_tab2 <- data.frame(row.names=samples1, dada_f=sapply(dada_forward, getN),
                           dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                           nonchim=rowSums(seqtab.nochim))
summary_tab3 <- merge(summary_tab1, summary_tab2, by=0, all=TRUE)
summary_tab3[is.na(summary_tab3)] <- 0                 # replace NA values
colnames(summary_tab3)[1] <- "Samples"
write.csv(summary_tab3, "SummaryTrimming.csv", sep= " ", quote = F, col.names = FALSE)
# Export the fasta file that contains all the ASVs sequences "AsvSequences.fa" 
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
asv_fasta
write(asv_fasta, "AsvSequences.fa")
# Export the table presenting ASVs abundance among samples: "AsvSequencesMetrics.csv
# Adding of the minimum Cutoff: ASV sequence has to be detected  at least 20 times within a sample to be considered as positive. 
asv_csv <- cbind(asv_headers, asv_seqs)   # table of ASVs name and sequences
asv_tab <- as.data.frame(as.table(seqtab.nochim)) # table of ASVs sequence and the number of reads corresponding to the sequence for each sample
colnames(asv_tab) <- c("Sample", "Sequence",  "Reads_number") 
asv_tab$Sample <- gsub("_L001_R1_001_filtered.fastq", "", asv_tab$Sample)
asv_tabfinal <- merge(asv_csv, asv_tab, by.x = "asv_seqs", by.y = "Sequence") # merging of the table
asv_tabfinal <- asv_tabfinal[,c(2,1,3,4)]
asv_tabfinal$Tot_reads <- NA
summary_tab4 <- subset(summary_tab3, summary_tab3$filtered!=0)
asv_tabfinal$Tot_reads <- summary_tab4$nonchim  # addition of total read number info
asv_tabfinal$Min_cutoff <- 20 # addition of the "Min_Cutoff" threshold (20 reads for UPDx analysis, see papers)
FinTab <- asv_tabfinal[with(asv_tabfinal, order(asv_tabfinal$Sample, asv_tabfinal$Reads_number, decreasing=TRUE)),]
write.table(FinTab, "AsvSequencesMetrics.csv", sep=" ", quote=F, row.names = FALSE)

########################################################################
#### ASVs taxonomic assignment 
########################################################################

# Blastn rule: maximum of 100 sequences per query file
asvseq <- read.fasta(file ="AsvSequences.fa")
seq_number <- length(asvseq)

if (seq_number > 100){
  file.split("AsvSequences.fa", size = 140, same.dir = TRUE, verbose = TRUE,
             suf = "part", win = TRUE)
  sub_files <- list.files(path=getwd(), pattern= ".fa_part" )
} else {
  sub_files <- "AsvSequences.fa"
}

# BLASTn execution (Blastn supposed to be on the $PATH)
for (nfile in sub_files) {
  output_name <- paste0("Blastn_RawOutput",nfile,".csv")
  blast_output <- system2(command = "blastn",                    #location of the BLastn bin on your computer
                          args = c("-db", "nt", "-remote",    # -db:database, nt:nucleotide database, -remote: NCBI server use
                            "-query", nfile,           #-query + path to input File
                            "-outfmt", 6,              # output file as csv tabular #6, standard fields --> qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                            "-evalue", "1e-6",         # minimum evalue
                            "-max_target_seqs", "1",   # maximum of subject
                            "-max_hsps", "1",          # maximum of alignment per subject
                            #"-ungapped",              # Perform ungapped alignment, not activated
                            "-out", output_name))      # Name of the output file
}

# Create Blast Output folder
file_csv <- list.files(pattern="*.fa.csv")
dir.create("../5_Blastn_Output")
file.move(file_csv,  "../5_Blastn_Output")
setwd("../5_Blastn_Output")

# Merging the different blast output tables
blast_outputR2 <- data.frame()
for (csvfile in file_csv) {
  dd <- read_table(csvfile, col_names = FALSE)
  blast_outputR2 <- rbind(blast_outputR2, dd)
}
blast_outputR2

# Export Blastn Output table
colnames(blast_outputR2) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                               "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                               "S.start", "S.end", "E", "Bits" )  
write.table(blast_outputR2, "Blastn_RawOutput.csv", sep="\t", quote=F, row.names=FALSE)

# Remove the single file
for (csvfile in file_csv) {
  file.remove(csvfile)
}

# Appending taxonomy : Complete LINEAGE "FULL TAXO"
# From the "SubjectID" provided by GenBank in blast results, we can get the sequence UID using the function genbank2uid() of Taxize package
seq_IDs <- blast_outputR2
seq_IDs_nest <- seq_IDs %>%
  nest(., SubjectID) %>%       
  mutate(., uids = map(data, possibly(genbank2uid, otherwise = NA, quiet = TRUE))) 
seq_IDs2 <- seq_IDs_nest %>% unnest(c(data, uids)) 
# From the "UIDs" we can get the full classification using the function classification() of Taxize package
seq_IDs3 <- seq_IDs2
seq_IDs3$uids <- as.numeric(seq_IDs3$uids)
n=1
for (i in seq_IDs3$uids) {        
  clasi <- classification(i, db = 'ncbi')   
  seq_IDs3[n,14] <- toString(cbind(clasi))   
  n=n+1
  Sys.sleep(0.4)  # timeout between queries required to avoid "too many request error)
}
# rename colonne 14
colnames(seq_IDs3)[14] <- "Full_taxonomy"
seq_IDs4 <- seq_IDs3
# Format the writing of the taxonomy
seq_IDs4$Full_taxonomy[grep(", ",seq_IDs4$Full_taxonomy)] <- gsub(", ","~",seq_IDs4$Full_taxonomy)
seq_IDs4$Full_taxonomy[grep(" ",seq_IDs4$Full_taxonomy)] <- gsub(" ","_",seq_IDs4$Full_taxonomy)
Full_taxo_table <- seq_IDs4
Full_taxo_table$QueryID[grep("ASV",Full_taxo_table$QueryID)] <- gsub("ASV",">ASV",Full_taxo_table$QueryID[grep("ASV",Full_taxo_table$QueryID)])

print("line 267")

# Appending taxonomy : SPECIES AND GENUS "SHORT TAXO"
# From the "SubjectID" provided by GenBank in blast results, we can get the sequence UID using the function genbank2uid() of Taxize package
seq_IDss <- blast_outputR2
seq_IDss_nest <- seq_IDss %>%
  nest(., SubjectID) %>%        
  mutate(., uids = map(data, possibly(genbank2uid, otherwise = NA, quiet = TRUE))) 
seq_IDss2 <- seq_IDss_nest %>% unnest(c(data, uids)) 
# From the "UIDs" we can get the full classification using the function classification() of Taxize package
seq_IDss2$uids <- as.numeric(seq_IDss2$uids)
seq_IDss_nest2 <- seq_IDss2 %>%
  nest(., uids) %>%
  mutate(., tax_info = map(data, possibly(ncbi_get_taxon_summary, otherwise = NA, quiet = TRUE)))
seq_IDss3 <- seq_IDss_nest2 %>% unnest(c(data, tax_info))
# Format the writing of the taxonomy
seq_IDss4 <- seq_IDss3
seq_IDss4$name[grep(" ",seq_IDss4$name)] <- gsub(" ","_",seq_IDss4$name)
Short_taxo_table <- seq_IDss4
Short_taxo_table$QueryID[grep("ASV",Short_taxo_table$QueryID)] <- gsub("ASV",">ASV",Short_taxo_table$QueryID[grep("ASV",Short_taxo_table$QueryID)])

print("line 288")
 
################################
############ Merging taxonomy info and ASVs metrics "AsvSequencesMetrics.csv"
Count_table <- read_delim("../4_Dada2_Output/AsvSequencesMetrics.csv", 
                          " ", escape_double = FALSE, trim_ws = TRUE)
## FULL TAXO
Final_full_taxo_table <- merge(Count_table,Full_taxo_table, by.x = "asv_headers", by.y =  "QueryID" , all = TRUE)
## SHORT TAXO
Final_short_taxo_table <- merge(Count_table,Short_taxo_table, by.x = "asv_headers", by.y =  "QueryID" , all = TRUE)
## Merging "full taxonomy" and "short taxonomy" tables
Final_taxo_table <- cbind(Final_full_taxo_table,Final_short_taxo_table[,c(20,21)])
colnames(Final_taxo_table)[19] <- "Full_taxonomy"

################################
############ Adding max cutoff for UPDx analysis (see "Readme.txt" file)
# Since the beginning, a minimal cutoff of 20 reads has been applied. Meaning that at least 20 reads of an ASVs are required to say this sequence has been detected.
# A sliding cutoff or "cutoff max" is also added and calculated according to negative control contamination which is run-dependant
# Data from negative control 
Bloodneg <- Final_taxo_table[grepl("Negative|NEG|neg|uRBC", Final_taxo_table[["Sample"]]),]
# Number of "parasitic" sequences identified in negative control (contamination)
# Parasitic sequences: selection of ASVs not identified as "Chordata","Fungus", "Fungi", "synthetic", "Bacteria", "Virus", "Plantae", "Arthropoda", "Viridiplantae", "uncultured_eukaryote"
Bloodneg_reads <- subset(Bloodneg, Bloodneg$Reads_number != 0)
Bloodneg_parasite <- Bloodneg_reads %>% filter(!grepl("Chordata", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Fungus", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Fungi", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("synthetic", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Bacteria", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Virus", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Plantae", Full_taxonomy))
#Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Arthropoda", Full_taxonomy))  #see above the specific case of the Ixodida order. If you analyse tick samples, you should unlock this line and put lines 320 to 323 in comments 
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Viridiplantae", Full_taxonomy))
Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("uncultured_eukaryote", Full_taxonomy))
# We want to remove Arthropoda sequences but some of Babesia spp. sequences in NCBI database are badly identified as Ixodes sequence. 
# We can remove arthropod sequences except the Ixodida ones as we are working on blood samples.
dd <- Bloodneg_parasite[grepl("Arthropoda[d+]Ixodida", Bloodneg_parasite$Full_taxonomy), ]  
if (dd <- FALSE) {
  Bloodneg_parasite <- Bloodneg_parasite %>% filter(!grepl("Arthropoda", Full_taxonomy))
}
neg_sample <- unique(Bloodneg_parasite$Sample)
cont_reads <- c()
for (i in neg_sample) {
  subseti <- subset(Bloodneg_parasite, Bloodneg_parasite$Sample == i)
  sum_readi<- sum(subseti$Reads_number)
  cont_readsi <- (sum_readi/unique(subseti$Tot_reads))
  cont_reads <- append(cont_reads, cont_readsi)
}
cont_reads
mean_prop_contamin_run <- mean(cont_reads) 
sd_prop_contamin_run <- sd(cont_reads)
# Maximum cutoff determination
Final_taxo_table$mean_prop_contamin_run <- mean_prop_contamin_run
Final_taxo_table$sd_prop_contamin_run <- sd_prop_contamin_run
Number_negative_sample <- length(unique(Bloodneg_parasite$Sample))
Final_taxo_table$Max_cutoff<- round(((mean_prop_contamin_run + (Number_negative_sample * sd_prop_contamin_run)) * Final_taxo_table$Tot_reads))
Final_taxo_table <- Final_taxo_table[,c(1:6,24,7:21)]

################################
############ Cutoff for positivity (see "Readme.txt" file)
# If cutoff max > cutoff min, cutoof max will be applied, if cutoff max < cutoff min, cutoff min will be applied
Final_taxo_table$Min_cutoff <- 20
Final_taxo_table$Results <- "NEGATIVE"
Final_taxo_table$Cutoff <- NA
Final_taxo_table <- Final_taxo_table[,c(1:7,24,23,8:22)]
for (i in 1:nrow(Final_taxo_table)){
  if (is.na(Final_taxo_table$Max_cutoff[i])){
    Final_taxo_table$Results[i][Final_taxo_table$Reads_number[i] > 20] <- "POSITIVE"
    Final_taxo_table$Cutoff[i] <- "Min"
  } else {
    if (Final_taxo_table$Min_cutoff[i] > Final_taxo_table$Max_cutoff[i]) {
      Final_taxo_table$Results[i][Final_taxo_table$Reads_number[i] > 20] <- "POSITIVE"
      Final_taxo_table$Cutoff[i] <- "Min"
    } else {
      Final_taxo_table$Results[i][Final_taxo_table$Reads_number[i] > Final_taxo_table$Max_cutoff[i]] <- "POSITIVE"
      Final_taxo_table$Cutoff[i] <- "Max"
    }
  }
}
Final_taxo_tablef <- subset(Final_taxo_table, Final_taxo_table$Reads_number > 0)
Final_taxo_tablef <-Final_taxo_tablef[order(Final_taxo_tablef$Sample, Final_taxo_tablef$Reads_number, decreasing = TRUE),]

################################
############ OUTPUTS

######Tables: Selection of parasitic ASV only
# Export "AllAsvSequenceTaxonomyMetrics.csv" table: ASVs metrics + taxonomy + cutoff and "positivity" results
write.table(Final_taxo_tablef, "AllAsvSequenceTaxonomyMetrics.csv", sep="\t", quote=F, row.names=FALSE)
# Export "ParasiticAsvSequenceTaxonomyMetrics.csv" table: ASVs metrics + taxonomy + cutoff fo only sequences identify as Parasite and positive
Pos_seq <- subset(Final_taxo_table, Final_taxo_table$Results == "POSITIVE")
Pos_seq2 <- Pos_seq %>% filter(!grepl("Chordata", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Fungus", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Fungi", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("synthetic", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Bacteria", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Virus", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Plantae", Full_taxonomy))
#Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Arthropoda", Full_taxonomy))  #see above the specific case of the Ixodida order. If you analyse tick samples, you should unlock this line and put lines 320 to 323 in comments 
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Viridiplantae", Full_taxonomy))
Pos_seq2 <- Pos_seq2 %>% filter(!grepl("uncultured_eukaryote", Full_taxonomy))
# We want to remove Arthropoda sequences but some of Babesia spp. sequences in NCBI database are badly identified as Ixodes sequence. 
# We can remove arthropod sequences except the Ixodida ones as we are working on blood samples.
ee <- Pos_seq2[grepl("Arthropoda[d+]Ixodida", Pos_seq2$Full_taxonomy), ]  
if (ee <- FALSE) {
  Pos_seq2 <- Pos_seq2 %>% filter(!grepl("Arthropoda", Full_taxonomy))
}
write.table(Pos_seq2, "ParasiticAsvSequenceTaxonomyMetrics.csv", sep="\t", quote=F, row.names=FALSE)

######Fasta files with parasitic sequences
# global fasta file export & export of Fasta files for each sample
Pos_seq3 <- Pos_seq2
Pos_seq3$Seq_name <- paste(Pos_seq3$asv_headers,Pos_seq3$Sample, Pos_seq3$Perc.Ident,Pos_seq3$name, sep="_")
Seq_final <- Pos_seq3[,c(25,2)]
Seq_final <- c(rbind(Pos_seq3$Seq_name, Pos_seq3$asv_seqs))
write(Seq_final, "AllParasiticAsvSequences.fa")
Pos_seq4 <- Pos_seq2
samples <- scan("../samples.txt", what="character")
samples
for (i in samples) {
  sample_Seq_finali <- subset(Pos_seq2, Pos_seq2$Sample == i )
  sample_Seq_finali$Seq_name <- paste(sample_Seq_finali$asv_headers,sample_Seq_finali$Sample, sample_Seq_finali$Perc.Ident,sample_Seq_finali$name, sep="_")
  sample_Seq_finali <- sample_Seq_finali[,c(25,2)]
  Seq_final <- c(rbind(sample_Seq_finali$Seq_name, sample_Seq_finali$asv_seqs))
  write(Seq_final, paste0(i,"_ParasiticAsvSequences.fa"))
}