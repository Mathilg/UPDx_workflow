#!/usr/bin/Rscript

########################################################################
########################################################################
## UPDx PHYLOGENETIC TREE (>1 file to analyze)

###################################################
##### Package/program Requirement
##################################################

##Packages required
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(dbscan)
library(gridExtra)  
library(stringr)
library(cluster)
library(phylogram) 
library(ggdendro)
library(msa)
library(ggtree)
library(colorspace)
library(ape)
library(seqinr)
library(ggtree)
library(Biostrings)
library(DECIPHER)
library(treeio)
library(taxize)
library(dplyr)
library(filesstrings)

#############################
####### IMPORTING DATA + Merging files + Phylogenetic tree
#############################

## import the reference database (fasta file) 
## intercept the parameters/arguments passed by shell
args <- commandArgs()
print(args)
UPDx_DBref <- args[6]
DBref <- readDNAStringSet(UPDx_DBref)

## First removed empty fasta file (for example Blood negative samples)
FastaFiles0 <- list.files(pattern= ".fa" )
# Use file.size() immediate, instead of file.info(docs)$size:
inds <- file.size(FastaFiles0) <= 2 
# Remove all documents with file.size = 1 from the directory
file.remove(FastaFiles0[inds])
##list the fasta files with parasitic sequences
FastaFiles <- list.files(pattern= ".fa" )
FastaFiles


############# Phylogenetic Tree == apply to each Fasta file (one fasta file per sample)

######### Within the same loop
##### Global tree = all the ASVs sequences found in the sample and Reference sequences
for (file_name in FastaFiles) {
  #sample <- readDNAStringSet(paste0(WD_path,"/",file_name))
  sample <- readDNAStringSet(file_name)
  mergedfile <- c(DBref, sample)
  mergedfile <- OrientNucleotides(mergedfile, reference = which.max(width(DBref)), type = "sequences")
  mergedfile_align <- msa(mergedfile)
  aln <- msaConvert(mergedfile_align, type="seqinr::alignment")
  d <- dist.alignment(aln, "identity",1)
  mergedfile_mat <- as.matrix(d) 
  mergedfile_mat_x <- as.hclust(agnes(x=mergedfile_mat, diss = TRUE, stand = TRUE, metric = "euclidean", method = "average"))
  mergedfile_tree <- as.phylo(mergedfile_mat_x) 
  p <- ggtree(mergedfile_tree, size = 1, layout = "circular") 
  p
  x <- p + geom_tiplab2(color = "black", size = 0.2, offset = 0.005) 
  x 
  
  tip_labels <- as.data.frame(mergedfile_tree$tip.label)
  tip_labels <- cbind(rownames(tip_labels), tip_labels)
  colnames(tip_labels) <- c("Tips","Specimen")
  group1 <- tip_labels[grep("ASV", tip_labels$Specimen), ]
  group1 <- c(rownames(group1))

  for (j in 1:length(group1)) {
    x <- x + geom_hilight(node=group1[j], fill="red", extend = 0.05, alpha = 0.7)
    Sys.sleep(0.4)
  }

  sample_name <- gsub('_ParasiticAsvSequences.fa','', file_name) 
  pdf(paste0("GlobalTree_",sample_name,".pdf"))
  print(x) 
  dev.off()
  
##### Subtree = Focus on one contig and its closely related reference sequences 
  #convert x to phylo class:
  x <- as.phylo(x)
  x
  
  ## Identify the node in which ASV is located
  node_asv <- grep("ASV", mergedfile_tree$tip.label)
  node_asv
  
  ## ASV name
  ASV_name <- mergedfile_tree$tip.label[grepl("ASV", mergedfile_tree$tip.label)]
  ASV_name
  
## One subtree per contig sequences
  n=1
        for (node in node_asv) {
          print(node)
          treei <- tree_subset(x, node, levels_back=4) 
          ASV_name_analyzed <- ASV_name[n]
          pi <- ggtree(treei, aes(color=group)) +
            scale_color_manual(values = c("black", "red")) +
            geom_tiplab() +  xlim(0, 2) + theme_tree2() +
            theme(legend.position="none") +
            ggtitle(paste0("Subtree_",ASV_name_analyzed)) +
            theme(plot.title = element_text(color = "firebrick4", size = 12, face = "bold"))
         pdf(paste0("Subtree_",ASV_name_analyzed,".pdf"))
         print(pi) 
         dev.off()
         Sys.sleep(0.4)
         n= n+1
        }
  
## All subtrees of one sample printed in one picture
  pi <- list()
  n=1
  m=1
        for (node in node_asv) {
          print(node)
          treei <- tree_subset(x, node, levels_back=4) 
          ASV_name_analyzed <- ASV_name[m]
          pi[[n]] <- ggtree(treei, aes(color=group)) +
            scale_color_manual(values = c("black", "firebrick3")) +
            geom_tiplab() +  xlim(0, 2) + theme_tree2() +
            theme(legend.position="none") +
            ggtitle(paste0("Subtree_",ASV_name_analyzed)) +
            theme(plot.title = element_text(color = "firebrick4", size = 12, face = "bold"))
          n <- n+1
          m <- m+1
        }
        pdf(paste0("Subtrees_", sample_name, ".pdf"), height=15, width=15)
        do.call(grid.arrange,pi)
        dev.off()
        Sys.sleep(0.4)
        
}

### move phylo files to a new folder
dir.create("../6_Phylo_Output")
movefile1 <- list.files(path=getwd(), pattern= "GlobalTree_*" )
movefile1
movefile2 <- list.files(path="../5_Blastn_Output", pattern= "Subtree_*" )
movefile2
file.move(movefile1, "../6_Phylo_Output")
file.move(movefile2, "../6_Phylo_Output")


