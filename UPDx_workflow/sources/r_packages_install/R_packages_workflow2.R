############################
##### R packages installation
############################


package_list1 <- c("ggplot2", "RColorBrewer", "reshape", "dbscan", "gridExtra", "stringr", "cluster", "phylogram",  "ggdendro", 
                   "colorspace", "ape", "Biostrings", "dplyr") 

package_list2 <- c("msa", "DECIPHER", "treeio") 

package_list3 <- c("ggtree") 

#Package installation
for (package in package_list1) {
  if(!require(package)){
    install.packages(package, character.only = TRUE, dependencies = TRUE)
    }
}

#Package installation via BiocManager
for (package in package_list2) {
  if(!require(package)){
        if (!requireNamespace("BiocManager", quietly=TRUE)) {
        install.packages("BiocManager")
        BiocManager::install(package)
    } else {
      BiocManager::install(package)
    }
  }
}

#Package installation via devtools
if (!require("devtools")) {
  install.packages("devtools")
} else {
  if(!require("ggtree")){
    devtools::install_github("YuLab-SMU/ggtree")
  }
}


# Package loading
for (package in package_list1) {
  library(package, character.only = TRUE)
} 
for (package in package_list2) {
  library(package, character.only = TRUE)
} 
library("ggtree")
library("Rtsne")


