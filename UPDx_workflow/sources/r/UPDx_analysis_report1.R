#!/usr/bin/Rscript

########################################################################
########################################################################
## UPDx Final report script

###################################################
##### Package/program Requirement
##################################################

##Packages required
library(tinytex)


## import the reference database (fasta file) 
## intercept the parameters/arguments passed by shell
args <- commandArgs()
print(args)
UPDx_Script <- args[6]

## generate report
dir.create("../7_Report_Output")
UPDx_report <- "../7_Report_Output"
samples <- scan("../samples.txt", what="character")
samples
rmarkdown::render(paste0(UPDx_Script, "/r/UPDx_analysis_report2.Rmd"), 
                  output_format = "pdf_document",
                  output_file = paste0("UPDX_report.pdf"),
                  output_dir = UPDx_report, 
                  knit_root_dir = getwd(),
                  clean = TRUE)

