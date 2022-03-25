############################
##### R packages installation
############################

#Package installation
package_list <- c("Rcpp", "readr", "tidyr", "tidyverse", "rlang", "taxize",  "filesstrings", "seqinr", "NCmisc")
for (package in package_list) {
  if(!require(package)){
    if(package == "taxize"){
      require(devtools)
      install_version("taxize", version = "0.9.5", repos = "http://cran.us.r-project.org")
      } else {
    install.packages(package, character.only = TRUE, dependencies = TRUE)
    }
  }
}
# Package loading
for (package in package_list) {
  library(package, character.only = TRUE)
} 

