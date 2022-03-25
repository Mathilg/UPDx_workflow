############################
##### R packages installation
############################

package_list1 <- c("tinytex", "knitr", "stringr", "dplyr", "readr", "tibble")

package_list2 <- c("kableExtra")

# Package installation
for (package in package_list1) {
  if(!require(package)){
    install.packages(package, character.only = TRUE, dependencies = TRUE)
    }
}



# Package installation via devtools
# kableExtra installation require libfontconfig1-dev package installion via apt
if (!require("devtools")) {
  install.packages("devtools")
} else {
  if(!require("kableExtra")){
    devtools::install_github("haozhu233/kableExtra")
  }
}

# Package installation via tinytex
if(!require("tinytex")){
    tinytex::install_tinytex()
  }


tinytex::install_tinytex()

# Loading
for (package in package_list1) {
  library(package, character.only = TRUE)
} 

for (package in package_list2) {
  library(package, character.only = TRUE)
}