###############################################################################
# install additional R packages
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise can fail silently

# repos <- c("http://cran.cnr.Berkeley.edu", "http://cran.mtu.edu","http://cran.r-studio.com")
repos <- c('http://cran.us.r-project.org')

required = c("tidyverse","gridExtra","ggmap",
"lubridate","knitr",
"igraph","ggnetwork",
"ape","phangorn","zoo",
"Rcpp")

install.packages(required,repos=repos,clean=TRUE)

q(save = "no")
