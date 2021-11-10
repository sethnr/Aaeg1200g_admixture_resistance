###############################################################################
# install additional R packages
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise can fail silently

repos <- c('http://cran.us.r-project.org')

required = c("Rcpp")
install.packages(required,repos=repos,clean=TRUE)

q(save = "no")
