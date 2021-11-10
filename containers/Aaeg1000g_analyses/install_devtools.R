###############################################################################
# install additional R packages
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise can fail silently

repos <- c('http://cran.us.r-project.org')

install.packages("devtools",repos=repos,clean=TRUE)

q(save = "no")
