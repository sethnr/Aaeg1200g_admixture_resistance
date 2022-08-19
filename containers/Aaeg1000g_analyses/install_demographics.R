###############################################################################
# install additional R packages
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise can fail silently

repos <- c('http://cran.us.r-project.org')

# demographics = c("related")
# install.packages(demographics,repos=repos,clean=TRUE)

devtools::install_github("petrelharp/local_pca/lostruct")
install.packages("genetics")
q(save = "no")
