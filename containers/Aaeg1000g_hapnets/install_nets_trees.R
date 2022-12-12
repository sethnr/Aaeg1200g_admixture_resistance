###############################################################################
# install additional R packages
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise can fail silently

repos <- c('http://cran.us.r-project.org')

nets_and_trees = c("igraph","ggnetwork","ape","phangorn","pegas","vcfR")

install.packages(nets_and_trees,repos=repos,clean=TRUE)

q(save = "no")
