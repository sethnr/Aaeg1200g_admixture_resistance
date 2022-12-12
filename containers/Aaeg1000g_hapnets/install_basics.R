###############################################################################
# install additional R packages
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise can fail silently

repos <- c('http://cran.us.r-project.org')

# required = c("gridExtra","zoo","devtools")
required = c("gridExtra","zoo","getopt","patchwork","pvclust")
install.packages(required,repos=repos,clean=TRUE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")


BiocManager::install("GenomicRanges",update=F)
BiocManager::install("Biostrings",update=F)
BiocManager::install("rtracklayer",update=F)
BiocManager::install("trackViewer",update=F)
install.packages("geneHapR")



q(save = "no")
