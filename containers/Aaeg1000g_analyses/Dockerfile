#FROM ubuntu:latest
FROM continuumio/miniconda3:latest
# FROM continuumio/miniconda3:4.10.3

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y gnupg2

#install google cloud keys
RUN wget https://packages.cloud.google.com/apt/doc/apt-key.gpg
RUN cat apt-key.gpg | apt-key add -

# auto-install google-cloud-sdk from online script
# see: https://cloud.google.com/sdk/docs/downloads-interactive
RUN curl -sSL https://sdk.cloud.google.com | bash


#install conda
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh
# RUN bash ./Miniconda3-py37_4.10.3-Linux-x86_64.sh -b -p /root/miniconda
# ENV PATH="${PATH}:/root/miniconda/bin"
# RUN conda init



# install conda packages
WORKDIR /root/
# RUN conda update --all
RUN conda install mamba -c conda-forge
RUN mamba install -c bioconda -c conda-forge -y snakemake
RUN mamba install -c bioconda snakemake-wrapper-utils

RUN mamba install -c bioconda -y vcftools

# RUN mamba install -c bioconda -y pindel
# RUN pip install lostruct

#python
RUN mamba install -c bioconda pandas
RUN mamba install -c bioconda numpy

#RUN mamba install -c bioconda -y scikit-learn
RUN mamba install -c conda-forge/label/main  -y scikit-allel
RUN mamba install -c anaconda -y scikit-learn


#BCFtools / HTSlib / Samtools
RUN mamba install -c bioconda -c conda-forge -y htslib'>=1.13,<1.14.0a0'
RUN mamba install -c bioconda -c conda-forge -y bcftools=1.13
RUN mamba install -c bioconda -c conda-forge -y samtools=1.11


############
# Rstats
############
RUN mamba install -c conda-forge  -y r-base
RUN mamba install -c conda-forge  -y r-tidyverse


# INSTALL R PACKAGES:

COPY install_basics.R  install_basics.R
RUN Rscript install_basics.R

COPY install_nets_trees.R  install_nets_trees.R
RUN Rscript install_nets_trees.R

COPY install_devtools.R  install_devtools.R
RUN Rscript install_devtools.R

COPY install_demographics.R  install_demographics.R
RUN Rscript install_demographics.R


ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"
ENV LD_LIBRARY_PATH="${PATH}"
