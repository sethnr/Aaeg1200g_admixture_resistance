FROM ubuntu:latest

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

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh
RUN bash ./Miniconda3-py37_4.12.0-Linux-x86_64.sh -b -p /root/miniconda
ENV PATH="${PATH}:/root/miniconda/bin"
RUN conda init


# install conda packages
WORKDIR /root/
RUN conda install  -c conda-forge mamba
RUN conda install -c bioconda -c conda-forge -y snakemake-minimal
#RUN conda install -c bioconda -c conda-forge -y snakemake-wrapper-utils
#RUN conda install -c bioconda/labels/cf201901  -c conda-forge -y bcftools
RUN conda install -c bioconda  -c conda-forge -y bcftools
RUN conda install -c bioconda -c conda-forge -y htslib'>=1.13,<1.14.0a0'

#phasing software
RUN conda install -c bioconda -c conda-forge whatshap nomkl
RUN conda install -c bioconda  -c conda-forge -y libgcc-ng
RUN conda install -c conda-forge -y boost-cpp
RUN conda install -c conda-forge/label/cf202003 boost
RUN conda install -c bioconda  -c conda-forge -y libstdcxx-ng
RUN conda install -c bioconda  -c conda-forge -y bzip2
RUN mamba install -c bioconda -c conda-forge -y shapeit4
ADD Eagle_v2.4.1.tar.gz eagle

#admixture software
ADD fs_4.1.1.tgz finestructure
ADD alder_v1.03.tar.gz alder
