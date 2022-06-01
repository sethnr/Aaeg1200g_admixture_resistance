FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y gnupg2


#alder / eagle / fs prerequisites:
RUN apt-get install -y fftw3
RUN apt-get install -y libgomp1
#RUN apt-get install -y liblapack3 libatlas3-base
#RUN apt-get install -y liblapack3gf liblapack-doc liblapack-dev


#install google cloud keys
RUN wget https://packages.cloud.google.com/apt/doc/apt-key.gpg
RUN cat apt-key.gpg | apt-key add -

# auto-install google-cloud-sdk from online script
# see: https://cloud.google.com/sdk/docs/downloads-interactive
RUN curl -sSL https://sdk.cloud.google.com | bash


#install conda (min py39)
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#RUN bash ./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p /root/miniconda
#ENV PATH="${PATH}:/root/miniconda/bin"
RUN bash ./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p /opt/miniconda

#add library paths (not found when run as singularity image)

ENV PATH="${PATH}:/opt/miniconda/bin:/opt/miniconda/lib"

RUN conda init

# install conda packages
#WORKDIR /root/
RUN conda install  -c conda-forge mamba
RUN mamba install -c bioconda -c conda-forge -y snakemake-minimal
RUN mamba install -c bioconda -c conda-forge -y bcftools
RUN mamba install -c bioconda -c conda-forge -y htslib'>=1.13,<1.14.0a0'

#phasing software
RUN mamba install -c bioconda -c conda-forge -y whatshap nomkl

RUN mamba install -c bioconda  -c conda-forge -y libgcc-ng
RUN mamba install -c conda-forge -y boost-cpp
RUN mamba install -c bioconda  -c conda-forge -y boost
RUN mamba install -c bioconda  -c conda-forge -y libstdcxx-ng
RUN mamba install -c bioconda  -c conda-forge -y bzip2
RUN mamba install -c bioconda  -c conda-forge -y shapeit4

#eagle / alder / fs prerequisites
#RUN mamba install -c conda-forge -y libgomp
#RUN mamba install -c conda-forge -y fftw
RUN mamba install -c conda-forge -y liblapack
RUN mamba install -c conda-forge -y openmp

RUN mamba install -c bioconda alder

RUN ln -s /opt/miniconda/lib/liblapack.so.3 /opt/miniconda/lib/liblapack.so.3gf
COPY eagle /usr/bin/eagle

#admixture software
COPY fs /usr/bin/fs
# COPY alder /usr/bin/alder
