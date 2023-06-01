import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.utils import min_version

from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

min_version("5.18.0")


##### Helper functions #####

# def get_remote_vcfs(wildcards):
#     """Get vcfs of given chrom"""
#     return [GS.remote(v) for v in vcfs.loc[(wildcards.chrom)].vcf]

def get_remote_bam_url(wildcards):
    bams = pd.read_table(config["bams"],dtype = str).set_index(["sample"])
    return bams.loc[wildcards.sample].bam

def get_remote_bai_url(wildcards):
    bams = pd.read_table(config["bams"],dtype = str).set_index(["sample"])
    return bams.loc[wildcards.sample].bai

def get_local_bam_url(wildcards):
    bams = pd.read_table(config["bams"],dtype = str).set_index(["sample"])
    return bams.loc[wildcards.sample].bam.replace("gs://","")

def get_local_bai_url(wildcards):
    bams = pd.read_table(config["bams"],dtype = str).set_index(["sample"])
    return bams.loc[wildcards.sample].bai.replace("gs://","")



def get_remote_bam(wildcards):
    bams = pd.read_table(config["bams"],dtype = str).set_index(["sample"])
    """Get bam of given sample """
    return GS.remote(bams.loc[wildcards.sample].bam, keep_local=False)

def get_remote_bai(wildcards):
    bams = pd.read_table(config["bams"],dtype = str).set_index(["sample"])
    """Get bam of given sample """
    return GS.remote(bams.loc[wildcards.sample].bai, keep_local=False)


def get_remote_vcf(wildcards):
    vcfs = pd.read_table(config["vcfs"],dtype = str,header=0).set_index(["chrom"])
    """Get processed/thinned vcf of given chrom """
    return GS.remote(vcfs.loc[wildcards.chrom].vcf, keep_local=True)

def get_remote_vcf_tbi(wildcards):
    vcfs = pd.read_table(config["vcfs"],dtype = str,header=0).set_index(["chrom"])
    """Get processed/thinned vcf of given chrom """
    vcf = vcfs.loc[wildcards.chrom].vcf
    return GS.remote(vcf+".tbi", keep_local=True)


def get_remote_vcf_raw(wildcards):
    vcf = pd.read_table(config["vcfs"],dtype = str).set_index(["chrom","block"])
    """Get vcf of given chrom and block """
    return GS.remote(vcf.loc[wildcards.chrom,wildcards.block].vcf)

def get_chrom_id(wildcards):
    chromtab = pd.read_table(config["chrommap"],dtype = str).set_index(["chrom"])
    """lookup ncbi chrom ID for chrom No"""
    return chromtab.loc[str(wildcards.chrom)].chrom_id



def get_chrom_blocks_string(wildcards):
    vcfs = pd.read_table(config["vcfs"],dtype = str).set_index("chrom")
    """Get block strings for given chrom"""
    return expand("chr{chrom}_blk{block}",chrom=wildcards.chrom,
                                          block=vcfs.loc[(wildcards.chrom)].block)



def get_bam_for_sample(wildcards):
    samples = pd.read_table(config["samples"],dtype = str).set_index(["sample"])
    """Get vcf of given chrom and block """
    return GS.remote(samples.loc[wildcards.sample].bam)

def get_bai_for_sample(wildcards):
    samples = pd.read_table(config["samples"],dtype = str).set_index(["sample"])
    """Get vcf of given chrom and block """
    return GS.remote(samples.loc[wildcards.sample].bai)


rule tabix_vcf:
    input:
        "results/{vcftype}/{vcfname}.vcf.gz"
    output:
        "results/{vcftype}/{vcfname}.vcf.gz.tbi",
    log:
        "logs/tabix/{vcftype}/{vcfname}.vcf.log",
    params:
        "-p vcf",
    wrapper:
        "0.74.0/bio/tabix"
