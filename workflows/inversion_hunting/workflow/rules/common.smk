import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.utils import min_version

from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

min_version("5.18.0")


container: "docker://sethnr/aaeg1000g_analysis:0.02"


###### Config file and sample sheets #####
configfile: "config_lostruct/config.yaml"




# wildcard_constraints:
#     chrom="|".join(vcfs.index.astype("string").unique())




##### Helper functions #####

# def get_remote_vcfs(wildcards):
#     """Get vcfs of given chrom"""
#     return [GS.remote(v) for v in vcfs.loc[(wildcards.chrom)].vcf]

def get_remote_vcf(wildcards):
    vcf = pd.read_table(config["vcfs"],dtype = str).set_index(["chrom","block"])
    """Get vcf of given chrom and block """
    return GS.remote(vcf.loc[wildcards.chrom,wildcards.block].vcf)

# def get_chrom_vcfs(wildcards):
#     """Get vcfs of given chrom"""
#     return vcfs.loc[(wildcards.chrom)].vcf

# def get_chrom_blocks(wildcards):
#     """Get blocks for given chrom"""
#     return vcfs.loc[(wildcards.chrom)].block

def get_chrom_blocks_string(wildcards):
    vcfs = pd.read_table(config["vcfs"],dtype = str).set_index("chrom")
    """Get block strings for given chrom"""
    return expand("chr{chrom}_blk{block}",chrom=wildcards.chrom,
                                          block=vcfs.loc[(wildcards.chrom)].block)

# def get_chrom_blocks_spp_string(wildcards):
#     vcfs = pd.read_table(config["vcfs"],dtype = str).set_index("chrom")
#     """Get block strings for given chrom"""
#     return expand("chr{chrom}_blk{block}_spp{spp}",chrom=wildcards.chrom,
#                                           block=vcfs.loc[(wildcards.chrom)].block,
#                                           spp=wildcards.spp)




rule tabix_vcf:
    input:
        "results_lostruct/{calltype}/{vcfname}.vcf.gz"
    output:
        "results_lostruct/{calltype}/{vcfname}.vcf.gz.tbi",
    log:
        "logs/tabix/{calltype}/{vcfname}.vcf.log",
    params:
        "-p vcf",
    wrapper:
        "0.74.0/bio/tabix"
