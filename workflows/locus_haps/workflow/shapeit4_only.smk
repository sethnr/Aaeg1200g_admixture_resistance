#from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
#GS = GSRemoteProvider()



container: "docker://sethnr/aaeg1000g_admixture:0.02"
#### Config file and sample sheets #####
configfile: "config/config.yaml"
#### reused methods for grabbing GS objects, etc ####
include: "rules/common.smk"

CHROMS=[1,2,3]

wildcard_constraints:
    chrom=CHROMS

rule all:
    input:
        expand("results/phase_shapeit/phased_shapeit8c_{chrom}.vcf.gz",chrom=CHROMS)


rule phase_shapeit:
    input:
        vcf="results/phase_whatshap/whatshap_{chrom}.vcf.gz",
        tbi="results/phase_whatshap/whatshap_{chrom}.vcf.gz.tbi",
        map=config['map']
    output:
        vcf=protected("results/phase_shapeit/phased_shapeit8c_{chrom}.vcf.gz"),
    resources:
        mem_mb=32000,
        disk_mb=1000000,
        h_rt='120:00:00',
        bamspace=0,
        cpus=8,
	#qsub_args='-l highmem -l nodes=10:srm'
    threads: 8
    params:
        chromid=lambda wildcards: get_chrom_id(wildcards)
    log:
        "logs/phase_shapeit/shapeit_chr{chrom}.log",
    shell:
        """
	#ulimit -v 60000000
        shapeit4 --input {input.vcf} \
		 --map {input.map} \
		 --effective-size 1000000 \
                 --use-PS 0.0001 \
		 --sequencing \
		 --region {params.chromid} \
                 --thread {threads} \
	         --output {output.vcf} >> {log}
        	 # --log {log} \
        """
