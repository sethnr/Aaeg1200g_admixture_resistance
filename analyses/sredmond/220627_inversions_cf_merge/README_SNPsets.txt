Aaeg1200G SNP sets (Jan 12 2022)

SNPs were called using ANGSD (snp pval 1e-6) across all non-African populations, and kept only if present at sites with >80% of individuals with data and minor allele present on greater than 10 chromosomes.  Then this full SNP set was thinned using PLINK with the following three logics, assuming that all or most SNPs segregating in non-African populations would also be segregating in African populations.s.

Three different SNP sets from PLINK. 
1. LD-pruned (n = 7,161,365, aka 7e6):  Unlinked SNPs (maf>0.01) within a window size of 100 SNPs, step size of 10, and cutoff or R2<0.1.
2. LD-pruned, thinned to 1M (n = 1,000,000 aka 1e6): The first set (7e6) randomly thinned down to 1e6 SNPs.
3. LD-pruned, thinned to 100k (n = 100,000 aka 100k): The first set (7e6) randomly thinned down to 100k SNPs.

We have a variety of files representing these SNPs.  
1) sites files for use in ANGSD or other.  These files include the major and minor allele at each site: 
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.100k.sites*
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.1e6.sites*
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.7e6.sites*

2) VCF files (including African and non-African populations, genotypes and genotype likelihoods, chromosomes labeled as 1,2,3):   
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr123.SNP.100k.vcf.gz
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr123.SNP.1e6.vcf.gz
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr123.SNP.7e6.vcf.gz

3) VCF files (including African and non-African populations, genotypes and genotype likelihoods, NC_ chromosome labels):   
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr123.SNP.100k.NC_labels.vcf.gz
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr123.SNP.1e6.NC_labels.vcf.gz
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr123.SNP.7e6.NC_labels.vcf.gz

4) Beagle-gl files (including African and non-African populations, only unrelated inds from Afr pops, genotype likelihoods, NC_ chromosome labels):   
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.non-AfrI3.AfrI3.unrel.n1206.100k.beagle.gz 
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.non-AfrI3.AfrI3.unrel.n1206.1e6.beagle.gz

5) PLINK files (including African and non-African populations, only unrelated inds from Afr pops, chromosomes labeled as 1,2,3):   
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.non-AfrI3.AfrI3.unrel.n1206.100k.ped (.map)
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.non-AfrI3.AfrI3.unrel.n1206.1e6.ped (.map) 
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/aegy.wgs.non-AfrI3.AfrI3.unrel.n1206.7es.ped (.map) 

6) BED interval files 
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.100k.intervals.bed
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.1e6.intervals.bed
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.7e6.intervals.bed

7) plink binary files. (only non-African populations)
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.100k.bed (and associated .bim etc. files)
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.1e6.bed (and associated .bim etc. files)
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.7e6.bed (and associated .bim etc. files)

8) VCF files (only non-African population, genotypes only (no likelihoods))
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.100k.vcf.gz
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.1e6.vcf.gz
gs://verily-aaeg1200g/vcfs/SNPs/pruned_sites/non-AfrI3.unlinked.maf01.7e6.vcf.gz


## ANGSD command for SNP calls
## Non-African pops
REF=/dglob/databases/GCF_002204515.2_AaegL5.0_genomic_Afr_male_altref_iter3.fna
THREADS=2
SITES=/dglob/data/sites/aegy.wgs.AfrI3.mapQ20.GoodSites.chr123.bcf.sites
OUTGENO=/dglob/data/snps/non-afr.regs/aegy.wgs.non-AfrI3.CL.jntSites
BAM_FILES=/dglob/data/non-AfrI3.CL.bam.lists/aegy.wgs.non-AfrI3.CL.bams

ls /dglob/data/regs/*.reg | parallel -j 64 angsd -bam $BAM_FILES -checkBamHeaders 0 -P $THREADS -fai $REF.fai -ref $REF -rf {} -sites $SITES -out $OUTGENO.{/} -minMapQ 20 -minQ 10 -remove_bads 1 -uniqueOnly 0 --ignore-RG 0 -gl 1 -domaf 1 -doMajorMinor 3 -SNP_pval 1e-6
