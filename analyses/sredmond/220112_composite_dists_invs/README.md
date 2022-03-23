# SNP genotypes

This resource contains SNP genotypes for samples grouped by sample set. 

`all` provides data for all non-N sites across the AgamP4 reference genome.

Samples have been filtered for quality. The order of samples in the arrays is consistent with the ordering in `metadata`.

The data are provided in Zarr format and are organised in the same way as the scikit-allel [vcf_to_zarr()](https://scikit-allel.readthedocs.io/en/stable/io.html#allel.vcf_to_zarr) function, where fields equivalent to VCF FORMAT (pertaining to samples) are stored in arrays within the `calldata` group.

`all/{sample_set}` provides data for a single sample set, with a `calldata/` group only. Arrays provided are:
`GT` genotype. (nSites x nSamples x ploidy (2))
`AD` allele depth. (nSites x nSamples x nAlleles (4))
`MQ` mapping quality. (nSites x nSamples)
`GQ` genotype quality. (nSites x nSamples)

To perform site filtering, the `site_filters` arrays from the `v3` release must be used.

For data that pertains to variants, which are common between sample sets, see `all/sites` in the `v3` release.

For further information or queries contact data@malariagen.net.
