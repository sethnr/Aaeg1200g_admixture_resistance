#!/usr/bin/env python

import allel
import numpy as np
import re
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--metafile', '--meta', '-m',  help='metadata file (must have "sample" column)')
parser.add_argument('--chromfile', '-C',  help='chrom file (must have "chrom,length,chrom_id" columns)')
parser.add_argument('--chrom', '-c', help='chromosome number')
parser.add_argument('--column', '-M', default='country', help='metadata column to search in')
parser.add_argument('--vcf', '-v', help='vcf file to test')
parser.add_argument('--block', '-b', default=1000000, help='block size to assess (variants)')
#parser.add_argument('--step', '-s', default=1000, help='step size for windowing (variants)')
parser.add_argument('--out', '-o',  help='out file prefix')

args = parser.parse_args()

metafile = args.metafile
chromfile = args.chromfile
column = args.column
chrom = int(args.chrom)
vcf = args.vcf
out = args.out
block = int(args.block)
#step = int(args.step)


meta = np.genfromtxt(metafile,delimiter='\t',names=True,dtype=None,encoding='utf-8')

chroms = np.genfromtxt(chromfile,delimiter='\t',names=True,dtype=None,encoding='utf-8')

clen = chroms['length'][chroms['chrom']==chrom][0]
cname = chroms['chrom_id'][chroms['chrom']==chrom][0]
print(chromfile,chroms['chrom'],chrom,clen,cname)
allsamps = meta['sample']


#make snp and geno files for eigenstrat
genofile = open(out+'.geno', 'a')
snpfile = open(out+'.snp', 'a')

for S in range(1,clen,block):
    region='{}:{}-{}'.format(cname,S,S+(block-1))
    # print(region)
    callset = allel.read_vcf(vcf, region=region, fields='*')
    #callset = allel.read_vcf(vcf, region=region, samples=allsamps, fields='*')
    if callset is not None:
        gtpop = allel.GenotypeArray(callset['calldata/GT'])
        acpop = gtpop.count_alleles()
        refcounts = gtpop.to_allele_counts(0).astype(str)[:,:,0]
        refcounts[gtpop.is_missing()]=" "
        np.savetxt(genofile, refcounts, fmt='%s',delimiter='')

        ziparray = [(n,c,m,p) for n,c,m,p in zip([str(chrom)+":"+str(I) for I in callset['variants/POS']],
                    [chrom]*len(callset['variants/POS']),
                    callset['variants/CM'][:,0],
                    callset['variants/POS']
                )]
        snparray = np.array(ziparray,dtype=[('snpid', 'U15'), ('chrom', 'int'), ('cM', 'int'), ('bp', 'int')])
        np.savetxt(snpfile, snparray, fmt='%s',delimiter='\t')

genofile.close()
snpfile.close()



#make sample file from meta file
sampfile = open(out+'.ind', 'w')

#alder failing due to long names!
shortsamples = [re.sub("Debug.*aegypti_","",S) for S in meta['sample']]

#sampzip = [(s,x,c) for s,x,c in zip(meta['sample'].tolist(),
sampzip = [(s,x,c) for s,x,c in zip(shortsamples,
                           ['U']*len(allsamps),
                           meta[column].tolist())]
samparray = np.array(sampzip,dtype=[('sample', 'U50'), ('sex', 'U2'), ('class', 'U30')])
np.savetxt(sampfile, samparray, fmt='%s',delimiter='\t')
sampfile.close()
