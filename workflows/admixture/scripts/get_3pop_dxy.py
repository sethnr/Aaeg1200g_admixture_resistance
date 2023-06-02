#!/usr/bin/env python


import allel
import numpy as np

import re
import sys

# from tabulate import tabulate
# import pandas as pd
# from plotnine import *


#vcf='phased_shapeit_1_200samples.vcf.gz'
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--pop1', '-p1', help='population 1')
parser.add_argument('--pop2', '-p2', help='population 2')
parser.add_argument('--pop3', '-p3', help='population 3 (admixed)')
parser.add_argument('--metafile', '--meta', '-m',  help='metadata file (must have "sample" column)')
parser.add_argument('--chrom', '-c', help='chromosome number')
parser.add_argument('--column', '-M', default='country', help='metadata column to search in')
parser.add_argument('--vcf', '-v', help='vcf file to test')
parser.add_argument('--out', '-o', help='vcf file to test')
parser.add_argument('--block', '-b', default=10000, help='block size to assess (variants)')
parser.add_argument('--step', '-s', default=1000, help='step size for windowing (variants)')

args = parser.parse_args()

pop1 = args.pop1
pop2 = args.pop2
pop3 = args.pop3
metafile = args.metafile
column = args.column
chrom = args.chrom
vcf = args.vcf
out = args.out
if out is None:
    out = "3pop_{}_{}_{}".format(pop1,pop2,pop3)
block = args.block
step = args.step


meta = np.genfromtxt(metafile,delimiter='\t',names=True,dtype=None,encoding='utf-8')

pop1samps = meta['sample'][meta[column]==pop1]
print("parsing {} samples from {}".format(pop1samps.size,pop1),file=sys.stderr)
callset1 = allel.read_vcf(vcf, samples=pop1samps, fields='*')
gtpop1 = allel.GenotypeArray(callset1['calldata/GT'])
acpop1 = gtpop1.count_alleles()

pop2samps = meta['sample'][meta[column]==pop2]
print("parsing {} samples from {}".format(pop2samps.size,pop2),file=sys.stderr)
callset2 = allel.read_vcf(vcf, samples=pop2samps, fields='*')
gtpop2 = allel.GenotypeArray(callset2['calldata/GT'])
acpop2 = gtpop2.count_alleles()

pop3samps = meta['sample'][meta[column]==pop3]
print("parsing {} samples from {}".format(pop3samps.size,pop3),file=sys.stderr)
callset3 = allel.read_vcf(vcf, samples=pop3samps, fields='*')
gtpop3 = allel.GenotypeArray(callset3['calldata/GT'])
acpop3 = gtpop3.count_alleles()

print("testing if {} is admixed from {} and {}".format(pop3,pop1,pop2),file=sys.stderr)

f3, f3sd, f3z, f3blk, f3jack = allel.average_patterson_f3(acpop3, acpop1, acpop2, block, normed=True)
#f3wins = allel.moving_patterson_f3(acpop3, acpop1, acpop2, block, step=step, normed=True)
dxy13, wins13 = allel.windowed_divergence(callset1['variants/POS'],acpop3, acpop2, block, step=step, normed=True)
dxy23, wins23 = allel.windowed_divergence(callset1['variants/POS'],acpop3, acpop1, block, step=step, normed=True)


poswins = np.vstack((
            [pop3] * len(f3wins),
            [pop1] * len(f3wins),
            [pop2] * len(f3wins),
            [chrom] * len(f3wins),
            wins13,
            dxy13,
            wins12,
            dxy12)).transpose()

np.savetxt(out+"_dxy.txt", poswins, fmt='%s', delimiter='\t')
# In[ ]:
