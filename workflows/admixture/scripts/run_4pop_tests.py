#!/usr/bin/env python


import allel
import numpy as np

import re
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--pop1', '-p1', help='population 1')
parser.add_argument('--pop2', '-p2', help='population 2')
parser.add_argument('--pop3', '-p3', help='population 3 (admixed)')
parser.add_argument('--pop4', '-p4', help='population 3 (admixed)')
parser.add_argument('--metafile', '--meta', '-m',  help='metadata file (must have "sample" column)')
parser.add_argument('--chrom', '-c', help='chromosome number')
parser.add_argument('--column', '-M', default='country', help='metadata column to search in')
parser.add_argument('--vcf', '-v', help='vcf file to test')
parser.add_argument('--out', '-o', help='vcf file to test')
parser.add_argument('--block', '-b', default=500000, help='block size to assess (variants)')
parser.add_argument('--step', '-s',  default=250000, help='step size for windowing (variants)')

args = parser.parse_args()

pop1 = args.pop1
pop2 = args.pop2
pop3 = args.pop3
pop4 = args.pop4
metafile = args.metafile
column = args.column
chrom = args.chrom
vcf = args.vcf
out = args.out
if out is None:
    out = "4pop_{}_{}_{}_{}".format(pop1,pop2,pop3,pop4)
block = int(args.block)
step = int(args.step)


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

pop4samps = meta['sample'][meta[column]==pop4]
print("parsing {} samples from {}".format(pop4samps.size,pop4),file=sys.stderr)
callset4 = allel.read_vcf(vcf, samples=pop4samps, fields='*')
gtpop4 = allel.GenotypeArray(callset4['calldata/GT'])
acpop4 = gtpop4.count_alleles()

print("testing if {} is admixed from {} and {} with outgroup {}".format(pop2,pop1,pop3,pop4),file=sys.stderr)

f4, f4se, f4z, f4blk, f4jack = allel.average_patterson_d(acpop1, acpop2, acpop3, acpop4, block)
#f4wins = allel.moving_patterson_f4(acpop1, acpop2, acpop3, acpop4, block, step=step, normed=True)


# poswins = np.vstack((
#             [pop3] * len(f4wins),
#             [pop1] * len(f4wins),
#             [pop2] * len(f4wins),
#             [chrom] * len(f4wins),
#             allel.moving_statistic(callset1['variants/POS'], min, block, step=step),
#             allel.moving_statistic(callset1['variants/POS'], np.average, block, step=step),
#             allel.moving_statistic(callset1['variants/POS'], max, block, step=step),
#             f3wins)).transpose()

sumfile=open(out+"_summary.txt",'w')
zlim=3

result = ""
signif=False
if(f4z>zlim):
    result = "{}->{}".format(pop1,pop3)
    signif=True
elif(f4z<(zlim*-1)):
    result = "{}->{}".format(pop2,pop3)
    signif=True
print("\t".join(map(str,[pop1,pop2,pop3,pop4,chrom,f4,f4se,f4z,signif,result])), file=sumfile)
sumfile.close()

#np.savetxt(out+"_blocks.txt", poswins, fmt='%s', delimiter='\t')
