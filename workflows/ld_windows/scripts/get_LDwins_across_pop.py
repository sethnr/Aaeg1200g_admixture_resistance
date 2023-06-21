#!/usr/bin/env python

import allel
import numpy as np
import scipy as sp

import re
import sys

import pandas as pd

import random

from math import isnan, isfinite

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vcf',     '-v', help='genotypes (vcf)')
parser.add_argument('--meta',    '-m', help='metadata file')
parser.add_argument('--out',     '-o', help='outfile name (bed)')
parser.add_argument('--block',   '-b', default=50000, help='block size to assess (variants)')
parser.add_argument('--step',    '-s', default=50000, help='step size for windowing (variants)')
parser.add_argument('--chrom',   '-c', help='chromosome id')
parser.add_argument('--pop',     '-p', help='population name')
parser.add_argument('--poplevel','-l', help='population level (in meta file)')

args = parser.parse_args()

vcffile = args.vcf
metafile = args.meta
outfile = args.out
pop = args.pop
poplevel = args.poplevel

chrom = int(args.chrom)

bsize=int(args.block)
bstep=int(args.step)

chrlens = {1:310827022,
          2:474425716,
          3:409777670}

meta = np.genfromtxt(metafile,delimiter='\t',names=True,dtype=None,encoding='utf-8')
pop1samps = meta['sample'][meta[poplevel]==pop]

print("getting callset for {}. n={}",format(country,str(len(pop1samps))),file=sys.stderr)

callset1 = allel.read_vcf(vcffile, samples=pop1samps, fields='*')
altcounts = callset1['calldata/GT'][:,:,0] + callset1['calldata/GT'][:,:,1]

print(altcounts.shape,file=sys.stderr)

#remove all non-variant sites, and singletons
hasvar = np.logical_not(np.any(np.vstack([np.all(altcounts==2,1),
                np.all(altcounts==1,1),
                np.all(altcounts==0,1),
                np.sum(altcounts,1)==1]),0))

goodposn = callset1['variants/POS'][hasvar]
goodcounts = altcounts[hasvar,:]

print("calculating ld across {} vars",format(str(len(goodposn))),file=sys.stderr)
#get median r2 in wins
r2, r2wins, r2n  = allel.windowed_r_squared(goodposn, goodcounts, 
                                            size=bsize, start=0, stop=chrlens[chrom], 
                                            step=bsize, percentile=50)

print("printing {} windows",format(str(len(r2))),file=sys.stderr)
r2tab = pd.DataFrame({
        "country":[pop]*len(r2),
        "chrom":[chrom]*len(r2),
        "start":r2wins[:,0],
        "end":r2wins[:,1],
        "r2samp":r2})
    
r2tab.to_csv(outfile, sep="\t", index=False)
