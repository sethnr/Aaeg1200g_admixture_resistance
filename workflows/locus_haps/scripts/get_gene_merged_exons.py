#!/usr/bin/env python

import allel
import numpy as np

import re
import sys
import os

import argparse


# gff = "refs/VectorBase-59_AaegyptiLVP_AGWG.refseq.gff"
# genes = allel.gff3_to_recarray(gff,attributes=["gene_id"])
# geneid = 'AAEL013277'
# outfile = "merged_cds.kdr.bed"

parser = argparse.ArgumentParser()
parser.add_argument('--gff', '-g', help='genes (gff3)')
parser.add_argument('--gene', '-n', help='gene name')
parser.add_argument('--out', '-i', help='outfile name (bed)')

args = parser.parse_args()
gff = args.gff
geneid = args.gene
outfile = args.out


exons = genes[np.logical_and(genes['type']=='exon', genes['gene_id']==geneid)]


ranges = [i for i in zip(exons['start'],exons['end'])]
j=0
s=0
merged=[]
for i in range(len(ranges)):
    if ranges[i][0]<ranges[j][1]:
        continue
    else:
        j = j + 1
        merged.append((ranges[s][0], ranges[j][1]))
        s = i

outbed = np.array([i for i in zip([exons['seqid'][1]] * len(merged),
            [i[0] for i in merged],
            [i[1] for i in merged],
            [exons['gene_id'][1]] * len(merged))],
        dtype=[("seqid","<U15"),("start",int),("end",int),("gene","<U15")])

np.savetxt(outfile,outbed,delimiter="\t",fmt="%s")


# In[ ]:
