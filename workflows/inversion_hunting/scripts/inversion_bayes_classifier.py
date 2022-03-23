#!/usr/bin/env python

import numpy as np
import numpy.lib.recfunctions as rf

from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_digits

import argparse
import sys
import re

import pandas as pd


parser = argparse.ArgumentParser(description='get allele numbers table')

parser.add_argument('-s','--snps', action="store", dest='snps', type=str, help='012 file of SNPs for classifier training', nargs='?', default=None)
parser.add_argument('-t','--classes', action="store", dest='targets', type=str, help='target values', nargs='?', default=None)
parser.add_argument('-p','--testprop', action="store", dest='testprop', type=int, help='save every Nth value for testing', nargs='?', default=10)
parser.add_argument('-b','--batch', action="store_true", dest='batch', help='output tab delimited text', default=False)
parser.add_argument('-o','--output', action="store", dest='predsout', type=int, help='output file (inversion predictions)', nargs='?', default=10)

args = parser.parse_args()


#
# class args:
#     pass
# args.snps = "inv_Senegal_n100_aim_snps.txt"
# args.targets = "inv_Senegal_n100_mean_calls_TRAINING.txt"
# args.testprop = 0.1
# args.batch = True
# args.predsout = "inv_Senegal_n100_nbayes_predictions.txt"



idvars = ['chrom','pos',
             'i','inv','country',
             'assoc','qual']

aimsnps = np.genfromtxt(args.snps,dtype=None,names=True)
samples = aimsnps.dtype.names[7:]
aims = aimsnps[idvars]
snps = np.genfromtxt(args.snps,dtype='int',usecols=np.arange(7,len(samples)+7),skip_header=1)
snps = np.transpose(snps)
snps.shape


#get inversion targets:
targets = np.genfromtxt(args.targets,dtype=None,names=True,
                        delimiter="\t")

#identify samples from training set in SNP table
train = np.isin(samples,targets['sample'])
test = np.logical_not(np.isin(samples,targets['sample']))


#make empty table for predictions
invpreds = np.empty((len(samples),len(targets.dtype.names[2:])),dtype="float")
invpreds[:] = -1

i=-1
for invcode in targets.dtype.names[2:]:
    i+=1
    invsnps = [re.sub('\D','',i.decode('utf-8')) for i in aims['inv']]
    strn = snps[train,:]
    ttrn = targets[invcode]

    #notNA = [i!='NA' for i in ttrn]
    #strn = strn[notNA,:]
    #ttrn = ttrn[notNA]

    stst = snps[test,:]
    #ttst = [i.decode('utf-8') for i in targets[invcode]]

    classifier = GaussianNB()
    classifier.fit(strn, ttrn)

    pred = classifier.predict(stst)
    invpreds[np.logical_not(np.isin(samples,targets['sample'])),i] = pred


# In[404]:



#save and plot in new R script...n
np.savetxt(args.predsout,invpreds)
