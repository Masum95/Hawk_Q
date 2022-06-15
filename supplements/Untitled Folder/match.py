import pandas as pd
import os

# df = pd.read_csv('pvals_top_merged.txt',sep='\t', header=None)
# df = pd.read_csv('out_wo_bonf.kmerDiff',sep='\t', header=None)
import subprocess

proc = subprocess.Popen("awk -F '\t' '{print $1}' out_wo_bonf.kmerDiff",     encoding='utf8', shell=True, stdout=subprocess.PIPE)
kmers = proc.communicate()[0].strip().split('\n')
# kmers = os.popen("awk -F '\t' '{print $1}' out_wo_bonf.kmerDiff").read() #df.iloc[:,0].values

with open('5perPass','r') as f:
    for line in f:
        line = line.strip()
        for kmer in kmers:
            kmer = kmer.strip()

            if line == kmer:
                print('match')


#!/bin/bash 
# pass=`cat 5perPass`
# bonf=`awk -F '\t' '{print $1}' out_wo_bonf.kmerDiff`
# for line in pass
# do

#     found=0
#     for fs in bonf
#     do
    
#     if [[ *$fs* == *$line* ]]; then
#       found=1
#           echo $line

#     fi

#     done;

# done;
