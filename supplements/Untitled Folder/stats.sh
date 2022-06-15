#!/bin/bash 


BONF_KMERDIFF_FILE="out_wo_bonf.kmerDiff"
UNIQUE_KMERDIFF_FILE="out_unique.kmerDiff"
TOTAL_KMERS_FILE="total_kmers.txt"
SORTED_FILE="sorted_files.txt"
PHENOVALUE_FILE="phenotype_values.txt"
TOTAL_KMER_CNT_FILE="total_kmer_counts.txt"

stats_File="stats"
echo "Hawk_q" > $stats_File
wc -l  $BONF_KMERDIFF_FILE >> $stats_File
wc -l $UNIQUE_KMERDIFF_FILE >> $stats_File



