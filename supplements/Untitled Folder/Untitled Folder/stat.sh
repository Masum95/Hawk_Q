BONF_KMERDIFF_FILE="out_wo_bonf.kmerDiff"
UNIQUE_KMERDIFF_FILE="out_unique.kmerDiff"
TOTAL_KMERS_FILE="total_kmers.txt"
SORTED_FILE="sorted_files.txt"
PHENOVALUE_FILE="phenotype_values.txt"
TOTAL_KMER_CNT_FILE="total_kmer_counts.txt"

echo "Hawk Q" > stats
echo  "Total number of kmers  -> " `wc -l $BONF_KMERDIFF_FILE` >>stats
echo  "Total number of kmers  -> " `wc -l $UNIQUE_KMERDIFF_FILE` >>stats
# echo "Ratio : `expr `cat $BONF_KMERDIFF_FILE | wc -l ` + `cat $BONF_KMERDIFF_FILE | wc -l ` `"
#printf "%.4f\n" $((`cat $BONF_KMERDIFF_FILE | wc -l` / `cat $UNIQUE_KMERDIFF_FILE | wc -l `))
echo "Ratio :  $(calc $((`cat $BONF_KMERDIFF_FILE | wc -l` / `cat $UNIQUE_KMERDIFF_FILE | wc -l `)) )"  >> stats
# echo "Hawk " >> stats 
# echo  "Total number of kmers  -> " `wc -l case_nout_wo_bonf.kmerDiff` >>stats
# echo  "Total number of kmers  -> " `wc -l control_out_wo_bonf.kmerDiff` >>stats


# echo  "Total number of kmers  -> " `wc -l case_out_w_bonf.kmerDiff` >>stats
# echo  "Total number of kmers  -> " `wc -l control_out_w_bonf.kmerDiff` >>stats
