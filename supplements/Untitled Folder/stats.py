import os
import re
import numpy as np
import math
import random
from matplotlib import pyplot as plt
BONF_KMERDIFF_FILE = "out_wo_bonf.kmerDiff"
UNIQUE_KMERDIFF_FILE = "out_unique.kmerDiff"
TOTAL_KMERS_FILE = "total_kmers.txt"
SORTED_FILE = "sorted_files.txt"
PHENOVALUE_FILE = "phenotype_values.txt"
TOTAL_KMER_CNT_FILE = "total_kmer_counts.txt"
MATCH_FILE = "compareFile4"


print(BONF_KMERDIFF_FILE)
os.system('./stats.sh')
stats_File = open("stats", "a")
float_regex = "[+-]?([0-9]*[.])?[0-9]+e[+-][0-9]+"

case_match = '[A-Z]+   ' + float_regex + "   " + float_regex
control_match = '[A-Z]+   ' + float_regex + "                  " + float_regex

total_line_cnt = 0
case_match_cnt = 0
control_match_cnt = 0

caseq_pval = []
case_pval = []
control_pval = []
controlq_pval = []
with open(MATCH_FILE) as f:

    for line in f:
        splt = line.split()

        total_line_cnt+=1
        if re.search(case_match, line):
            case_match_cnt+=1
            caseq_pval.append( float(splt[1]) ) 
            case_pval.append( float(splt[2]) ) 

        elif re.search(control_match, line):
            control_match_cnt+=1
            controlq_pval.append( float(splt[1]) ) 
            control_pval.append( float(splt[2]) ) 
            # print(line)
    stats_File.write('Total kmer in ' + MATCH_FILE + ' = ' + str(total_line_cnt) + '\n' )
    stats_File.write("Matched with case/control " + str(case_match_cnt + control_match_cnt) + '\n' ) 
    stats_File.write( '% of  kmer match with case  = ' + str( round( case_match_cnt/(case_match_cnt + control_match_cnt), 5) ) + '\n'  ) 
    stats_File.write( '% of  kmer match with control  = ' + str( round( control_match_cnt/(case_match_cnt + control_match_cnt), 5 ) ) + '\n' )
    stats_File.write( '% of  kmer with no match  = ' +  str(total_line_cnt - (case_match_cnt + control_match_cnt))  + '\n' )

random_size = 1000
random_sample = random.sample(list(zip(caseq_pval,case_pval)), min(random_size, len(case_pval) ))
xx,yy = list(zip(*random_sample))
cse = plt.scatter(np.negative(np.log(xx)), np.negative(np.log(yy)) , marker='x', c='r' )

random_sample = random.sample(list(zip(controlq_pval,control_pval)), min(random_size, len(control_pval)))
xx,yy = list(zip(*random_sample))
cntrl = plt.scatter(np.negative(np.log(xx)), np.negative(np.log(yy)) , marker='o', c='g' )
plt.legend((cse, cntrl), ('case','cotrol') )
plt.xlabel('-log(pval) [hawkQ]')
plt.ylabel('-log(pval) [hawk]')
plt.show()
# # plt.savefig('foo.png')

