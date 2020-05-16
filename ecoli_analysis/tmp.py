import pandas as pd 

lst = []
with open('deletedFiles') as dlt:
	for line in dlt:
		lst.append( line.strip())

print(lst)

toKeep = ''
toKeepCnt = ''
with open('sorted_files.txt') as sorted, open('total_kmer_counts.txt','r') as cntFile:
	for line,cnt in zip(sorted, cntFile):
		fle = line.split('/')[-2]
		if fle not in lst:
			toKeep += line 
			toKeepCnt += cnt

with open('toKeep','w') as tokeep, open('toKeepCnt','w') as toknt:
	tokeep.write(toKeep)
	toknt.write(toKeepCnt)

 
