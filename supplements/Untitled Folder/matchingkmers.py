#!/usr/bin/python3  
import sys  
import collections

# print('Number of arguments:', len(sys.argv), 'arguments.') 
# print('Argument List:', str(sys.argv))
if sys.argv.__len__() < 5:
    print("Files?")
    exit()
file_q=sys.argv[1]
file_case=sys.argv[2]
file_cntl=sys.argv[3]
outf=sys.argv[4]
# print(file_q
#)
# print(file_case)
i=0
def_val='            '
kmers_map=collections.defaultdict(lambda : (def_val,def_val,def_val))
kmers_match=collections.defaultdict(lambda : ('-1.0','-1.0','-1.0'))
with open(file_q) as f:
    # content=f.read()
    for line in f:
        # print(line.split('\t')[0])
        kmers_map[line.split('\t')[0]]=(line.split('\t')[1],def_val,def_val)
print(kmers_map.__len__())
print(file_q+'_read')
with open(file_case) as f:
    # content=f.read()
    for line in f:
        if line.split('\t')[0] in kmers_map:
            kmers_map[line.split('\t')[0]]=(kmers_map[line.split('\t')[0]][0],line.split('\t')[3],def_val)
        # print(line.split('\t')[0])
        # kmers_map[line.split('\t')[0]]=line.split('\t')[2]
print(file_case+'_read')
with open(file_cntl) as f:
    # content=f.read()
    for line in f:
        if line.split('\t')[0] in kmers_map:
            kmers_map[line.split('\t')[0]]=(kmers_map[line.split('\t')[0]][0],kmers_map[line.split('\t')[0]][1],line.split('\t')[3])
        # print(line.split('\t')[0])
        # kmers_map[line.split('\t')[0]]=line.split('\t')[2]
print(file_cntl+'_read')
# for kmer, pvals in kmers_map.items():
#     i+=1
#     print(kmer, pvals)
#     if i>4:
#         break

with open(outf,"w") as f:
    f.write(file_q+'  '+file_case+'  '+file_cntl+'\n')
    for kmer, pvals in kmers_map.items():
        # if pvals[1]!=def_val or pvals[2]!=def_val:
            f.write(kmer+'   '+pvals[0]+'   '+pvals[1]+'   '+pvals[2]+'\n')
            # f.write('{0: <31}'.format())