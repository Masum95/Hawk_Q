#ifndef CONST_H

#define CONST_H



#define MAX_REC_LEN 10240

int noSamples;

int kmerLength = 31;
double logPI = 0.9189385;


int NUM_THREADS = 16; // every ith thread is responsible for calculating (i+num_thread) kmer's value
/*-------------------------------------------------File Dependency------------------------------------------------------*/

const char *EIGENGENO_FILE = "gwas_eigenstratX.geno";
const char *EIGENSNP_FILE = "gwas_eigenstratX.snp";
const char *BONF_KMERDIFF_FILE = "out_wo_bonf.kmerDiff";
const char *UNIQUE_KMERDIFF_FILE = "out_unique.kmerDiff";
const char *TOTAL_KMERS_FILE = "total_kmers.txt";
const char *SORTED_FILE = "sorted_files.txt";
const char *PHENOVALUE_FILE = "phenotype_values.txt";
const char *TOTAL_KMER_CNT_FILE = "total_kmer_counts.txt";

/*-------------------------------------------------File Dependency End------------------------------------------------------*/

pthread_mutex_t readFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t outFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t hashTable_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t eigenFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t countKmer_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t filewrite_mutex = PTHREAD_MUTEX_INITIALIZER;


FILE *outFile;
double pValThreshold;

FILE *eigenGenoFile;
FILE *eigenSNPFile;
#endif