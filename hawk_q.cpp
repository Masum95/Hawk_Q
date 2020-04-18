//
//  main.cpp
//  kmer
//

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>
using namespace std;
#include <bits/stdc++.h>
#include <pthread.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "kmer.h" // different constat values are defined there

#include "specialfunctions.h"

#define MAX_REC_LEN 10240

int noSamples;

int kmerLength = 31;

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

FILE *eigenGenoFile;
FILE *eigenSNPFile;

#pragma pack(push, 1)
class Kmer
{
public:
	long long int kmer;
	unsigned short *counts;
	double pVal;

	char forPCA;
	char forPval;
	Kmer(int noSamples);
	~Kmer();
	void show();
	void freeMemory();
};
#pragma pack(pop)

class KeyVal
{
public:
	long long int val;
	int count;
};

class HashTable
{
public:
	long long int totalKmers;
	long long int totalTests;
	vector<Kmer *> kmers[HASH_TABLE_LENGTH];
	pthread_mutex_t hashTableBuckets_mutex[HASH_TABLE_LENGTH];

	long long int *totalKmerCounts;
	double *phenotypeValues;
	double likelihoodNULL;

	HashTable();
	void insertKmer(long long int val, int count, int sampleNo);
	void show();
	void computeNullLikelihood();
	void computePvalues();
	void computeLikelihoodRatios();
	void dumpKmers(double sigLevel);
};

HashTable *ht;

class Factorials
{
public:
	double *factorials;
	Factorials(int max);
	double getFactorial(int number);
};

Factorials *f;

long long int getInt(char *s)
{
	long long int val = 0;
	int i = 0;
	char ch;
	while (s[i])
	{
		val = val << 2;

		ch = s[i];
		if (ch == 'A')
		{
			val = val | 0;
		}
		else if (ch == 'C')
		{
			val = val | 1;
		}
		else if (ch == 'G')
		{
			val = val | 2;
		}
		else
		{
			val = val | 3;
		}
		i++;
	}
	return val;
}

char *getKmer(long long int val, char *kmer, int kmerLength)
{

	int temp = 0;
	for (int i = kmerLength - 1; i >= 0; i--)
	{
		temp = val & 3;
		val = val >> 2;
		if (temp == 0)
			kmer[i] = 'A';
		else if (temp == 1)
			kmer[i] = 'C';
		else if (temp == 2)
			kmer[i] = 'G';
		else if (temp == 3)
			kmer[i] = 'T';
	}
	kmer[kmerLength] = '\0';
	return kmer;
}

unsigned long int getHash(unsigned long long int key)
{
	/*
     key = (~key) + (key << 18); // key = (key << 18) - key - 1;
     key = key ^ (key >> 31);
     key = key * 21; // key = (key + (key << 2)) + (key << 4);
     key = key ^ (key >> 11);
     key = key + (key << 6);
     key = key ^ (key >> 22);
     */
	return (unsigned long int)key;
}

Kmer::Kmer(int noSamples)
{
	counts = new unsigned short[noSamples];
	for (int i = 0; i < noSamples; i++)
		counts[i] = 0;
	//	caseCounts=(char *)calloc(noCases,sizeof(char));

	pVal = 0;
}

Kmer::~Kmer()
{
	//    delete [] caseCounts;
	//    delete [] controlCounts;
}

void Kmer::freeMemory()
{
	delete[] counts;

}

void Kmer::show()
{
	char kmerString[100];

	cout << getKmer(kmer, kmerString, 31) << " ";

	for (int k = 0; k < noSamples; k++)
	{
		cout << (int)counts[k] << " ";
	}
	cout << pVal << endl;
}

HashTable::HashTable()
{
	totalKmerCounts = new long long int[noSamples];
	phenotypeValues = new double[noSamples];
	for (int i = 0; i < noSamples; i++)
		totalKmerCounts[i] = 0;

	for (int i = 0; i < HASH_TABLE_LENGTH; i++)
	{
		pthread_mutex_init(&hashTableBuckets_mutex[i], NULL);
	}

	totalKmers = 0;
	totalTests = 0;
}

void HashTable::insertKmer(long long int val, int count, int sampleNo)
{
	unsigned long int index = getHash(val) % HASH_TABLE_LENGTH;
	int found = 0;

	pthread_mutex_lock(&hashTableBuckets_mutex[index]);

	for (int i = 0; i < kmers[index].size(); i++)
	{
		if (val == kmers[index][i]->kmer)
		{

			kmers[index][i]->counts[sampleNo] = count;
			found = 1;
			break;
		}
	}
	if (found == 0)
	{
		Kmer *km = new Kmer(noSamples);
		km->kmer = val;

		km->counts[sampleNo] = count;
		kmers[index].push_back(km);

		pthread_mutex_lock(&countKmer_mutex);
		totalKmers++;
		pthread_mutex_unlock(&countKmer_mutex);
		printf("%lld\n", totalKmers);
	}
	pthread_mutex_unlock(&hashTableBuckets_mutex[index]);
}

// regresion with input : x = normalized ith kmer_count for each sample and y = phenotype_value
// return correlation coefficient (covariance / variance ) and linear regressed value
vector<double> regress(int noSamples, double *x, double *y)
{
	vector<double> result;
	double a, b;
	double x_mean = 0, y_mean = 0;
	double s_x = 0, s_y = 0, s_xy = 0;

	for (int i = 0; i < noSamples; i++)
	{
		x_mean += x[i];
		y_mean += y[i];
	}
	x_mean = x_mean / noSamples;
	y_mean = y_mean / noSamples;

	for (int i = 0; i < noSamples; i++)
	{
		s_x += (x[i] - x_mean) * (x[i] - x_mean);
		s_y += (y[i] - y_mean) * (y[i] - y_mean);
		s_xy += (x[i] - x_mean) * (y[i] - y_mean);
	}

	b = s_xy / s_x;			 // slope
	a = y_mean - b * x_mean; // y-intercept
	// regression using covariance and correlation between ith kmer count and phenotype value
	result.push_back(a);
	result.push_back(b);
	return result;
}

double normalProb(double x, double m, double s)
{
	double logpi = 0.9189385;

	return -log(s) - logpi - (x - m) * (x - m) / (2 * s * s);
}

double logPI = 0.9189385;

void HashTable::computeNullLikelihood()
{
	double *y = phenotypeValues;
	double y_mean = 0;

	for (int k = 0; k < noSamples; k++)
	{
		y_mean += y[k];
	}
	y_mean = y_mean / noSamples;

	double e_null = 0;

	for (int k = 0; k < noSamples; k++)
	{
		e_null += (y[k] - y_mean) * (y[k] - y_mean);
	}
	e_null = e_null / noSamples;
	// e_null = sqrt(e_null / noSamples); // variance of phenotype value w.r.t mean

	likelihoodNULL = 0;

	for (int k = 0; k < noSamples; k++) // first 2 parts should be out of loop ?
	{
		likelihoodNULL += (-log(e_null)/2.0 - logPI - (y[k] - y_mean) * (y[k] - y_mean) / (2 * e_null ));
	}
	//cout << "Null Likelihood " << likelihoodNULL << endl;
}
// compute likelihood ratio for every kmer( i.e. how better this kmer fit the phenotype values )
long long int done = 0;

void *likelihoodRatio_thread(void *threadid)
{
	long tid;
	tid = (long)threadid;
	double likelihoodNull, likelihoodAlt, likelihoodRatio;

	double *y;

	double *x = new double[noSamples];

	double e_alt, y_p;

	int presentCount = 0;
	double presentRatio = 0;
	//cout << "starts here" << endl;
	done = 0;
	for (int i = tid; i < HASH_TABLE_LENGTH; i += NUM_THREADS)
	{
		//cout << i << " %%% " << ht->kmers[i].size() << endl;
		for (int j = 0; j < ht->kmers[i].size(); j++) // loop over all kmers that's been hashed to ith index of hash-table
		{

			presentCount = 0;
			for (int k = 0; k < noSamples; k++) // loop over sample for jth kmer
			{
				x[k] = ht->kmers[i][j]->counts[k] / (double)ht->totalKmerCounts[k];

				if (ht->kmers[i][j]->counts[k] > 0)
				{
					presentCount++;
				}
			}
			y = ht->phenotypeValues;

			vector<double> result = regress(noSamples, x, y);

			presentRatio = presentCount / (double)(noSamples);
			if (presentRatio >= 0.01 && presentRatio <= 0.99)
			{
				ht->kmers[i][j]->forPCA = 'y';
			}
			else
			{
				ht->kmers[i][j]->forPCA = 'n';
			}
			if (presentRatio >= 0.05)
			{
				ht->kmers[i][j]->forPval = 'y';
			}
			else
			{
				ht->kmers[i][j]->forPval = 'n';
			}

			e_alt = 0;

			for (int k = 0; k < noSamples; k++)
			{
				y_p = result[0] + result[1] * x[k];
				e_alt += (y[k] - y_p) * (y[k] - y_p);
			}

			e_alt = e_alt/noSamples;
			// e_alt = sqrt(e_alt / noSamples); #variance needed, not stdev

			likelihoodNull = ht->likelihoodNULL, likelihoodAlt = 0;
			//CHK
			double likelihoodAlt2 = 0;
		
			for (int k = 0; k < noSamples; k++)
			{
				double y_p = y[k] - (result[0] + result[1] * x[k]);
				likelihoodAlt += (-log(e_alt)/2 - logPI - (y_p * y_p ) / (2 * e_alt ));
			}
			// likelihoodAlt = (-log(e_alt)/2.0 * noSamples - logPI * noSamples - noSamples/2.0 );
			// if( abs(likelihoodAlt - likelihoodAlt2) >= 0.00005){
			// 	cout<<likelihoodAlt<<" - "<<likelihoodAlt2<<endl;
			// }
			likelihoodRatio = likelihoodAlt - likelihoodNull;

			if (likelihoodRatio < 0)
			{
				likelihoodRatio = 0;
			}

			double pVal = alglib::chisquarecdistribution(1, 2 * likelihoodRatio);

			ht->kmers[i][j]->pVal = pVal;
			done++;
		}
	}

	delete[] x;

	pthread_exit(NULL);
}

void HashTable::computeLikelihoodRatios()
{

	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for (t = 0; t < NUM_THREADS; t++)
	{
		rc = pthread_create(&threads[t], NULL, likelihoodRatio_thread, (void *)t);
		if (rc)
		{
			exit(-1);
		}
	}

	for (t = 0; t < NUM_THREADS; t++)
	{
		rc = pthread_join(threads[t], &status);
		if (rc)
		{
			exit(-1);
		}
	}
}

FILE *outFile;
double pValThreshold;

void *dump_thread(void *threadid)
{
	//printf("In dump Thread -> tid,  HASH_TABLE_LENGTH\n");
	long tid;
	tid = (long)threadid;

	char kmerString[100];

	for (int i = tid; i < HASH_TABLE_LENGTH; i += NUM_THREADS)
	{
		for (int j = 0; j < ht->kmers[i].size(); j++)
		{
			//cout<<ht->kmers[i][j]->forPCA<<"--"<<ht->kmers[i][j]->pVal<<"**"<<pValThreshold<<endl;
			// This part was missing which generates SNP and genoType file for eigenstrat. Not sure whether this is necessary.
			if (ht->kmers[i][j]->forPCA == 'y' && rand() / ((double)(RAND_MAX) + 1) < 0.01)
			{
				pthread_mutex_lock(&eigenFile_mutex);

				fprintf(eigenSNPFile, "%s\t%d\t%lf\t%d\n", getKmer(ht->kmers[i][j]->kmer, kmerString, 31), 1, 0.0, 0);

				for (int k = 0; k < noSamples; k++)
				{
					fprintf(eigenGenoFile, "%d\t", ht->kmers[i][j]->counts[k] > 0 ? 1 : 0);
				}

				fprintf(eigenGenoFile, "\n");
				pthread_mutex_unlock(&eigenFile_mutex);
			}

			// print all kmer to BONF_KMERDIFF_FILE for which we reject null hypothesis. That is significat kmer
			if (ht->kmers[i][j]->pVal <= pValThreshold && ht->kmers[i][j]->forPval == 'y')
			{
				pthread_mutex_lock(&outFile_mutex);

				fprintf(outFile, "%s\t%e\t", getKmer(ht->kmers[i][j]->kmer, kmerString, 31), ht->kmers[i][j]->pVal);

				for (int k = 0; k < noSamples; k++)
				{
					fprintf(outFile, "%d\t", ht->kmers[i][j]->counts[k]);
				}
				fprintf(outFile, "\n");
				pthread_mutex_unlock(&outFile_mutex);
			}
		}
	}
	// printf("Progress %0.2lf%\n", done * 1.0 / ht->totalKmers);

	for (int i = tid; i < HASH_TABLE_LENGTH; i += NUM_THREADS)
	{
		for (int j = 0; j < ht->kmers[i].size(); j++)
		{
			ht->kmers[i][j]->freeMemory();
			delete ht->kmers[i][j];
		}
		ht->kmers[i].clear();
	}

	pthread_exit(NULL);
}

void HashTable::dumpKmers(double sigLevel)
{
	outFile = fopen(BONF_KMERDIFF_FILE, "a");

	eigenGenoFile = fopen(EIGENGENO_FILE, "a");
	eigenSNPFile = fopen(EIGENSNP_FILE, "a");
	//CHECK
	//pValThreshold=sigLevel/(double)CUTOFF;
	pValThreshold = sigLevel / 100;

	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for (t = 0; t < NUM_THREADS; t++)
	{
		rc = pthread_create(&threads[t], NULL, dump_thread, (void *)t);
		if (rc)
		{
			exit(-1);
		}
	}

	for (t = 0; t < NUM_THREADS; t++)
	{
		rc = pthread_join(threads[t], &status);
		if (rc)
		{
			exit(-1);
		}
	}

	fclose(outFile);
	fclose(eigenSNPFile);
	fclose(eigenGenoFile);
}

void HashTable::show()
{
	/*    
    for(int i=0;i<HASH_TABLE_LENGTH;i++)
    {
        for(int j=0;j<kmers[i].size();j++)
        {
            kmers[i][j]->show();
        }
        
    }
*/
	// cout << totalKmers << endl;
	// cout << totalTests << endl;
}

void getKeyVal(char *s, KeyVal *kv)
{
	long long int val = 0;
	int countVal = 0;
	int i = 0;
	char ch;
	while (i < KMER_LENGTH)
	{

		ch = s[i++];

		val = val << 2 | bases[ch];
	}
	i++;
	kv->val = val;
	while (1)
	{
		ch = s[i];
		if (ch == '\0' || ch == '\n')
		{
			break;
		}
		countVal = countVal * 10 + ch - '0';
		i++;
	}
	kv->count = countVal;
}

FILE **kmerFiles;
long long int *vals;
unsigned short int *counts;

struct ThreadArg
{
	long long int valBar;
	int threadID;
};

void *readSamples(void *threadid)
{
	ThreadArg *ta = (ThreadArg *)threadid;
	long long int valBar = ta->valBar;
	int threadNo = ta->threadID;
	long long int val;
	int count;

	char *line = new char[MAX_REC_LEN];
	int MAX_FILE_READ = MAX_REC_LEN / sizeof(line[0]);
	for (int i = threadNo; i < noSamples; i += NUM_THREADS)
	{
		if (vals[i] != -1 && vals[i] < valBar)
		{
			ht->insertKmer(vals[i], counts[i], i);
			vals[i] = -1;
			counts[i] = -1;
		}
		if (vals[i] == -1)
		{
			while (fscanf(kmerFiles[i], "%lld %d\n", &val, &count) != EOF)
			{

				if (val < valBar)
				{
					ht->insertKmer(val, count, i);
				}
				else
				{
					vals[i] = val;
					counts[i] = count;
					break;
				}
			}
		}
	}

	delete[] line;
	pthread_exit(NULL);
}

int main(int argc, const char *argv[])
{

	noSamples = atoi(argv[1]);

	ht = new HashTable();

	FILE *countsFile = fopen(TOTAL_KMER_CNT_FILE, "r");
	FILE *phenotypeFile = fopen(PHENOVALUE_FILE, "r");

	for (int i = 0; i < noSamples; i++)
	{
		fscanf(countsFile, "%lld\n", &ht->totalKmerCounts[i]);
	}
	for (int i = 0; i < noSamples; i++)
	{
		fscanf(phenotypeFile, "%lf\n", &ht->phenotypeValues[i]);
	}
	fclose(countsFile);
	fclose(phenotypeFile);

	ht->computeNullLikelihood();

	FILE *outFile = fopen(BONF_KMERDIFF_FILE, "w");

	fclose(outFile);
	FILE *eigenGenoFile = fopen(EIGENGENO_FILE, "w");
	FILE *eigenSNPFile = fopen(EIGENSNP_FILE, "w");

	fclose(eigenSNPFile);
	fclose(eigenGenoFile);

	kmerFiles = new FILE *[noSamples];
	char *kmerFilename;
	kmerFilename = new char[5000];

	vals = new long long int[noSamples];
	counts = new unsigned short int[noSamples];
	char *line = new char[MAX_REC_LEN];
	char *line2 = new char[MAX_REC_LEN];
	int MAX_FILE_READ = MAX_REC_LEN / sizeof(line[0]);

	char *temp;
	char kmerString[100];
	unsigned short count;

	long long int valMax = 0x3FFFFFFFFFFFFFFF;
	long long int valBar = 0;
	long long int val = 0;
	long long int valInc = 0x0010000000000000;
	//	long long int valInc=0x0000000000100000;

	FILE *sortedFile = fopen(SORTED_FILE, "r");

	for (int i = 0; i < noSamples; i++)
	{
		fscanf(sortedFile, "%s\n", kmerFilename);
		//	sprintf(kmerFilename,"%d_case_kmers_sorted.txt",(i+1));
		kmerFiles[i] = fopen(kmerFilename, "r");

		if (kmerFiles[i] == NULL)
		{
			cout << kmerFilename << " file doesn't exist" << endl;
		}

		vals[i] = -1;
		counts[i] = -1;
	}
	ThreadArg *thArgs[NUM_THREADS];

	while (valBar < valMax)
	{
		valBar += valInc;

		pthread_t threads[NUM_THREADS];
		int rc;
		void *status;
		for (int i = 0; i < NUM_THREADS; i++)
		{
			thArgs[i] = new ThreadArg;
			thArgs[i]->valBar = valBar;
			thArgs[i]->threadID = i;
			rc = pthread_create(&threads[i], NULL, readSamples, (void *)thArgs[i]);
			if (rc)
			{
				exit(-1);
			}
		}

		for (int i = 0; i < NUM_THREADS; i++)
		{
			rc = pthread_join(threads[i], &status);
			delete thArgs[i];
			if (rc)
			{
				exit(-1);
			}
		}
		ht->computeLikelihoodRatios();
		ht->dumpKmers(SIG_LEVEL);
	}
	ht->show();

	outFile = fopen(BONF_KMERDIFF_FILE, "r");

	FILE *fileOut = fopen(UNIQUE_KMERDIFF_FILE, "w");

	FILE *totalKmerFile = fopen(TOTAL_KMERS_FILE, "w");
	fprintf(totalKmerFile, "%lld\n", ht->totalKmers);

	int count1, count2;
	double pVal;
	// in bonferonni correcton we divide sig_level by number of test
	// sort out those kmers which pass bonferroni correction test  and write them to UNIQUE_KMERDIFF_FILE
	while (fgets(line, MAX_FILE_READ, outFile) != NULL)
	{
		strcpy(line2, line);

		temp = strtok(line, "\t\n ");
		strcpy(kmerString, temp);

		temp = strtok(NULL, "\t\n ");
		pVal = atof(temp);

		if (pVal < SIG_LEVEL / ht->totalKmers) // divided by number of unique kmers
		{
			fputs(line2, fileOut);
		}
	}

	return 0;
}
