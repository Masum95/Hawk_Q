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
#include <unistd.h>
#include "kmer.h" // different constat values are defined there
#include "hashTable.hpp" 
#include "kmerClass.hpp"
#include "specialfunctions.h"
#include "const.hpp"


#define eps 1e-30

#pragma pack(push, 1)



#pragma pack(pop)




class Factorials
{
public:
	double *factorials;
	Factorials(int max);
	double getFactorial(int number);
};

Factorials *f;





// compute likelihood ratio for every kmer( i.e. how better this kmer fit the phenotype values )







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

long long int findSignificantP(long long int totkmer){
	long long int m = totkmer;

	FILE *inFile = fopen(BONF_KMERDIFF_FILE, "r");


	char *line= new char[MAX_REC_LEN];
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	char kmerString[200];

	long long int k = 1;
	long long int maxk = 0;

	//FILE *in=fopen("pvals_control_top_merged_sorted.txt","r");
	vector<double> pvalList;
	char *temp;
	double pVal;
	while (fgets(line, MAX_FILE_READ, inFile) != NULL)
	{

		temp = strtok(line, "\t\n ");
		strcpy(kmerString, temp);

		temp = strtok(NULL, "\t\n ");
		pVal = atof(temp);

		if (pVal >= 0.05)
			continue;

		pvalList.push_back(pVal);
	}

	sort(pvalList.begin(), pvalList.end()); 

	for(k=1; k<=pvalList.size(); k++){

		pVal = pvalList[k-1];
		double comp = 0.05*(((double)k)/(double)m);
		if(pVal<comp || fabs(pVal-comp)<eps){
			maxk = max(k,maxk);
		}
	}

				// fputs(line2, fileOut);

	fclose(inFile);
	return maxk;
}

int main(int argc, char * const argv[])
{	
    if (argc < 3)
    {
        fprintf(stderr, "hawk_q -n $noSamples -t $NUM_THREADS\n");
        exit(EXIT_FAILURE);
    }
    int opt;
    while ((opt = getopt(argc, argv, "n:t:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            noSamples = atoi(optarg);
            break;
        case 't':
            NUM_THREADS = atoi(optarg);
            break;
        }
    }
    if (noSamples <= 0)
    {
        fprintf(stderr, "Provide a valid $noSamples\n");
        exit(EXIT_FAILURE);
    }
    if (NUM_THREADS <= 0)
    {
        fprintf(stderr, "Provide a valid $noSamples\n");
        exit(EXIT_FAILURE);
    }
    printf("#Samples: %d\n", noSamples);
    printf("#Threads: %d\n", NUM_THREADS);

	ht = new HashTable(noSamples);

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


	FILE *totalKmerFile = fopen(TOTAL_KMERS_FILE, "w");
	fprintf(totalKmerFile, "%lld\n", ht->totalKmers);

	fclose(totalKmerFile);

	FILE *fileOut = fopen(UNIQUE_KMERDIFF_FILE, "w");

	// int count1, count2;
	double pVal;
	// // in bonferonni correcton we divide sig_level by number of test
	// // sort out those kmers which pass bonferroni correction test  and write them to UNIQUE_KMERDIFF_FILE
	
	double pvalThreshold = findSignificantP(ht->totalKmers); // Benjamini and Hochberg test 

	while (fgets(line, MAX_FILE_READ, outFile) != NULL)
	{
		strcpy(line2, line);

		temp = strtok(line, "\t\n ");
		strcpy(kmerString, temp);

		temp = strtok(NULL, "\t\n ");
		pVal = atof(temp);

		if( pVal < pvalThreshold ) // divided by number of unique kmers
		{
			fputs(line2, fileOut);
		}
	}
	
	fclose(fileOut);

	return 0;
}

