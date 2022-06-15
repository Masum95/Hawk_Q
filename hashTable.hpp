

#ifndef HASH_TABLE_CLASS_HPP

#define HASH_TABLE_CLASS_HPP


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
#include "kmerClass.hpp"
#include "specialfunctions.h"

#include "const.hpp"
#include "utilFunctions.hpp"

long long int done = 0;





class HashTable
{
public:
	long long int totalKmers;
	long long int totalTests;
	int noSamples;
	vector<Kmer *> kmers[HASH_TABLE_LENGTH];
	pthread_mutex_t hashTableBuckets_mutex[HASH_TABLE_LENGTH];

	long long int *totalKmerCounts;
	double *phenotypeValues;
	double likelihoodNULL;

	HashTable(int);
	void insertKmer(long long int val, int count, int sampleNo);
	void show();
	void computeNullLikelihood();
	void computePvalues();
	void computeLikelihoodRatios();
	void dumpKmers(double sigLevel);
	// void *likelihoodRatio_thread(void*);

};

HashTable *ht;


HashTable::HashTable(int noSamples)
{
	this->noSamples = noSamples;
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

			e_alt = e_alt / noSamples;
			// e_alt = sqrt(e_alt / noSamples); #variance needed, not stdev

			likelihoodNull = ht->likelihoodNULL, likelihoodAlt = 0;
			//CHK
			double likelihoodAlt2 = 0;

			// for (int k = 0; k < noSamples; k++)
			// {
			// 	double y_p = y[k] - (result[0] + result[1] * x[k]);
			// 	likelihoodAlt += (-log(e_alt)/2 - logPI - (y_p * y_p ) / (2 * e_alt ));
			// }
			likelihoodAlt = (-log(e_alt) / 2.0 * noSamples - logPI * noSamples - noSamples / 2.0);

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


void *dump_thread(void *threadid)
{
	//printf("In dump Thread -> tid,  HASH_TABLE_LENGTH\n");
	long tid;
	tid = (long)threadid;

	char kmerString[100];

	eigenGenoFile = fopen(EIGENGENO_FILE, "a");
	eigenSNPFile = fopen(EIGENSNP_FILE, "a");

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
			//if (ht->kmers[i][j]->forPval == 'y')
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

	for (int k = 0; k < noSamples; k++) 
	{
		likelihoodNULL += (-log(e_null) / 2.0 - logPI - (y[k] - y_mean) * (y[k] - y_mean) / (2 * e_null));
	}
	//cout << "Null Likelihood " << likelihoodNULL << endl;
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
}

#endif

