

#ifndef KMER_CLASS_HPP

#define KMER_CLASS_HPP


#include "const.hpp"
#include "utilFunctions.hpp"


class Kmer
{
public:
	long long int kmer;
	unsigned short *counts;
    int noSamples;
	double pVal;

	char forPCA;
	char forPval;
	Kmer(int noSamples);
	~Kmer();
	void show();
	void freeMemory();
};

Kmer::Kmer(int noSamples)
{
	counts = new unsigned short[noSamples];
    this->noSamples = noSamples;
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

#endif