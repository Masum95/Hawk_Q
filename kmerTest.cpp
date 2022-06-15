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

#include <fstream>
#include "specialfunctions.h"

#define MAX_REC_LEN 10240

string kmerToTest;
double pval;

int noSamples;

/*-------------------------------------------------File Dependency------------------------------------------------------*/

const char *GWAS_IND_FILE = "gwas_eigenstratX.ind";
const char *TOTAL_KMERS_FILE = "gwas_eigenstratX.total";
const char *KMER_FILE = "kmers_matchWithBonf_count";

/*-------------------------------------------------File Dependency End------------------------------------------------------*/


long long int *totalKmerCounts;
double *phenotypeValues;
double likelihoodNULL;
int *counts;

void init()
{
	totalKmerCounts = new long long int[noSamples];
	phenotypeValues = new double[noSamples];
	counts = new int[noSamples];
	for (int i = 0; i < noSamples; i++)
		totalKmerCounts[i] = 0;
}
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

double logPI = 0.9189385;

void computeNullLikelihood()
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
	cout<< "y_mean " << y_mean <<" y_var "<< e_null<<endl;
	likelihoodNULL = 0;

	// for (int k = 0; k < noSamples; k++) // first 2 parts should be out of loop ?
	// {
	// 	likelihoodNULL += (-log(e_null) / 2.0 - logPI - (y[k] - y_mean) * (y[k] - y_mean) / (2 * e_null));
	// }
	double likelihoodNULL2 = (-log(e_null) / 2.0 * noSamples - logPI * noSamples - noSamples / 2.0);
	likelihoodNULL = (-log(e_null) / 2.0 * noSamples );

	cout << "Null Likelihood " << likelihoodNULL << "likelihoodNul one liner = "<< likelihoodNULL2<<endl;
}
long long int done = 0;

void likelihoodRatio()
{
	long tid;
	double likelihoodAlt, likelihoodRatio;

	double *y;

	double *x = new double[noSamples];

	double e_alt, y_p;

	int presentCount = 0;
	double presentRatio = 0;
	done = 0;

	for (int k = 0; k < noSamples; k++) // loop over sample for jth kmer
	{
		x[k] = counts[k] / (double)totalKmerCounts[k];
		// cout<< counts[k] <<" "<< (double)totalKmerCounts[k] <<" "<< x[k]<<endl;
	}
	y = phenotypeValues;

	vector<double> result = regress(noSamples, x, y);

	cout<<"m " <<result[1] <<" "<< "b "<< result[0]<<endl;
	e_alt = 0;

	for (int k = 0; k < noSamples; k++)
	{
		y_p = result[0] + result[1] * x[k];
		e_alt += (y[k] - y_p) * (y[k] - y_p);
		//  cout<<y_p<<endl;
	}
	e_alt = e_alt / noSamples;
	cout<<"e_alt " << e_alt<<endl;
	// e_alt = sqrt(e_alt / noSamples); #variance needed, not stdev

	likelihoodAlt = 0;
	//CHK
	double likelihoodAlt2 = 0;

	// for (int k = 0; k < noSamples; k++)
	// {
	// 	double y_p = y[k] - (result[0] + result[1] * x[k]);
	// 	likelihoodAlt += (-log(e_alt)/2 - logPI - (y_p * y_p ) / (2 * e_alt ));
	// }
	// likelihoodNULL = (-log(e_null) / 2.0 * noSamples - logPI * noSamples - noSamples / 2.0);
	likelihoodAlt = (-log(e_alt) / 2.0 * noSamples );

	// likelihoodAlt = (-log(e_alt) / 2.0 * noSamples - logPI * noSamples - noSamples / 2.0);
	cout<< "likelihoodAlt " <<likelihoodAlt<<endl;
	likelihoodRatio = likelihoodAlt - likelihoodNULL;
	cout<<"likelihood Ratio  "<<likelihoodRatio<<endl;

	if (likelihoodRatio < 0)
	{
		likelihoodRatio = 0;
	}
	double pVal = alglib::chisquarecdistribution(1, 2 * likelihoodRatio);
	cout<<"pval list"<<endl;
	for(int i=0;i<100;i++)
		cout<<i<<" --->"<< alglib::chisquarecdistribution(1,  i)<<endl;
	cout<<"Pval = " <<pVal<<endl;

}



int main(int argc, const char *argv[])
{

	noSamples = atoi(argv[1]);

	init();


	std::ifstream infile(KMER_FILE);

	infile >> kmerToTest;
	infile >> pval;


	string tmp, tmp2;
	int  i = 0;
	while (infile >> counts[i++]);
	infile.close();

	infile.open(GWAS_IND_FILE);
	i = 0;
	while (infile >> tmp >> tmp2 >> phenotypeValues[i++]);
	infile.close();

	infile.open(TOTAL_KMERS_FILE);
	i = 0;
	while (infile >> totalKmerCounts[i++]);

	computeNullLikelihood();
	cout<<likelihoodNULL<<endl;
	likelihoodRatio();
}
