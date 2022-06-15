
#ifndef UTIL_FUNCTION_H

#define UTIL_FUNCTION_H

int 
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



template <typename _T>
string NumberToString(_T Number)
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}



double normalProb(double x, double m, double s)
{
	double logpi = 0.9189385;

	return -log(s) - logpi - (x - m) * (x - m) / (2 * s * s);
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


class KeyVal
{
public:
	long long int val;
	int count;
};

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






#endif