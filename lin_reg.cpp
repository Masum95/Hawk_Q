#include <cstring>
#include <pthread.h>
#include <semaphore.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
#include "lr.h"
#include "specialfunctions.h"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

#define eps 1e-30
pthread_mutex_t printMutex = PTHREAD_MUTEX_INITIALIZER;

#define PCA_COUNT 10
//#define NULL_MODEL_FEATURE_COUNT 3
//#define ALT_MODEL_FEATURE_COUNT 4
#define MAX_LINE_LENGTH 512
#define CLASS_NAME_LENGTH 8
#define CLASS_NAME1 "Case"
#define CLASS_NAME2 "Control"
#define CHUNK_SIZE 2000
// #define DEBUG

using namespace std;
double ln2pi = 1.837877;
// variables, semaphores and mutex used to sync and control thread life
pthread_mutex_t done_count_lock;
int done_count;
int thread_exit_signal;
sem_t all_start;
sem_t all_done;

ifstream con_file;
ifstream feature_z_file;
ifstream ind_file;
ifstream total_file;
ifstream cov_file;

/*-------------------------------------------------File Dependency------------------------------------------------------*/

const char *EIGEN_IND_FILE = "gwas_eigenstratX.ind";
const char *BONF_KMERDIFF_TOP_FILE = "out_unique.kmerDiff";
const char *PCS_EVEC_FILE = "pcs.evec";
const char *TOTAL_KMER_CNT_FILE = "gwas_eigenstratX.total";

/*-------------------------------------------------File Dependency End------------------------------------------------------*/

/*------------------------------------------------Helper Function---------------------------------------------------------*/
vector<string> stringTokenize(string str, string delim)
{
	//cout<<str<<"*"<<endl;
	vector<string> wordVector;
	std::stringstream ss(str);
	std::string line;
	while (std::getline(ss, line))
	{
		std::size_t prev = 0, pos;
		while ((pos = line.find_first_of(delim, prev)) != std::string::npos)
		{
			if (pos > prev)
				wordVector.push_back(line.substr(prev, pos - prev));
			prev = pos + 1;
		}
		if (prev < line.length())
			wordVector.push_back(line.substr(prev, std::string::npos));
	}
	//cout<<wordVector[0]<<"-"<<wordVector[1]<<"-"<<wordVector[2]<<endl;
	ss.str("");
	return wordVector;
}
/*------------------------------------------------Helper Function---------------------------------------------------------*/

// These variables are global to ease the passing to multiple threads
int start_indx;
int num_of_thread;
int PC;
int cov_count;
string covfile;
int read_row_count;

int NULL_MODEL_FEATURE_COUNT;
int ALT_MODEL_FEATURE_COUNT;

int mx_iter;
double learn_rate;

double null_likelihood;
unsigned int nrow;
// std::vector<std::vector<double> > Z;
MatrixXd Z;
// std::vector<std::vector<double> > C;
MatrixXd C;
// std::vector<double> Y;
MatrixXd Y;
std::vector<unsigned long long int> totals;
std::vector<std::vector<unsigned long long int>> kmercounts;
std::vector<double> output;
// std::vector<std::vector<double> > global_features_NULL;
MatrixXd global_features_NULL;
// std::vector<std::vector<double> > global_features_ALT;
MatrixXd global_features_ALT;

VectorXd null_model;

struct thread_info
{
	int thread_no;
};

int open_file_connection();
int find_row_count();
void init_sync_primitve();
void *worker_thread_func(void *thread_no_as_ptr);

void printmodel(vector<double> model)
{
	for (size_t i = 0; i < model.size(); ++i)
	{
		cout << model[i] << " ";
	}
	cout << "\n";
}

// run with [-t number of worker] [-c covariate file name] [-p number of principle components]
int main(int argc, char **argv)
{
	num_of_thread = 1;
	PC = 2;
	covfile = "";
	cov_count = 0;

	for (int i = 0; i < argc; ++i)
	{
		if (strcmp(argv[i], "-t") == 0)
		{
			num_of_thread = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-c") == 0)
		{
			int ll = strlen(argv[i + 1]);
			for (int j = 0; j < ll; ++j)
			{
				covfile.push_back(argv[i + 1][j]);
			}
		}
		else if (strcmp(argv[i], "-p") == 0)
		{
			PC = atoi(argv[i + 1]);
		}
	}
	// cout << "done with args" << endl;

	if (open_file_connection())
	{
		cout << "Error in opening file" << std::endl;
		return 0;
	}

	learn_rate = 0.1;
	mx_iter = 25;

	nrow = find_row_count(); // = sample size

	//nrow and nrow are equal
	// Z = std::vector<std::vector<double> >(nrow, std::vector<double>(PCA_COUNT, 0)); //2d matrix with row equals to sample_size and column = 	10 (PCAs)							  			//PCA_cnt( all initialized to 0)
	Z.resize(nrow, PCA_COUNT);
	// Y = std::vector<double>(nrow); //phenotype_value (now case or control //CHK)
	Y.resize(nrow, 1);
	std::vector<double> RAWC;

	totals = std::vector<unsigned long long int>(nrow); // total kmers for each sample

	for (unsigned int l = 0; l < nrow; l++)
	{
		char buf[MAX_LINE_LENGTH];
		char class_name[CLASS_NAME_LENGTH];
		for (int l1 = 0; l1 < PCA_COUNT; l1++)
		{
			feature_z_file >> Z(l, l1); // pcs.evec file ( all samples components along new pcs)
		}

		/*
         * because each line's 3rd string (which is last) is of importance
         * just take three string input from the line and keep the last
         */

		ind_file.getline(buf, MAX_LINE_LENGTH - 1);

		Y(l, 0) = atof(stringTokenize(string(buf), "\t")[2].c_str());
		// cout<<Y[l]<<endl;

		total_file >> totals[l];
		// cout<<totals[l]<<endl;
	}
	//release the file connections
	feature_z_file.close();
	ind_file.close();
	total_file.close();

	//reading covariate file...
	//like Z, i dunno how much PC is there
	if ((int)covfile.size() > 0)
	{
		double cc;
		while (cov_file >> cc)
		{
			RAWC.push_back(cc);
		}
		int sz = (int)RAWC.size();
		// C = std::vector<std::vector<double> >(nrow, std::vector<double>(sz / nrow, 0));
		C.resize(nrow, sz / nrow);
		int k = 0;
		for (int i = 0; i < nrow; ++i)
		{
			for (int j = 0; j < (sz / nrow); ++j)
			{
				C(i, j) = RAWC[k];
				k++;
			}
		}
		cov_count = sz / nrow;
		cov_file.close();
	}

#ifdef DEBUG
	cout << "nrow : " << nrow << std::endl;
	cout << "Z" << std::endl;
	for (int l = 0; l < nrow; l++)
	{
		for (int l1 = 0; l1 < PCA_COUNT; l1++)
		{
			cout << Z(l, l1) << ' ';
		}
		cout << std::endl;
	}

	cout << "Y" << std::endl;
	for (int l = 0; l < nrow; l++)
	{
		cout << Y[l] << ' ';
	}
	cout << std::endl;

	cout << "totals" << std::endl;
	for (int l = 0; l < nrow; l++)
	{
		cout << totals[l] << ' ';
	}
	cout << std::endl;
#endif

	/*
	 * below matrix creation is done for fitting glm using glm function
	 * 4th column of matrix will be different for each sample (as per understanding)
	 */
	int chunk_size = CHUNK_SIZE;

	NULL_MODEL_FEATURE_COUNT = 1 + PC + cov_count + 1;		// 1 for total_kmer_cnt for ith sample
	ALT_MODEL_FEATURE_COUNT = 1 + NULL_MODEL_FEATURE_COUNT; // 1 extra for normalized kmer counts

	// global_features_NULL = std::vector<std::vector<double> >(nrow, std::vector<double>(NULL_MODEL_FEATURE_COUNT));
	// global_features_ALT = std::vector<std::vector<double> >(nrow, std::vector<double>(ALT_MODEL_FEATURE_COUNT));

	global_features_NULL.resize(nrow, NULL_MODEL_FEATURE_COUNT);
	global_features_ALT.resize(nrow, ALT_MODEL_FEATURE_COUNT);
	for (unsigned int l = 0; l < nrow; l++)
	{
		global_features_NULL(l, 0) = 1;
		global_features_ALT(l, 0) = 1;
		for (unsigned int z = 0; z < PC; ++z)
		{
			global_features_NULL(l, z) = Z(l, z);
			global_features_ALT(l, z) = Z(l, z);
		}
		for (unsigned int c = 0; c < cov_count; ++c)
		{
			global_features_NULL(l, PC + c) = C(l, c);
			global_features_ALT(l, PC + c) = C(l, c);
		}
		global_features_NULL(l, PC + cov_count) = totals[l];
		global_features_ALT(l, PC + cov_count) = totals[l];
	}
	// cout << "****" << endl;
	// cout<< global_features_NULL <<endl;
	null_model = linear_regression(global_features_NULL, Y);

	null_likelihood = 0;

	MatrixXd y_hat = global_features_NULL * null_model; //X * B' // N X f * f X 1 = N X 1
	MatrixXd y_dif_y_hat = Y - y_hat;					// N X 1 								  // 4x1 matrix
	MatrixXd cov = covariance_matrix(y_dif_y_hat);

	double determinant = cov.determinant();
	// cout<<determinant<<endl;

	// cout<<y_hat.rows() <<" " <<y_hat.cols()<<endl;
	// cout<<y_dif_y_hat.rows() <<" " <<y_dif_y_hat.cols()<<endl;
	// cout<<cov.rows()<<"  "<<cov.cols()<<endl;
	double variance = (y_dif_y_hat.transpose() * y_dif_y_hat)(0, 0) / nrow;

	null_likelihood = -(ln2pi + log(variance) + 1) * nrow / 2.0;
	// cout<<null_likelihood<<"--"<<variance<< endl;
	output = std::vector<double>(CHUNK_SIZE);
	int chunkread = 0;

	/*
     * First we will init sync primitive to create a signal for threads 
     * to start
     */
	init_sync_primitve();
	// Then create the thread(s)
	std::vector<pthread_t> thread_list(num_of_thread);
	for (int l = 0; l < num_of_thread; l++)
	{
		thread_info *info = new thread_info;
		info->thread_no = l;
		pthread_create(&thread_list[l], NULL, worker_thread_func, (void *)info);
	}

	while (true)
	{
		char buf[MAX_LINE_LENGTH];
		kmercounts.clear();
		for (read_row_count = 0; read_row_count < chunk_size; read_row_count++)
		{
			con_file >> buf >> buf;

			if (con_file.eof())
			{
				break;
			}

			kmercounts.push_back(std::vector<unsigned long long int>(nrow));
			for (unsigned int l = 0; l < nrow; l++)
			{
				con_file >> kmercounts[read_row_count][l];
			}
		}

		//kmercounts er size protibar CHUNK_SIZE kore bartese
		//kmercounts er each row te 15 ta column

#ifdef DEBUG
		/*
         * loop to see extraction from con_file is done correctly
         */
		cout << "Portion of kmercounts : " << std::endl;
		for (int l = 0; l < chunk_size; l++)
		{
			for (int l1 = 0; l1 < nrow; l1++)
			{
				cout << kmercounts(l, l1) << ' ';
			}
			cout << std::endl;
		}
#endif
		//read_row_count is equal to CHUNK_SIZE
		//shudhu sesh bar ektu kom hoite pare
		/*
         * Below for loop can be done in parallel. The plan is to divide the  
         * available iteration to multiple thread. If one thread is complete 
         * it will increase a signal variable from sequence of signal and go to sleep. 
         * If all signals are marked than main thread will write the result in stdout then  
         * read more data and signal the worker thread to restart. Main thread will act 
         * as a watcher.
         * 
         * Every thread will do interleaved reading from kmercounts and write to 
         * particular loaction exclusive to thread. As it's write operation is not in 
         * same memory address for different thread no synchrnization needed
         */

		pthread_mutex_lock(&done_count_lock);
		/* 
         * done_count holds how many thread has completed their part
         * done_count is 0. So that, last thread to complete this iteration 
         * can read the done_count and signal main thread appropriately
		 */
		done_count = 0;
		pthread_mutex_unlock(&done_count_lock);

		sem_init(&all_done, 0, 0);
		// sem_init(&all_start, 0, num_of_thread-1);
		for (int l = 0; l < num_of_thread; l++)
			sem_post(&all_start);

		sem_wait(&all_done);
		// all_done semphore is posted, means all outputs are prepared.
		// So, dump them to stdout

		// for (int l = 0; l < read_row_count; l++)
		// {
		// 	cout << output[l] << endl;
		// }
		//cout<<"write"<<endl;
		cout.flush();
		//cout<<"flush"<<endl;

		if (read_row_count < chunk_size)
		{
			thread_exit_signal = 1;
			for (int l = 0; l < num_of_thread; l++)
				sem_post(&all_start);
			break;
		}
		start_indx += read_row_count;
	}

	return 0;
}

int cntpr = 0;
int open_file_connection()
{
	con_file.open(BONF_KMERDIFF_TOP_FILE);
	feature_z_file.open(PCS_EVEC_FILE);
	ind_file.open(EIGEN_IND_FILE);
	total_file.open(TOTAL_KMER_CNT_FILE);
	if ((int)covfile.size() > 0)
	{
		char cvv[200];
		for (int i = 0; i < (int)covfile.size(); ++i)
		{
			cvv[i] = covfile[i];
			cvv[i + 1] = '\0';
		}
		cov_file.open(cvv);
		if (!cov_file)
		{
			cout << covfile << " not found";
			return 1;
		}
	}

	if (!con_file)
	{
		cout << "out_wo_bonf_top.kmerdiff not found";
		return 1;
	}
	if (!feature_z_file)
	{
		cout << "pcs.evec not found";
		return 1;
	}
	if (!ind_file)
	{
		cout << "gwas_eigenstratX.ind not found";
		return 1;
	}
	if (!total_file)
	{
		cout << "total_kmer_counts.txt not found";
		return 1;
	}
	return 0;
}

// find the number of sample
int find_row_count()
{
	int nrow = 0;
	char buf[MAX_LINE_LENGTH];

	total_file.getline(buf, MAX_LINE_LENGTH - 1);
	while (!total_file.eof())
	{
		nrow++;
		total_file.getline(buf, MAX_LINE_LENGTH - 1);
	}
	total_file.clear(); // needed before seekg if not c++11
	total_file.seekg(0, ios::beg);

	return nrow;
}

void init_sync_primitve()
{
	sem_init(&all_start, 0, 0);
	sem_init(&all_done, 0, 0);
	pthread_mutex_init(&done_count_lock, 0);
	thread_exit_signal = 0;
}

void *worker_thread_func(void *arg)
{
	int thread_no = ((thread_info *)arg)->thread_no;
	int interleave = num_of_thread;

	std::vector<double> counts(nrow); // normalized count for each k_mer for each sample
	// std::vector<std::vector<double> > thread_local_features_ALT(global_features_ALT);
	MatrixXd thread_local_features_ALT = global_features_ALT;
	while (true)
	{
		sem_wait(&all_start);
		if (thread_exit_signal == 1)
		{
			break;
		}

		for (int l = thread_no; l < read_row_count; l += interleave) // loops over kmers
		{
			for (unsigned int l1 = 0; l1 < nrow; l1++) //loops over samples
			{
				counts[l1] = kmercounts[l][l1] / (double)totals[l1]; // ith kmer of jth individual is normalized
																	 // pthread_mutex_lock(&printMutex);

				// cout << kmercounts[l][l1] << " ";
				// pthread_mutex_unlock(&printMutex);
			}

			//create the fourth column of matrix
			for (unsigned int l1 = 0; l1 < nrow; l1++)
			{
				thread_local_features_ALT(l1, ALT_MODEL_FEATURE_COUNT - 1) = counts[l1]; // make ith kmer as one of feature
																						 // others are pc and covariates
			}

			// for multivariate lin_reg Y = B * X

			MatrixXd B = linear_regression(thread_local_features_ALT, Y); // get coefficients of multivariate linear regression
																		  // B =  feature X 1 matrix
			// cout<<B<<endl;
			// MatrixXd cov = covariance_matrix(thread_local_features_ALT);
			// pthread_mutex_lock(&printMutex);
			// cout << endl;
			// cout << "---------" << endl;

			// cout << cov << endl;
			// cout << "+++++++++" << endl;
			// cout << thread_local_features_ALT << endl;
			// pthread_mutex_unlock(&printMutex);
			MatrixXd y_hat = thread_local_features_ALT * B; //X * B' // N X f * f X 1 = N X 1
			MatrixXd y_dif_y_hat = Y - y_hat;				// N X 1 								  // 4x1 matrix
			// double determinant = cov.determinant();
			// cout<<determinant<<"()()"<<endl;
			// double determinant = 1;
			// log-likelihood of alternate model
			double alt_likelihood = 0;
			double variance = (y_dif_y_hat.transpose() * y_dif_y_hat)(0, 0) / nrow;
			alt_likelihood = -(ln2pi + log(variance) + 1) * nrow / 2.0;
			// likelihoodAlt = -(ln2pi  + log(variance) + 1  ) * nrow/2.0;

			// cout<<"alt_likelihood"<<alt_likelihood<<" "<<variance<< endl;
			// for (int dat = 0; dat < thread_local_features_ALT.rows(); ++dat)
			// {														// loop over sample
			// 	MatrixXd data(thread_local_features_ALT.cols(), 1); //4x1
			// 	for (int j = 0; j < thread_local_features_ALT.cols(); ++j)
			// 	{ // size = 4/5
			// 		data(j, 0) = thread_local_features_ALT(dat, j);
			// 	}
			// }
			// pthread_mutex_lock(&printMutex);
			// cout << endl;
			// cout << "---------" << endl;
			// cout<<alt_likelihood<<endl;

			// pthread_mutex_unlock(&printMutex);
			// cov = cov.block(0, 0, NULL_MODEL_FEATURE_COUNT, NULL_MODEL_FEATURE_COUNT);
			// 	pthread_mutex_lock(&printMutex);
			// cout << "---------" << endl;
			// cout<<cov<<endl;

			// pthread_mutex_unlock(&printMutex);

			// for (int dat = 0; dat < global_features_NULL.rows(); ++dat)
			// { // loop over sample

			// 	MatrixXd data(global_features_NULL.cols(), 1); //3x1
			// 	for (int j = 0; j < global_features_NULL.cols(); ++j)
			// 	{ // size = 4/5
			// 		data(j, 0) = global_features_NULL(dat, j);
			// 	}
			// 	MatrixXd mu = global_features_NULL.colwise().mean().transpose(); //3x1
			// 	MatrixXd X_mu = data - mu;										 // 3x1 matrix
			// 	double determinant = cov.determinant();
			// 	// log-likelihood of alternate model
			// 	null_likelihood += -(NULL_MODEL_FEATURE_COUNT * ln2pi + log(determinant) + (X_mu.transpose() * (cov.inverse() * X_mu)).determinant()) / 2.0;
			// }
			double log_likelihood_ratio = alt_likelihood - null_likelihood;
			output[l] = alglib::chisquarecdistribution(1, 2 * log_likelihood_ratio); // lth kmer significance

			pthread_mutex_lock(&printMutex);
			// cout << alt_likelihood << "  " << null_likelihood<<endl;

			cout << output[l] << endl;
			pthread_mutex_unlock(&printMutex);
		}

		pthread_mutex_lock(&done_count_lock);
		done_count++;
		if (done_count == num_of_thread)
		{
			sem_post(&all_done);
		}
		pthread_mutex_unlock(&done_count_lock);
	}

	// delete ((thread_info *)arg);
}
