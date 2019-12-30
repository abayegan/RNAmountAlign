#ifndef EVD_H_
#define EVD_H_


#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <algorithm>  

using namespace std;

class EVD {
public:
	EVD(vector<double>, int C);
	virtual ~EVD();
	double getStats();
	int fit();
	double getLambda(){return lambda;}
	double getMu(){return mu;}
	int getSize(){return size;}
	int getNumIter(){return num_iterations;}
	double getSum(){return sum_of_elements;}
	

private:
	vector<double> data;
	double lambda;
	double mu;
	int size;
	double sum_of_elements;
	int num_iterations;
	int CENSORED; //Ignore the left tail in fitting
	vector< vector<double> >  histogram;
	void makeHistogram();
	void printHistogram(double min, double max, double bin_size);
	int censoredFit(vector<double> data, double c, int z);
	int totalFit(vector<double> data);
	double pvalue(double x, double mu, double lambda);
};

#endif /* EVD_H_ */
