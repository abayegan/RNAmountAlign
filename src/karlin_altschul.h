#ifndef KARLIN_ALTSCHUL_H_
#define KARLIN_ALTSCHUL_H_
#include <string.h>
#include <vector>
#include <map>
#include <limits>
#include "aux.h"

using namespace std;
typedef map<double,double>::const_iterator  mapit;

struct score_prob{
	map<double,double> scmap;
	void checkProb(){
		float sum=0.0;
		for ( mapit it = scmap.begin(); it != scmap.end(); ++it ){
			sum+=it->second;
		}
		if (!equal(1,sum)){
			reportErr("Sum of score probabilities is not 1!");
		}
	}
	double expectedValue(){
		float exp=0.0;
		float sc,pr;
		for ( mapit it = scmap.begin(); it != scmap.end(); ++it ){
			sc=it->first;pr=it->second;
			exp+= sc*pr;
		}
		return exp;	
	}
	double standard_deviation(double mu){
		float std=0.0;
		float sc,pr;
		for ( mapit it = scmap.begin(); it != scmap.end(); ++it ){
			sc=it->first;pr=it->second;
			std+= pr*(sc-mu)*(sc-mu);
		}
		return sqrt(std);
	}
	void print() const{
		float sc,pr;
		printf("#SCORE\tPROB\n");
		for ( mapit it = scmap.begin(); it != scmap.end(); ++it ){
			sc=it->first;pr=it->second;
			printf("%f\t%f\n",sc,pr);
		}
	}
};


//#define MAXIT 200;
//#define PER 2;//(PERCISION):round to PER decimal places

void karlin_altschul(const score_prob &sc_pr, double & K, double &lambda );
double computeLambda(const score_prob &sc_pr);
double computeK(const score_prob &sc_pr, double lam);
double computeH(const score_prob &sc_pr,double lam);
#endif /* KARLIN_ALTSCHUL_H_ */
