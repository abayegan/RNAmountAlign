#ifndef STATS_H_
#define STATS_H_
#include <string.h>
#include <math.h>

#include "aux.h"
#include "mountain_height.h"
#include "arguments.h"
#include "pair_align.h"
#include "evd.h"
#include "karlin_altschul.h"


void absDiffBpProbability(const double &bp1, const double &bp2, score_prob &sc_pr);
void nucMutationProb(const vector<double> &npr1, const vector<double> &npr2 , score_prob &sc_pr, const float mat[4][4]);
void karlinStats(const vector<double> &npr1,double **bppr1, int n1, const vector<double> &npr2, double **bppr2, int n2, const float subMat[4][4],float gamma,double &K, double &lambda);
void computeStats_bin_gc( Arguments & args, vector<float> & gc_vec, vector<vector<double> > & gc_mu_lam_vec,const char dist);
double computeEvalue_EVD(double x, double lambda, double mu);
double computeEvalue(double x, double lambda, double K, int N);
double computeK(double lambda, double mu,int n);
double comuteMu(double lambda, double K, int N);
void get_params_for_gc(float gc, vector<vector<double> > & gc_mu_lam_vec, double & lambda , double & mu );
void printStats(int RND_NUM, int randlen,vector<vector<double> > gc_mu_lam_vec);
double computeVectorMean(vector<double> vec);
double computeVectorStd(vector<double> vec, double mu);
double computePvalueNormal(double x, double sigma, double mu);
double cdf_normal(double x);
double Eval2Pval(double e);
double Pval2Eval(double p);
double KAevalueRegressionTransform(double KA);
#endif /* STATS_H_ */
