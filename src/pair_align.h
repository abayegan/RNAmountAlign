#include <vector>
#include <string.h>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <stdio.h>
#include "aux.h"
#include "arguments.h"

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

using namespace std;

const int INF = numeric_limits<int>::max();

class Alignment {

public:
	Alignment();
	Alignment(string,string,const vector<double> &, const vector<double> &,
				double **, double **,const Arguments & args);
	Alignment(string,string,const vector<double> &, const vector<double> &,double , double,const Arguments & args);
	void printAlignment();
	string sprintAlignment();
	void computeQuality();
	void align();
	virtual ~Alignment();

	char getType();
	string getSeq1();
	string getSeq2();
	vector<int> getAlignmentPositions();
	int getLength();
	bool operator < (const Alignment &) const;
	bool operator > (const Alignment &) const;
	double getScore();
	
	//Alignment& operator=(const Alignment& that);
private:
	char aln_type;
	string seq1_raw;
	string seq2_raw;
	int len1_raw;
	int len2_raw;
	vector<double> h1_raw;
	vector<double> h2_raw;
	string id1,id2;
	
	string seq1_aln;
	string seq2_aln;
	int aln_start[2];
	int aln_stop[2];
	vector<double> h1_aln;
	vector<double> h2_aln;
	double ** bppr1;
	double ** bppr2;
	
	int len_aln;
	vector<double> prob_profile;
	float gap_init;
	float gap_ext;
	float gap_init_str;
	float gap_ext_str;
	float weighted_gap_init;
	float weighted_gap_ext;
	float str_weight;
	float shift_str;
	float mult_seq;
	
	vector<vector<double> > scoreMatrix_p;
	vector<vector<double> > scoreMatrix_q;
	vector<vector<double> > scoreMatrix_r;
	
	vector<vector<char> > pathMatrix_p;
	vector<vector<char> > pathMatrix_q;
	vector<vector<char> > pathMatrix_r;
	
	const char * sub_mat_name;
	float seq_similarity[4][4];
	double max_score;
	string quality;
	float gap(int);
	string outformat;
	bool fixed_gamma;
	double bp1,bp2; //base pairing probability in seq1,seq2


	void buildNeedlemanWunschScoreMatrix();
	void backtrackNeedlemanWunsch();
	void buildSmithWatermanScoreMatrix();
	void backtrackSmithWaterman();
	void buildSemiGlobalScoreMatrix();
	void backtrackSemiGlobal();
	void scaleScores();
	void computeMeanStdSeqScores(vector<double> p,vector<double> q, double & mu, double & std);
	void computeMeanStdStrScores(vector<double> p,vector<double> q, double & mu, double & std);
	double similarity(int,int,bool);
	void getBPprobs();
};

#endif /* ALIGNMENT_H_ */
