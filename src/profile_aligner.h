/*
 * ProfileAligner.h
 *
 *  Created on: Dec 18, 2016
 *      Author: Amir Bayegan
 *
 */

#ifndef ENERGY_PARAMS_PROFILEALIGNER_H_
#define ENERGY_PARAMS_PROFILEALIGNER_H_
#include<vector>
#include<math.h>
#include "aux.h"
#include "profile.h"
using namespace std;

class ProfileAligner {
public:
	ProfileAligner(const Profile &pr1, const Profile &pr2, const Arguments & args, double shift, double mul);
	void printAlignment(string);
	void computeQuality();
	Profile align();
	char getType();
	int getLength();
	int* getStart();
	int* getStop();
	double getScore();
	ProfileAligner();
	virtual ~ProfileAligner();
private:
	char aln_type;
	Profile p1;
	Profile p2;
	Profile p; /*profile after alignment*/
	int aln_start[2];
	int aln_stop[2];
	int len_aln;
	float gap_init;
	float gap_ext;
	float gap_init_str;
	float gap_ext_str;
	float weighted_gap_init;
	float weighted_gap_ext;
	float str_weight;
	double mult_seq, shift_str;
	int cut_off;
	vector<vector<double> > scoreMatrix_p;
	vector<vector<double> > scoreMatrix_q;
	vector<vector<double> > scoreMatrix_r;
	vector<vector<char> > pathMatrix_p;
	vector<vector<char> > pathMatrix_q;
	vector<vector<char> > pathMatrix_r;
	const char * sub_mat_name;
	float seq_similarity[5][5];
	//double ** similarity;
	double max_score;
	string quality;
	float gap(int);
	bool calc_prob;

	void buildNeedlemanWunschScoreMatrix();
	void backtrackNeedlemanWunsch();
	void buildSmithWatermanScoreMatrix();
	void backtrackSmithWaterman();
	double computeMaxHeight();
	double ** buildSimilarityMatrix();
	double similarity(int i,int j);
	/*extends all sequences of the profile alignment by one character */
	void extendProfile(vector<string> &, vector<vector<double> > &, int , int , int );
	/*reverse every sequence in the profile--move this to profile class*/
	void reverseProfile();
	double getSumOfPairsScore(int ,int );

	/*keeps track of the alignment positions for every sequence in profiles(important for local alignment)*/
	vector<vector<int> > computeAlignmentPositions();
	vector<vector<int> >  mergePositions(vector<vector<int> > ,vector<vector<int> >);
	
};

#endif /* ENERGY_PARAMS_PROFILEALIGNER_H_ */
