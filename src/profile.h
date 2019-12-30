/*
 * Profile.h
 *
 *  Created on: Dec 18, 2016
 *      Author: Amir Bayegan
 *
 * Class contains a profile of aligned sequences.
 * All the sequences in the profile must have the same length.
 * This is used in ProfileAligner to build multiple or pairwise alignments.
 */

#ifndef PROFILE_H_
#define PROFILE_H_

#include <vector>
#include "aux.h"

using namespace std;


class Profile {
public:
	Profile(const vector<string> &,
			const vector<string> &,
			const vector<vector<double> > &, vector< vector<int> >, double);
	virtual ~Profile();
	Profile();
	Profile(const Profile &);
	//Profile& operator=(const Profile &);
	int getLength();
	void setLength(int);
	int getNumberOfSeqs();
	vector<string> getSeqs();
	vector<vector<double> > getHeights();
	vector<string> getIds();
	vector< vector<double> > getPSSM();
	vector<vector<int> > getPositions();
	void reverseProfile();
	vector<double> computeSumOfHeights();
	vector<double> getSumOfHeights();
	void printProfile(string);
	double getBp();
	vector< vector<double> > pssm;/*5 x m vector*/

private:
	vector<string> p; /*sequences*/
	vector<vector<double> > h; /*for each sequence there is an row vector of heights*/
	//vector< vector<double> > pssm;/*5 x m vector*/
	vector<double> h_sum;
	int n; /*length of the alignments in the profiles*/
	int m; /*number of sequences in each profile*/
	double bp; /*average base pairing probability of sequences in the profile*/
	vector<string> ids; /*sequence ids in the same order as sequences*/
	vector<vector<int> > positions;/*m by 2 vector of 0-indexed (start,stop) positions for each seq*/
	int checkLength();
	void checkPositions();
	void computePSSM();
};

#endif /* PROFILE_H_ */
