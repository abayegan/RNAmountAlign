/*
 * multialign.h
 *
 *  Created on: Dec 15, 2016
 *      Author: Amir Bayegan
 */


#ifndef MULTIALIGN_H_
#define MULTIALIGN_H_

#include <vector>
#include <string>
#include <cmath>
#include "pair_align.h"
#include "aux.h"
#include "profile_aligner.h"
#include "mountain_height.h"
#include "profile.h"

using namespace std;

#define MAX(x,y)  ( (x)>(y) ? (x) : (y) )

vector<string> multialign(const Arguments &args);
double pairwiseDistance(string ,string);

struct Node
{
   int leaf;   /* leaf = 1 if true, else 0 */
   int depth;   /* depth of leaf is 0 */
   float  dist;    /* dist is distance between left, right clusters */
   int index;   /* name of species, if node is leaf, else -1 */
   int size;   /* size of current cluster */
   struct Node *left;
   struct Node *right;
   Profile profile;
   //Node();
   //Node(const Node & that);
};
vector<string> multi_align(const Arguments &args);
double seqIdentity(string seq1,string seq2);
void BuildPairwiseDistanceMatrix(vector<vector<double> > heights, vector<double> bpvec, double shift_str, double mul_seq ,int SEQNUM,  Arguments args, double ** pairwiseDist, vector<vector<Profile> > &profvec );
Node * UPGMA(double ** D, vector<vector<Profile> > &profvec, vector<double> bpvec,double shift_str, double mul_seq, int N, int weighted, const vector<vector<double> > & heights ,const Arguments &args );
void printTree( Node * p, int n);
void scaleSimilarityScores(const Arguments &args , const vector<double> & bpvec,double &shift_str, double &mul_seq);
#endif /* MULTIALIGN_H_ */
