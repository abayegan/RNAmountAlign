#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <string>
#include <iterator>

#include "aux.h"

using namespace std;

#ifndef ARGUMENTS_H_
#define ARGUMENTS_H_

const int MAXSEQNUM=10000; //maximum number of seqs in fasta file
class Arguments{
public:
	float gap_init;
	float gap_ext;
	float gap_init_str;
	float gap_ext_str;
	float str_weight;
	int percentile;
	bool print_prob;
	bool normal;
	string eps_name;
	string fasta_name;
	string query_name;
	string target_name;
	string output_name;
	vector<string> seqs;
	vector<string> ids;
	char aln_type;
	const char * matrix;
	float subMat[4][4];
	string outformat;
	int window_size;
	int step_size;
	int rnd_num;
	float gc_seg;
	bool fixed_gamma;
	bool modified_height;
	char * exec_path;
	bool stats;
	bool evd;
	bool alifold;
	bool qsearch;
	char * argv0;
	
	Arguments();
	void parseArgs(const int argc, char * argv[]);
	void printUsage(const char * argv0);
	void printArgs();
	string sprintArgs();
	void readFasta(const string fname, vector<string> & seqs, vector<string> & ids);
	void loadRibosumMatrix(const char *);
};





#endif /* ARGUMENTS_H_ */
