#ifndef MOUNTAINHEIGHT_H_
#define MOUNTAINHEIGHT_H_

#include <math.h>
#include <stdlib.h>
#include <vector>
#include<iostream>
#include "aux.h"
extern "C"{
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "RNAstruct.h"
#include "cofold.h"
#include "part_func_co.h"
#include "read_epars.h"
#include "alifold.h"
}

using namespace std;


void BasePairProbabilities(char * seq , int n, double * mfe, char* mfeStr, double *ensEng,char* argv0, double ** bppr);
vector<double> expectedModifiedHeight(double ** bppr, int n);
vector<double> expectedHeight(double **bppr , int n);
vector<double> heightDiff(vector<double> h);
void readEnergyFile(char * argv0);
double basePairProbability(double **bppr , int n);
void printHeights(const vector<double> &h);
string sprintHeights(const vector<double> & h);
void runAlifold(const char ** seqs , char * str, double &mfe);
#endif /* MOUNTAINHEIGHT_H_ */
