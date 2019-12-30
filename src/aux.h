#ifndef AUX_H_
#define AUX_H_
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>  
#include <iterator>

#include "pair_align.h"


using namespace std;
#ifdef _WIN32
	#include <windows>
#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix)
#include <unistd.h>
#include <limits.h>
#include <math.h>
#endif

#define EPSILON 0.0001
const int LINLEN=60;
extern bool verbose;
double round_percision(double x,int percision);
extern char outtxt[1000000];
bool equal(double a, double b);
double readBasePairProb(string fname,double pa,double pc,double pg,double pu);
vector<string> mergeVectors(vector<string> v1,vector<string> v2);
double ** allocate2DdoubleArray(int row, int col);
char ** allocate2DcharArray(int row, int col);
void free2DdoubleArray(double ** array,int row);
void free2DcharArray(char ** array,int row);
void free2DAlignmentArray(Alignment **,int);
double max(double p, double q, double r, char * ptr );
template<typename T>
void printMatrix0(T ** mat,int n,int m);
void printMatrix(vector<vector<char> > mat,int n,int m);
void printMatrix1(double ** mat,int n,int m);
void printMatrix2(vector<vector<pair<int,int> > > mat,int n,int m);
void print2Dvec(vector<vector<double> > vec);
double getPercentileFromSortedVector(std::vector<double> vec,int percentile);
char* getExecPath(char* argv0);
void reportErr(const char * msg);
void reportWarning(const char * msg);
int char2num(char);
Alignment ** allocate2DAlignmentArray(int , int );
void getNucleotideProbabilities(string seq, vector<double>  & nucpr);
void generateRandomRNA(char * seq, int length, float pa,float pc,float pg,float pu);
int rouletteWheel(float *p, int n);
double sre_random();
float getGCcontent(string seq);
void writeVec(char * fname, vector<double> vec);

#endif /* AUX_H_ */
