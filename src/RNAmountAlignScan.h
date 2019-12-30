#ifndef RNAMOUNTALIGNSCAN_H_
#define RNAMOUNTALIGNSCAN_H_ 
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include "arguments.h"
#include "pair_align.h"
#include "mountain_height.h"
#include "aux.h"
#include "evd.h"
#include "stats.h"


struct Window{
	vector<int> aln_pos;
	int index; //starts from 0
	int start; //1-indexed
	int end; //1-indexed
	float gc;
	double lambda;
	double mu;
	double karlin_k;
	double karlin_lambda;
	double karlin_evalue;
	double evalue;
	double score;
	string seq;
	void print();
};

void Window::print(){
	printf("window:%d,(%d,%d),GC:%.2f,lambda:%.2f,mu:%.2f,eval:%e\n",index,start,end,gc,lambda,mu,evalue);
}
//	overload operator < to compare objects by Evalue
bool operator < ( const Window & l, const Window & r){ 
   return l.evalue < r.evalue;
   }


#endif /* RNAMOUNTALIGNSCAN_H_ */
