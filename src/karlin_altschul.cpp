#include "karlin_altschul.h"

void karlin_altschul(const score_prob &sc_pr, double & K, double &lambda ){
	lambda = computeLambda(sc_pr);
	
	K = computeK(sc_pr,lambda);
	//printf("Lambda=%f\n",lambda);
	//printf("K=%f\n",K);
}

double computeH(const score_prob &sc_pr,double lam){
	double H=0.0;
	double sc,pr;
	for ( mapit it = sc_pr.scmap.begin(); it != sc_pr.scmap.end(); ++it ){
		sc=it->first;pr=it->second;
		H+=sc*pr*exp(sc*lam);
	}
	return H*lam;
}

double computeK(const score_prob &sc_pr, double lam){
	//printf("#COMPUTING K...\n");
	//sc_pr.print();
	int MAXIT=50;
	double TH=0.001;
	int PER=2;
	score_prob round_scpr,prev_scpr;
	int i,j,m,k;
	double numerator=0.0,denom,Cstar;
	double sc,pr,prevsc,prevpr;
	mapit it,prev;
	
	if(equal(lam,0)){
		return numeric_limits<double>::quiet_NaN();
	}
	
	denom=computeH(sc_pr,lam);
	//printf("Denome=%f\n",denom);
	for ( it = sc_pr.scmap.begin(); it != sc_pr.scmap.end(); ++it ){
		prev_scpr.scmap[round_percision(it->first,PER)] += it->second;
		if(it->first<0)
			numerator +=  it->second*exp(it->first*lam);
		else
			numerator +=  it->second;
	}
	for (j=2;j<MAXIT;++j){
		//printf("***j=%d\n",j);
		double sum1=0.0,sum2=0.0;
		score_prob next_scpr;
		for ( it = sc_pr.scmap.begin(); it != sc_pr.scmap.end(); ++it ){
			sc=round_percision(it->first,PER);pr=it->second;
			for ( prev = prev_scpr.scmap.begin(); prev != prev_scpr.scmap.end(); ++prev ){
				prevsc=prev->first;prevpr=prev->second;
				next_scpr.scmap[round_percision(prevsc + sc,PER)] += prevpr*pr;
			}
		}
		//next_scpr.print();
		//next_scpr.checkProb();
		//printf("exp=%f\n",next_scpr.expectedValue());
		for ( it = next_scpr.scmap.begin(); it != next_scpr.scmap.end(); ++it ){
			sc=it->first;pr=it->second;
			if(sc<0)
				sum1 += pr*exp(sc*lam);
			else
				sum2 += pr;
		}
		numerator+=(sum1+sum2)/j;
		//printf("sum1:%f,sum2:%f,numerator:%f\n",sum1,sum2,numerator);
		//if(exp(-2.0*numerator)/denom-Cstar < TH){
			//Cstar= exp(-2.0*numerator)/denom;
			//break;
		//}
		if(sum1+sum2 < TH){
			//printf("sum1+sum2=%f,TH=%f\n",sum1+sum2,TH);
			Cstar= exp(-2.0*numerator)/denom;
			break;
		}
		Cstar= exp(-2.0*numerator)/denom;
		prev_scpr= next_scpr;
	}
	return Cstar;
}

double computeLambda(const score_prob &sc_pr){
	double TH=0.00001;
	int MAXIT=100;
	double sum=0.0,up=0.5,point;
	int i;
	double sc,pr;
	mapit it;
	//find an upper bound for lambda
	while(sum<=1){
		sum=0.0;
		up = 2*up;
		for ( it = sc_pr.scmap.begin(); it != sc_pr.scmap.end(); ++it ){
			sc=it->first;pr=it->second;
			sum+=pr*exp(up*sc);
		}
	}
	//find lambda using binary search
	double low=0.0;
	for (i=0;i<MAXIT;++i){
		point = (up+low)/2.;
		sum=0.0;
		for (  it = sc_pr.scmap.begin(); it != sc_pr.scmap.end(); ++it ){
			sc=it->first;pr=it->second;
			sum+=pr*exp(point*sc);
		}
		if((abs(up-low) < TH)) break;
		if(sum>1)
			up=point;
		else
			low=point;
		//print("up=%f;low=%f;point=%f\n",up,low,point)
	}
	return low;
}
