#include "stats.h"

double computeEvalue_EVD(double x, double lambda, double mu ){
	return exp(-1.0 * lambda * (x - mu));
}

//This function is adapted from John D. Cook
double cdf_normal(double x){
	//A&S refers to Handbook of Mathematical Functions by Abramowitz and Stegun
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

double computePvalueNormal(double x, double sigma, double mu){
	double z=(x-mu)/sigma;
	return 1.0 - cdf_normal(z);
}

//EML = 0.1764 + 0.7991 * E_KA
double KAevalueRegressionTransform(double KA){
	//return 0.1764 + 0.7991 * KA;
	return  0.7991 * KA;
}

double computeEvalue(double x, double lambda, double K, int N){
	return K*N*exp(-1.0*lambda*x);
}

double computeK(double EVDlambda, double EVDmu,int n){
 return exp(EVDlambda*EVDmu)/n;
}

double comuteMu(double lambda, double K, int N){
	return log(K*N)/lambda;
}

double Eval2Pval(double e){
	return 1-exp(-e);
}

double Pval2Eval(double p){
	return -log(1-p);
}

double computeVectorMean(vector<double> vec){
	double sum=0.0;
	for (int i=0;i<vec.size();++i){
		sum+= vec[i];
	}
	return sum/vec.size();
}
double computeVectorStd(vector<double> vec, double mu){
	double sum=0.0;
	for (int i=0;i<vec.size();++i){
		sum+= (vec[i]-mu)*(vec[i]-mu);
	}
	return sqrt(sum/(vec.size()-1.0));
}


//given base pairing probability compute the probability of -|h_a'(i)-h_b'(j)| for two 
//positions i,j in sequences a,b resp. If only the MFE str is considered and NOT ensemble
//h' \in{1,0,-1}.
void absDiffBpProbability(const double &bp1, const double &bp2,score_prob &sc_pr){
	double pr,qr,pu,qu;
	pr=bp1/2.;
	qr=bp2/2.;
	pu=1.-bp1;
	qu=1.-bp2;
	
	//in both positions the same "(" , ")" or "."
	sc_pr.scmap[0]=pr*qr + pr*qr + pu*qu;
	
	//"(" and "." or ")" and "."
	sc_pr.scmap[-1]=pr*qu + qr*pu + pr*qu + qr*pu;
	
	//"(" and ")"
	sc_pr.scmap[-2]=pr*qr + pr*qr;
	sc_pr.checkProb();
}

void nucMutationProb(const vector<double> &npr1, const vector<double> &npr2 , score_prob &sc_pr, const float subMat[4][4] ){
	const char *nuc="ACGU";
	int i,j;
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
			if (i!=j){
				sc_pr.scmap[subMat[char2num(nuc[i])][char2num(nuc[j])]]= npr1[char2num(nuc[i])]*npr2[char2num(nuc[j])] + npr1[char2num(nuc[j])]*npr2[char2num(nuc[i])];
			}
			else{
				sc_pr.scmap[subMat[char2num(nuc[i])][char2num(nuc[i])]]=npr1[char2num(nuc[i])]*npr2[char2num(nuc[i])];
			}
		}
	}
	sc_pr.checkProb();
}


void karlinStats(const vector<double> &npr1,double **bppr1, int n1, const vector<double> &npr2, double **bppr2, int n2, const float subMat[4][4],float gamma,double &K, double &lambda){
	//printf("Beginning of karlinStats\n");
	double bpprob1,bpprob2;
	score_prob str_prob,seq_prob,comb_pr;
	double mu_seq=0.0,std_seq=0.0,mu_str=0.0,std_str=0.0;
	double mult_seq=1.0,shift_str=1.0;	
	int i,j;
	mapit str_it,seq_it;
	double seq_sc,seq_pr,str_sc,str_pr;
	
	bpprob1 = basePairProbability(bppr1,n1);
	bpprob2 = basePairProbability(bppr2,n2);
	//printf("BP probs: %f,%f\n",bpprob1,bpprob2);
	absDiffBpProbability(bpprob1,bpprob2,str_prob);
	//str_prob.print();
	nucMutationProb(npr1,npr2,seq_prob,subMat);
	//seq_prob.print();
	mu_str = str_prob.expectedValue();
	std_str = str_prob.standard_deviation(mu_str);
	
	mu_seq = seq_prob.expectedValue();
	std_seq = seq_prob.standard_deviation(mu_seq);
	
	//printf("museq=%f;stdseq=%f\n",mu_seq,std_seq);
	//printf("mustr=%f;stdstr=%f\n",mu_str,std_str);
	if (!equal(std_seq,0.0) and !equal(std_str,0.0)){
		mult_seq = std_str/std_seq;
	}
	shift_str = mult_seq*mu_seq - mu_str;
	//printf("shift_str=%f,mult_seq=%f\n",shift_str,mult_seq);
	for ( seq_it = seq_prob.scmap.begin(); seq_it != seq_prob.scmap.end(); ++seq_it ){
			seq_sc=seq_it->first;seq_pr=seq_it->second;
		for ( str_it = str_prob.scmap.begin(); str_it != str_prob.scmap.end(); ++str_it ){
			str_sc=str_it->first;str_pr=str_it->second;
			if(gamma<0){
				gamma = (bpprob1+bpprob2)/2.0;
			}
			comb_pr.scmap[gamma*(str_sc+shift_str) + (1.0-gamma)*(seq_sc*mult_seq)]= seq_pr*str_pr;
		}
	}
	//comb_pr.checkProb();
	//comb_pr.print();
	//printf("expScore:%f\n",comb_pr.expectedValue());
	karlin_altschul(comb_pr,K,lambda);
}

void computeStats_bin_gc( Arguments & args, vector<float> & gc_vec, vector<vector<double> > & gc_mu_lam_vec,const char distr){
	int i;
	
	float pa,pc,pg,pu;
	
	//int randlen = 2*qlen;
	double lambda,mu,K;
	float gc_max,gc_min,gc,gc_low,gc_high;
	int qlen = args.seqs[0].size();
	int randlen = args.window_size;
	vector<double> topscore_vec;
	vector<double> row;
	char fname[100];
	if (!args.qsearch)
		randlen=args.seqs[1].length();
	//variables for computing the expected heights
	double MFEq,ensEngq,MFEr,ensEngr;
	double ** bpprq = allocate2DdoubleArray(qlen,qlen+1);
	double ** bpprref = allocate2DdoubleArray(randlen,randlen+1);
	char * sq = new char[qlen+1];
	char * sr = new char[randlen+1];
	char * mfeStrq = new char[qlen+1];
	char * mfeStrr = new char[randlen+1];
	char * randseq = new char[randlen+1];
	vector<double> hq,hr;

	strcpy(sq,args.seqs[0].c_str());
	BasePairProbabilities(sq,qlen,&MFEq,mfeStrq,&ensEngq,args.argv0,bpprq);
	//compute heights for the query sequence
	if(args.modified_height){
		hq = expectedModifiedHeight(bpprq,qlen);
	}
	else{
		hq = expectedHeight(bpprq,qlen);
	} 
	hq=heightDiff(hq);
	sort(gc_vec.begin(),gc_vec.end());
	gc_min = gc_vec.front();
	gc_max = gc_vec.back();
	for (gc_low=gc_min;gc_low<=gc_max;gc_low+=args.gc_seg){
		gc_high = gc_low+ args.gc_seg;
		if(gc_high>gc_max) 
			gc_high=gc_max;
		gc = (gc_high + gc_low)/2.0;	//take the middle point in each GC range

		pc = pg = 0.5 *gc;
		pa = pu = 0.5 *(1-gc);
		
		row.clear();
		topscore_vec.clear();
		for(i=0 ; i < args.rnd_num; ++i){ //start sampling random RNAs
			generateRandomRNA(randseq,randlen,pa,pc,pg,pu);
			//fold the RNA and get heights
			BasePairProbabilities(randseq,randlen,&MFEr,mfeStrr,&ensEngr,args.argv0,bpprref);
			if(args.modified_height){
				hr = expectedModifiedHeight(bpprref,randlen);
			}
			else{
				hr = expectedHeight(bpprref,randlen);
			}
			hr = heightDiff(hr);
			//perform the local alignment
			Alignment aligner(args.seqs[0],randseq,hq,hr,bpprq,bpprref,args);
			aligner.align();
			topscore_vec.push_back(aligner.getScore());
			//aligner.printAlignment();
		}
		//write random scores to output
		//sprintf(fname,"%s_%.2f_%d_random_scores",args.fasta_name.c_str(),gc,args.rnd_num);
		//writeVec(fname,topscore_vec);
		
		if (distr=='N'){
			mu = computeVectorMean(topscore_vec);
			lambda = computeVectorStd(topscore_vec,mu);
		}
		else if (distr=='E'){
			//fit scores to EVD distributions and get the parameters
			EVD g(topscore_vec,true);
			if (!g.fit()) reportErr("Error in max likelihood fit of scores to EVD distribution\n");
			lambda = g.getLambda();
			K =  computeK(lambda,g.getMu(),randlen);
			mu = comuteMu(lambda,K,randlen);//IMPORTANT: don't mix up the EVD mu with Karlin/Altschul mu.
		}
		row.push_back(gc_low);row.push_back(gc_high);row.push_back(mu);row.push_back(lambda);
		gc_mu_lam_vec.push_back(row);
	}
	printStats( args.rnd_num,  randlen, gc_mu_lam_vec);	
	free2DdoubleArray(bpprq,qlen);
	free2DdoubleArray(bpprref,randlen);
	delete [] sq;
	delete [] sr;
	delete [] mfeStrq;
	delete [] mfeStrr;
}


//gc_vec: 
void get_params_for_gc(float gc, vector<vector<double> > & gc_mu_lam_vec, double & lambda , double & mu ){
	for (int i = 0; i < gc_mu_lam_vec.size(); ++i){
		if(gc_mu_lam_vec[i][0]<= gc && gc < gc_mu_lam_vec[i][1]){
			mu = gc_mu_lam_vec[i][2];
			lambda = gc_mu_lam_vec[i][3];
			return;
		}
		if((gc-gc_mu_lam_vec[i][1])<0.00001){ 		//The last bin is inclusive for both beg and end
			mu = gc_mu_lam_vec[i][2];
			lambda = gc_mu_lam_vec[i][3];
			return;
		}
	}
	reportErr("Error in get_params_for_gc(). GC content not found!\n");
}

void printStats(int rnd_num, int randlen,vector<vector<double> > gc_mu_lam_vec){
	int n = gc_mu_lam_vec.size();
	printf("#Stats:\n");
	
	printf("GC Bins: ");
	for (int i=0; i <n-1; ++i){
		printf("[%.2f-%.2f),",gc_mu_lam_vec[i][0],gc_mu_lam_vec[i][1]);
	}
	printf("[%.2f-%.2f]\n",gc_mu_lam_vec[n-1][0],gc_mu_lam_vec[n-1][1]);
	
	printf("%d random seqs of size %d generated for each each GC bin.\n",rnd_num,randlen);
	
	printf("\nFitting to EVD:\nGC_Content\tLocation_Param\tScale_Param\n");
	for (int i=0; i <n ; ++i){
		printf("%.3f\t%.2f\t%.2f\n",(gc_mu_lam_vec[i][0]+gc_mu_lam_vec[i][1])/2.0 , gc_mu_lam_vec[i][2],gc_mu_lam_vec[i][3]);
		}
	printf("\n");
}


