#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <sstream>
#include "arguments.h"
#include "pair_align.h"
#include "multi_align.h"
#include "mountain_height.h"
#include "aux.h"
#include "stats.h"

char outtxt[1000000];

int main(int argc , char * argv[]){
	int m,i;
	//parse the input parameters
	Arguments args;
	args.parseArgs(argc,argv);
	m =args.seqs.size();
	vector<string> aln;
	string a = args.sprintArgs().c_str();
	strcat(outtxt,a.c_str());
	if(args.seqs.size()==2){//pairwise alignment
		double evalue,pval;
		int n1 = args.seqs[0].length();
		int n2 = args.seqs[1].length();
		
		vector<vector<double> > gc_mu_lam_vec; //stores mu and lambda for each gc bin
		//compute expected heights
		double ** bppr1 = allocate2DdoubleArray(n1,n1+1);
		double ** bppr2 = allocate2DdoubleArray(n2,n2+1);
		double MFE1,ensEng1,MFE2,ensEng2;
		char * s1 = new char[n1+1];
		char * s2 = new char[n2+1];
		char * mfeStr1 = new char[n1+1];
		char * mfeStr2 = new char[n2+1];
		strcpy(s1,args.seqs[0].c_str());
		strcpy(s2,args.seqs[1].c_str());
		vector<double> h1,h2,raw_h1,raw_h2;
		BasePairProbabilities(s1,n1,&MFE1,mfeStr1,&ensEng1,argv[0],bppr1);
		BasePairProbabilities(s2,n2,&MFE2,mfeStr2,&ensEng2,argv[0],bppr2);
		if(args.modified_height){
			raw_h1 = expectedModifiedHeight(bppr1,n1);
			raw_h2 = expectedModifiedHeight(bppr2,n2);
		}
		else{
			raw_h1 = expectedHeight(bppr1,n1);
			raw_h2 = expectedHeight(bppr2,n2);
		}

		h1 = heightDiff(raw_h1);
		h2 = heightDiff(raw_h2);
		
		//Align the sequences
		Alignment aligner(args.seqs[0],args.seqs[1],h1,h2,bppr1,bppr2,args);
		aligner.align();
		
		//Statistics
		if(args.stats and args.aln_type=='l'){
			if(args.evd){
				vector<float> gc_vec;
				gc_vec.push_back(getGCcontent(args.seqs[0]));
				computeStats_bin_gc(args,gc_vec,gc_mu_lam_vec,'E');
				evalue = computeEvalue_EVD(aligner.getScore(),gc_mu_lam_vec[0][3],gc_mu_lam_vec[0][2]);
				pval = Eval2Pval(evalue);
			}
			else{
				vector <double> nucpr1, nucpr2;
				double k,lambda;
				getNucleotideProbabilities(s1,nucpr1);
				getNucleotideProbabilities(s2,nucpr2);
				karlinStats(nucpr1,bppr1,n1,nucpr2,bppr2,n2,args.subMat,args.str_weight,k,lambda);
				evalue = computeEvalue(aligner.getScore(),lambda,k,n1*n2);
				evalue = KAevalueRegressionTransform(evalue);
				pval = Eval2Pval(evalue);
			}
		}
		else if(args.stats){
			vector<float> gc_vec;
			gc_vec.push_back(getGCcontent(args.seqs[0]));
			computeStats_bin_gc(args,gc_vec,gc_mu_lam_vec,'N');
			pval = computePvalueNormal(aligner.getScore(),gc_mu_lam_vec[0][3],gc_mu_lam_vec[0][2]);
			evalue = Pval2Eval(pval);
		}
		
		// Output
		if(verbose){
			sprintf(outtxt+strlen(outtxt),"#%s:\n%s\n%s (%f)\n",args.ids[0].c_str(),args.seqs[0].c_str(),mfeStr1,MFE1);
			sprintf(outtxt+strlen(outtxt),"#%s:\n%s\n%s (%f)\n",args.ids[1].c_str(),args.seqs[1].c_str(),mfeStr2,MFE2);
			sprintf(outtxt+strlen(outtxt),"\n#Ensemble expected heights:\n#%s:",args.ids[0].c_str());
			strcat(outtxt,sprintHeights(raw_h1).c_str());
			sprintf(outtxt+strlen(outtxt),"#%s:",args.ids[1].c_str());
			strcat(outtxt,sprintHeights(raw_h2).c_str());
			sprintf(outtxt+strlen(outtxt),"\n#Ensemble incremental heights:\n#%s:",args.ids[0].c_str());
			strcat(outtxt,sprintHeights(h1).c_str());
			sprintf(outtxt+strlen(outtxt),"#%s:",args.ids[1].c_str());
			strcat(outtxt,sprintHeights(h2).c_str());
		}
		if(args.stats and args.aln_type=='l'){
			if(args.evd){
				sprintf(outtxt+strlen(outtxt),"\n#Statistics:\nNumber of random alignments: %d\n",args.rnd_num);
				sprintf(outtxt+strlen(outtxt),"Random alignment scores: mean= %f;std= %f\n",gc_mu_lam_vec[0][2],gc_mu_lam_vec[0][3]);
				sprintf(outtxt+strlen(outtxt),"Extreme value distribution E-value: %f \n",evalue);
				sprintf(outtxt+strlen(outtxt),"Extreme value distribution P-value: %f \n",pval);
			}
			else{
				sprintf(outtxt+strlen(outtxt),"\n#Statistics:\nExtreme value distribution E-value: %f \n",evalue);
				sprintf(outtxt+strlen(outtxt),"Extreme value distribution P-value: %f \n",pval);
			}
		}
		else if(args.stats){
			sprintf(outtxt+strlen(outtxt),"\n#Statistics:\nNumber of random alignments: %d\n",args.rnd_num);
			sprintf(outtxt+strlen(outtxt),"Random alignment scores: mean= %f;std= %f\n",gc_mu_lam_vec[0][2],gc_mu_lam_vec[0][3]);
			sprintf(outtxt+strlen(outtxt),"Normal distribution E-value: %f \n",evalue);
			sprintf(outtxt+strlen(outtxt),"Normal distribution P-value: %f \n",pval);
		}
		a=aligner.sprintAlignment();
		strcat(outtxt+strlen(outtxt),a.c_str());
		aln.push_back(aligner.getSeq1());
		aln.push_back(aligner.getSeq2());
		
		
		free2DdoubleArray(bppr1,n1);
		free2DdoubleArray(bppr2,n2);
		delete [] s1;
		delete [] s2;
		delete [] mfeStr1;
		delete [] mfeStr2;
	}
	else{ //multiple alignment
		aln = multi_align(args);
	}
	
	if(args.alifold){
		double cons_mfe;
		int n = aln[0].length();
		char * cons_str = new char[n+1];
		vector<const char*> cstring_aln;
		cstring_aln.reserve(aln.size()+1);
		for(int i = 0; i < aln.size(); ++i)
			cstring_aln.push_back(aln[i].c_str());
		cstring_aln[aln.size()]="\0";
		runAlifold(&cstring_aln[0], cons_str, cons_mfe);
		string cons(cons_str);
		const char * id="CONS_STR";
		if (args.outformat=="clustal")
		{
			int c = n / LINLEN;
			int d = n % LINLEN;
			int start=0;
			for(int i=0;i<c;i++){
				sprintf(outtxt+strlen(outtxt),"%-15s %-6s\n",id,cons.substr(start,LINLEN).c_str());
				start+=LINLEN;
			}
			if(d>0){
				sprintf(outtxt+strlen(outtxt),"%-15s %-6s  (%-5.2f)\n",id,cons.substr(start,d).c_str(),cons_mfe);
			}
		}
		else if(args.outformat=="fasta"){
			int c = n / LINLEN;
			int d = n % LINLEN;
			int start=0;
			sprintf(outtxt+strlen(outtxt),">CONS_STR\n");
			for(int i=0;i<c;i++){
				sprintf(outtxt+strlen(outtxt),"%s\n",cons.substr(start,LINLEN).c_str());
				start+=LINLEN;
			}
			if(d>0){
				sprintf(outtxt+strlen(outtxt),"%s\n",cons.substr(start,LINLEN).c_str());
			}
		}
		delete cons_str;
	}
	if(args.output_name==""){
			cout << outtxt;
		}
		else{
			ofstream outfile;
			outfile.open(args.output_name.c_str());
			outfile << outtxt;
			outfile.close();
	}
	return 0;
}



