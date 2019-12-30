#include "RNAmountAlignScan.h"


 
int main(int argc , char * argv[]){
	double ** bppr;
	char * mfeStr;
	double MFE,ensEng;
	char * s;
	int n,m,i,w;
	vector<Window> window_vec;
	Window window;
	string curr_seq;
	double lambda,mu,K,sp;
	vector<double> nucpr1,nucpr2;
	//parse the input parameters
	Arguments args;
	args.parseArgs(argc,argv);
	m =args.seqs.size();
	
	args.printArgs();
	
	
	if(args.seqs.size()==2){//only 2 sequences in the input
		int n1 = args.seqs[0].length();
		int n2 = args.seqs[1].length();
		if (args.window_size>n2){
			args.window_size = n2;
		}
		
		//variables for computing expected heights
		double ** bppr1 = allocate2DdoubleArray(n1,n1+1);
		double MFE2,ensEng2,MFE1,ensEng1;		
		char * s1 = new char[n1+1];
		char * mfeStr1 = new char[n1+1];
		strcpy(s1,args.seqs[0].c_str());
		vector<double> h1,h2;
		
		
		if (n2<n1){
			printf("Query sequence must be shorter than the target\n");
			exit(1);
		}
		if (args.window_size<n1){
			printf("Windows size must be higher than the query size\n");
			exit(1);
		}
		
		BasePairProbabilities(s1,n1,&MFE1,mfeStr1,&ensEng1,argv[0],bppr1);
		if(args.modified_height){
			h1 = expectedModifiedHeight(bppr1,n1);
		}
		else{
			h1 = expectedHeight(bppr1,n1);
		} 
		h1= heightDiff(h1);
		getNucleotideProbabilities(s1,nucpr1);
		
		//-----------------------Perform the alignment in sliding window
		for (i=0;i<n2;i+=args.step_size){
			w = args.window_size;
			if(i+w-1 >n2-1){
				w = n2-i;
				if(w<n1){break;}//last window is too small				
			}
			double ** bppr2 = allocate2DdoubleArray(w,w+1);
			char * s2 = new char[w+1];
			char * mfeStr2 = new char[w+1];
			curr_seq= args.seqs[1].substr(i,w);
			strcpy(s2,curr_seq.c_str());
			BasePairProbabilities(s2,w,&MFE2,mfeStr2,&ensEng2,argv[0],bppr2);
			if(args.modified_height){
				h2 = expectedModifiedHeight(bppr2,w);	
			}
			else{
				h2 = expectedHeight(bppr2,w);
			}
			h2 = heightDiff(h2);		
			window.seq = curr_seq;
			Alignment aligner(args.seqs[0],curr_seq,h1,h2,bppr1,bppr2,args);
			aligner.align();
			//aligner.printAlignment();
			window.score = aligner.getScore();
			window.aln_pos = aligner.getAlignmentPositions();
			window.index = i/args.step_size;
			window.start = i+1;
			window.end = i+w;
			window.gc = getGCcontent(curr_seq);
			if(args.aln_type=='l'){
				getNucleotideProbabilities(s2,nucpr2);
				karlinStats(nucpr2,bppr2,w,nucpr1,bppr1,n1,args.subMat,args.str_weight,window.karlin_k,window.karlin_lambda);
			}
			window_vec.push_back(window);
			//window.print();
			free2DdoubleArray(bppr2,w);
			delete [] s2;
			delete [] mfeStr2;
		}
		
		//------------------compute statistics from random sequences
		if(args.stats){
			if(args.aln_type!='l'){
				vector<float> gc_vec; //stores gc content of each window in reference
				vector<vector<double> > gc_mu_lam_vec; //stores mu and lambda for each gc bin
				for(i=0;i<window_vec.size();++i){
					gc_vec.push_back(window_vec[i].gc);
				}
				computeStats_bin_gc(args,gc_vec,gc_mu_lam_vec,'N');
				for(i=0;i<window_vec.size();++i){
					get_params_for_gc(window_vec[i].gc,gc_mu_lam_vec,window_vec[i].lambda,window_vec[i].mu);
					window_vec[i].evalue = computePvalueNormal(window_vec[i].score, window_vec[i].lambda,window_vec[i].mu);
				}
				
				sort(window_vec.begin(),window_vec.end());
				//-----------------print results
				printf("#Window Size:%d; Step Size:%d\n\n#Hits:\n",args.window_size,args.step_size);
				printf("%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n","Score","Window_Begin","Window_End","Hit_Begin","Hit_end","Pvalue(ND)");
				for (vector<Window>::iterator it=window_vec.begin(); it!= window_vec.end(); ++it){
					printf("%-12.2f\t%-12d\t%-12d\t%-12d\t%-12d\t%-12e\n",it->score,it->start,it->end,it->aln_pos[2]+it->start-1,it->aln_pos[3]+it->start-1,it->evalue);
				}
			}
			else{//local
				if(args.evd){
					vector<float> gc_vec; //stores gc content of each window in reference
					vector<vector<double> > gc_mu_lam_vec; //stores mu and lambda for each gc bin
					for(i=0;i<window_vec.size();++i){
						gc_vec.push_back(window_vec[i].gc);
					}
					computeStats_bin_gc(args,gc_vec,gc_mu_lam_vec,'E');
					for(i=0;i<window_vec.size();++i){
						get_params_for_gc(window_vec[i].gc,gc_mu_lam_vec,window_vec[i].lambda,window_vec[i].mu);
						window_vec[i].evalue = computeEvalue_EVD(window_vec[i].score, window_vec[i].lambda,window_vec[i].mu);
					}
					
					sort(window_vec.begin(),window_vec.end());
					//-----------------print results
					printf("#Window Size:%d; Step Size:%d\n\n#Hits:\n",args.window_size,args.step_size);
					printf("%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n","Score","Window_Begin","Window_End","Hit_Begin","Hit_end","E-value(EVD)");
					for (vector<Window>::iterator it=window_vec.begin(); it!= window_vec.end(); ++it){
						printf("%-12.2f\t%-12d\t%-12d\t%-12d\t%-12d\t%-12e\n",it->score,it->start,it->end,it->aln_pos[2]+it->start-1,it->aln_pos[3]+it->start-1,it->evalue);
					}	
				}
				else{//karlin-altschul
					for(i=0;i<window_vec.size();++i){
						sp = n1*(window_vec[i].end-window_vec[i].start+1);
						window_vec[i].karlin_evalue = computeEvalue(window_vec[i].score, window_vec[i].karlin_lambda,window_vec[i].karlin_k, sp);
						//printf("Eval=%df\n",window_vec[i].karlin_evalue);
					}
					//-----------------print results
					printf("#Window Size:%d; Step Size:%d\n\n#Hits:\n",args.window_size,args.step_size);
					printf("%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n","Score","Window_Begin","Window_End","Hit_Begin","Hit_end","Evalue(KA)");
					for (vector<Window>::iterator it=window_vec.begin(); it!= window_vec.end(); ++it){
						printf("%-12.2f\t%-12d\t%-12d\t%-12d\t%-12d\t%-12e\n",it->score,it->start,it->end,it->aln_pos[2]+it->start-1,it->aln_pos[3]+it->start-1,it->karlin_evalue);
					}
				}
			}
		}	
		else{//output only alignment scores
			printf("#Window Size:%d; Step Size:%d\n\n#Hits:\n",args.window_size,args.step_size);
				printf("%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n","Score","Window_Begin","Window_End","Aln_Begin","Aln_end");
				for (vector<Window>::iterator it=window_vec.begin(); it!= window_vec.end(); ++it){
					printf("%-12.2f\t%-12d\t%-12d\t%-12d\t%-12d\n",it->score,it->start,it->end,it->aln_pos[2]+it->start-1,it->aln_pos[3]+it->start-1);
				}	
		}
		
		delete [] s1;
		delete [] mfeStr1;
		free2DdoubleArray(bppr1,n1);
	}
	return 0;
}



