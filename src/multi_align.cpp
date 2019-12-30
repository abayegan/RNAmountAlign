/*
 * multialign.cpp
 *
 *  Created on: Dec 15, 2016
 *      Author: Amir Bayegan
 */

#include "multi_align.h"
//#define DEBUG

vector<string> multi_align(const Arguments &args){
	int i,j,n,l,MAXL;
	const int SEQNUM = args.seqs.size();
	//---get length of the longest sequence
	MAXL = -1;
	for (i = 0; i < SEQNUM; ++i) {
		l = args.seqs[i].length();
		if (l>MAXL) {
			MAXL = l;
		}
	}
	
	char * s = new char[MAXL+1];
	char * mfeStr = new char[MAXL+1];
	double ** bppr = allocate2DdoubleArray(MAXL+1,MAXL+1) ;
	double MFE,ensEng;
	double ** pairwiseSim= allocate2DdoubleArray(SEQNUM,SEQNUM);
	vector<double> h;
	vector<vector<double> > inc_h;
	vector<vector<double> > raw_h;
	vector<double> bpvec;//average base pairing probability of sequences
	
	double bp;
	double shift_str=1,mul_seq=1;

	vector<vector<Profile> > profvec;
	for(i = 0; i < SEQNUM; ++i){
		profvec.push_back(vector<Profile>(SEQNUM));
		for (j = 0; j < SEQNUM; ++j){
			profvec[i].push_back(Profile());
		}
	}
	

	for (i = 0; i < SEQNUM; ++i) {
		n = args.seqs[i].length();
		//fold and get base pairing probabilities
		strcpy(s,args.seqs[i].c_str());
		BasePairProbabilities(s,n,&MFE,mfeStr,&ensEng,args.argv0,bppr);
		
		bp=0;
		for(j=0;j<n;++j)
			bp+=bppr[j][n];
		bp/=n;
		
		//compute expected incremental heights
		if(args.modified_height){
			h = expectedModifiedHeight(bppr,n);
		}
		else{
			h = expectedHeight(bppr,n);
		}
		raw_h.push_back(h);
		h = heightDiff(h);
		inc_h.push_back(h);
		bpvec.push_back(bp);
	}
	scaleSimilarityScores(args,bpvec,shift_str,mul_seq);
	shift_str=0.5;
	//printf("scaling done...\n");
	BuildPairwiseDistanceMatrix(inc_h, bpvec, shift_str,mul_seq, SEQNUM, args,pairwiseSim,profvec);
	//printf("pairwise matrix done...\n");
	//flush(cout);
	Node * gtree = UPGMA(pairwiseSim,profvec, bpvec, shift_str, mul_seq ,SEQNUM, 0, inc_h, args);
	//clustalAlignment(gtree,args);
	
	delete [] s;
	delete [] mfeStr;
	free2DdoubleArray(bppr,n);
	return gtree->profile.getSeqs();
}

double seqIdentity(string seq1,string seq2){
	int i,n;
	int ident=0;
	if(seq1.length()!=seq2.length()){
		reportErr("Error in pairwiseIdentity()! Sequences do not have the same length");
	}
	n = seq1.length();
	for (i = 0; i < n; ++i) {
		if(seq1[i]==seq2[i] && seq1[i]!='-' && seq2[i]!='-'){
			ident+=1;
		}
	}
	return double(ident)/n;
}


void BuildPairwiseDistanceMatrix(vector<vector<double> > heights, vector<double> bpvec, double shift_str, double mul_seq ,int SEQNUM,  Arguments args, double ** pairwiseSim, vector<vector<Profile> > &profvec ){
	int i,j,k,l;
	

	vector<vector<double> > htmp;
	vector<string> stmp;
	vector<string> idtmp;
	vector<vector<int> > pos;

	for(i = 0; i < SEQNUM; ++i){
		for (j = i+1; j < SEQNUM; ++j) {
			pos.clear();htmp.clear();idtmp.clear();stmp.clear();
			htmp.push_back(heights[i]);
			stmp.push_back(args.seqs[i]);
			idtmp.push_back(args.ids[i]);
			vector<int> posi;
			posi.push_back(0);
			posi.push_back(args.seqs[i].length());
			pos.push_back(posi);
			Profile pi(stmp,idtmp,htmp,pos,bpvec[i]);

			pos.clear();htmp.clear();idtmp.clear();stmp.clear();
			htmp.push_back(heights[j]);
			stmp.push_back(args.seqs[j]);
			idtmp.push_back(args.ids[j]);
			vector<int> posj;
			posj.push_back(0);
			posj.push_back(args.seqs[i].length());
			pos.push_back(posj);
			Profile pj(stmp,idtmp,htmp,pos,bpvec[j]);
			args.aln_type='g';
			ProfileAligner pairaln= ProfileAligner(pi,pj,args,shift_str,mul_seq);
			profvec[i][j] = pairaln.align();
			
			
			
			//pairwiseDist[i][j] =  1.0-seqIdentity(profvec[i][j].getSeqs()[0],profvec[i][j].getSeqs()[1]);
			pairwiseSim[i][j] =  pairaln.getScore();
			pairwiseSim[j][i] = pairwiseSim[i][j];
			//printf("D(%d,%d)=%f\n",i,j,pairaln.getScore());
			//pairaln.printAlignment("clustal");
		}
	}
	//for(int i=0;i<profvec.size();i++){
			//for(int j=0;j<profvec[i].size();j++){
				//printf("%d,%d\n",i,j);
				//profvec[i][j].printProfile("fasta");
			//}
	//}
	return;
}

Node * UPGMA(double ** D, vector<vector<Profile> > & profvec, vector<double> bpvec, double shift_str, double mul_seq, int N, int weighted, const vector<vector<double> > & heights ,const Arguments &args ){
		/*------------------------------------------------------
	Declarations
	-------------------------------------------------------*/
	int ch;
	int i_size, j_size;
	/* size of i-th and j-th clusters to be amalgamated*/

	int num;    /* number of clusters */
	int i,i_index;    /* i_index index of new cluster to amalgamate */
	int j,j_index;    /* j_index index of new cluster to amalgamate */
	float d,dmax;   /* used to find closest clusters to amalgamate */
	vector<Node *> c(N);  /* c[i] is tree of i-th cluster */
	Node * p;

	double ** E = allocate2DdoubleArray(N,N);
	ProfileAligner *f; /*used to save and print the final alignment*/
	/*-------------------------------------------
	   E[i][j] distance between i-th and j-th cluster,
	   updated upon amalgamation of clusters.
	   copy D into E, E[i][j] distance between i-th
	   and j-th cluster.
	   -------------------------------------------*/
	E = D;
	/*-------------------------------------------
	   Initialize clusters and cluster sizes.
	   -------------------------------------------*/
	for (i=0;i<N;i++){
		c[i] = new Node();
		c[i]->leaf = 1;
		c[i]->depth = 0;
		c[i]->dist = 0;
		c[i]->index = i;
		c[i]->size = 1;
		c[i]->left = NULL;
		c[i]->right = NULL;
		vector<vector<double> > htmp;
		vector<string> stmp;
		vector<string> idtmp;
		vector<vector<int> > pos;
		htmp.push_back(heights[i]);
		stmp.push_back(args.seqs[i]);
		idtmp.push_back(args.ids[i]);
		vector<int> pos0;
		pos0.push_back(0);
		pos0.push_back(args.seqs[i].length());
		pos.push_back(pos0);
		Profile p0(stmp,idtmp,htmp,pos,bpvec[i]);
		c[i]-> profile = p0;
	}
	
	/*-------------------------------------------
	   Using E[i][j] find the closest 2 clusters, amalgamate,
	   update array of clusters, size of clusters, and E[i][j]'s.
	   Repeatedly do this, until there is only 1 cluster left.
	    -------------------------------------------*/
	num = N;
	while ( num > 1 )
	{
		dmax = -INF;
		for (i=0;i<N-1;i++)
			for (j=i+1;j<N;j++)
			{
				if (( c[i] != NULL ) && ( c[j] != NULL ))
				{
					d = E[i][j];
					if ( d > dmax )
					{
						dmax = d;
						i_index = i;
						j_index = j;
					}
				}
			}
		/*-------------------------------------------
	   Create new node, left child is c[i], right child is c[j].
	   Value of p node is distance between c[i] and c[j] clusters.
	   c[i] is now p, the amalgamated cluster (i<j), and index = 0
	   for non-leaf nodes.
	    -------------------------------------------*/
		i = i_index;
		j = j_index;
		p = new Node;
		p->leaf = 0;
		p->depth = 1 + MAX(c[i]->depth,c[j]->depth);
		p->dist = E[i][j];
		p->index = -1; /* important */
		i_size = c[i]->size;    /* size of i-th cluster */
		j_size = c[j]->size;   /* size of j-th cluster */
		p->size = i_size + j_size;
		p->left = c[i];
		p->right = c[j];
		
		
		if(c[i]->profile.getNumberOfSeqs()==1 && c[j]->profile.getNumberOfSeqs()==1 && args.aln_type=='g'){//reuse the pairwise alignments
			p->profile=profvec[i][j];
			//p->profile.printProfile("clustal");
			//for (int i=0; i< p->profile.getSumOfHeights().size();++i)
				//printf("%.2f ",p->profile.getSumOfHeights()[i]);
			//printf("\n");
		}
		else{
			f = new ProfileAligner(c[i]->profile, c[j]->profile,args, shift_str, mul_seq);
			p->profile = f->align();
			//p->profile.printProfile("clustal");
			//for (int i=0; i< p->profile.getSumOfHeights().size();++i)
				//printf("%.2f ",p->profile.getSumOfHeights()[i]);
			//printf("\n");
		}
		c[i] = p;
		c[j] = NULL;
		num--;
		/*-------------------------------------------
	   Update E.
	    -------------------------------------------*/
		for (j=0;j<N;j++)
			if ( j != j_index )
			{
				if ( weighted == 0 ) /* UPGMA */
					E[i_index][j] = (i_size * E[i_index][j] + j_size * E[j_index][j])/(float)(i_size + j_size);
				else   /* WPGMA */
					E[i_index][j] = (E[i_index][j] + E[j_index][j])/2;
			}
		E[i_index][i_index] = 0;
		for (j=0;j<N;j++)
			E[j][i_index] = E[i_index][j];   /* symmetrify */
		for (i=0;i<N;i++)
			E[i][j_index] = 0;
		for (j=0;j<N;j++)
			E[j_index][j] = 0;   /* remove j_index column and row */
		//delete p;
		//delete f;
	}    /* end of while loop */
	//delete p;
	f->printAlignment(args.outformat);
	return c[i_index];
}



//void MultiAlign::WriteNewick(Node * p){
//	p.
//
//}

/*-------------------------------------------
print phylogeny tree. When call, n should be depth
of tree.
-------------------------------------------*/
void printTree( Node * p, int n)
{
	int TABSIZE=8;
	int i;
	int depth( Node * );

	if ( p->leaf == 1 )
	{
		for (i=0;i<n-1;i++) putchar(9);   /* 9 is tab code */
		for (i=0;i<TABSIZE;i++) printf("-");
		printf("%d\n",p->index);
	}
	else
	{
		printTree(p->right,n+1);
		/*----------------- vertical bar towards parent ---*/
		for (i=0;i<n;i++) putchar(9);   /* 9 is tab code */
		printf("|\n");

		for (i=0;i<n-1;i++) putchar(9);   /* 9 is tab code */
		for (i=0;i<TABSIZE;i++) printf("-");
		printf("%2.2f\n",p->dist/(float) 2);
		/* ancestor is half distance from children */

		/*----------------- vertical bar towards parent ---*/
		for (i=0;i<n;i++) putchar(9);   /* 9 is tab code */
		printf("|\n");

		printTree(p->left,n+1);
	}
}

void scaleSimilarityScores(const Arguments &args, const vector<double> & bpvec,double &shift_str, double &mul_seq){
	double bpl,bp;
	double sum=0;
	int n=bpvec.size();
	double nucpr[4]={0,0,0,0};
	int m=0;
	double mu_seq,std_seq,mu_str,std_str;

	//compute average base pairing probability over all sequences
	for (int i = 0; i < n; ++i) {
		sum+=bpvec[i];
	}
	bp = sum/n;

	//compute nucleotide probability assuming a single long sequence is given
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < args.seqs[i].length(); ++j) {
			nucpr[char2num(args.seqs[i][j])]++;
			m++;
		}
	}
	for (int i = 0; i < 4; ++i){
		nucpr[i]/=m;
	}

	//compute mean and std of sequence similarity scores
	for (int nuc0 = 0; nuc0 < 4; ++nuc0) {
			for (int nuc1 = 0; nuc1 < 4; ++nuc1) {
				mu_seq += nucpr[nuc0]*nucpr[nuc1]*args.subMat[nuc0][nuc1];
			}
	}
	for (int nuc0 = 0; nuc0 < 4; ++nuc0) {
			for (int nuc1 = 0; nuc1 < 4; ++nuc1) {
				std_seq+= nucpr[nuc0]*nucpr[nuc1]*(args.subMat[nuc0][nuc1]-mu_seq)*(args.subMat[nuc0][nuc1]-mu_seq);
			}
	}
	std_seq = sqrt(std_seq);

	//compute mean and std of structure similarity scores
	bpl=bp/2.0;//base pairing to the left or right
	mu_str = -1.0*(4*bpl*(1.0-bp))+ //#"(" with "." or ")" with "."
					-2.0*(2*bpl*bpl); //	#"(" with ")"


	std_str = (-1.0 - mu_str)*(-1.0 - mu_str)*(4*bpl*(1.0-bpl))+ //#"(" with "." or ")" with "."
			(-2.0 - mu_str)*(-2.0 - mu_str)*(2*bpl*bpl); //	#"(" with ")"

	std_str= sqrt(std_str);

	if (!equal(std_seq,0.0) and !equal(std_str,0.0)){
	//		mult_seq = (stdStr/stdSeq)*(str_weight/(1-str_weight));
			mul_seq = std_str/std_seq;
		}
		shift_str = mul_seq*mu_seq - mu_str;
		sprintf(outtxt+strlen(outtxt),"#Scaling factors:multiplication=%f,shift=%f\n",mul_seq,shift_str);
		sprintf(outtxt+strlen(outtxt),"#Mean: structure:%f,sequence=%f\n",mu_str,mu_seq);
		sprintf(outtxt+strlen(outtxt),"#Std: structure=%f,sequence=%f\n\n",std_str,std_seq);
}








