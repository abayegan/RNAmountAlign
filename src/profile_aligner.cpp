/*
 * ProfileAligner.cpp
 *
 *  Created on: Dec 18, 2016
 *      Author: Amir Bayegan
 */

#include "profile_aligner.h"
#include <time.h> 
ProfileAligner::~ProfileAligner() {
}

ProfileAligner::ProfileAligner(const Profile &pr1, const Profile &pr2, const Arguments & args, double shift, double mul) {

	p1 = pr1;
	p2 = pr2;
	len_aln=0;
	aln_start[0]=0;aln_start[1]=0;
	aln_stop[0]=p1.getLength();aln_stop[1]=p2.getLength();

	/* alignment params*/
	max_score=-INF;
	gap_init = args.gap_init;
	gap_ext = args.gap_ext;
	gap_init_str = args.gap_init_str;
	gap_ext_str = args.gap_ext_str;
	str_weight = args.str_weight;
	//cut_off = percentile;
	aln_type = args.aln_type;
	calc_prob=false;
	for (int i=0;i<=4;i++){
		for (int j=0;j<=4;j++){
			if (i!=4 && j!=4) {
				seq_similarity[i][j]=args.subMat[i][j];
			} else if (i==4 or j==4) {
				seq_similarity[i][j]=args.gap_ext; /* in computing SPS sim(X,'-')=Beta where X in{A,C,G,U}.*/
			} else{
				seq_similarity[i][j] =0; /* sim('-','-')=0 */
			}
		}
	}

	weighted_gap_init = (1-str_weight)*gap_init + str_weight*(gap_init_str);
	weighted_gap_ext = (1-str_weight)*gap_ext + str_weight*(gap_ext_str);
	
	/*auxiliary matrices*/
	vector< vector<double> > tmp0(p1.getLength()+1, vector<double>(p2.getLength()+1,0));
	scoreMatrix_p = tmp0;
	scoreMatrix_q = tmp0;
	scoreMatrix_r = tmp0;
	vector< vector<char> > tmp1(p1.getLength()+1, vector<char>(p2.getLength()+1,'N'));
	pathMatrix_p =tmp1;
	pathMatrix_q =tmp1;
	pathMatrix_r =tmp1;
	mult_seq=mul; shift_str=shift;
	//similarity = buildSimilarityMatrix();
//	printf("*********similarity matrix***********");
//	printMatrix1(similarity,p1.getLength(),p2.getLength());
//	for (int i=0;i<p1.getLength();++i)
//		printf("%f ",p1.getSumOfHeights()[i]);
//
//	printf("\n");
//	for (int i=0;i<p2.getLength();++i)
//			printf("%f ",p2.getSumOfHeights()[i]);
//	printf("\n");
}

int* ProfileAligner::getStart(){
	return aln_start;
}
int* ProfileAligner::getStop(){
	return aln_stop;
}

Profile ProfileAligner::align(){
	if (aln_type=='g'){//global alignment
			buildNeedlemanWunschScoreMatrix();
			backtrackNeedlemanWunsch();
		}
		else{//local alignment
			clock_t t;
			buildSmithWatermanScoreMatrix();
			backtrackSmithWaterman();
		}
	return p;
}

void ProfileAligner::buildNeedlemanWunschScoreMatrix(){
	//-----------initialize matrices
	int n1 = p1.getLength();
	int n2 = p2.getLength();
	scoreMatrix_p[0][0]=0;
	scoreMatrix_q[0][0]=0;
	scoreMatrix_r[0][0]=0;
	for(int i=1;i<=n1;i++){
		scoreMatrix_p[i][0]= gap(i);
		scoreMatrix_q[i][0]= -INF;
		scoreMatrix_r[i][0]= -INF;
		pathMatrix_p[i][0]='P';
		pathMatrix_q[i][0]='N'; //N is not used
		pathMatrix_r[i][0]='N';
	}
	for(int i=1;i<=n2;i++){
			scoreMatrix_p[0][i]= -INF;
			scoreMatrix_q[0][i]= gap(i);
			scoreMatrix_r[0][i]= -INF;
			pathMatrix_q[0][i]='Q';
			pathMatrix_p[0][i]='N';
			pathMatrix_r[0][i]='N';
	}
	//------Fill the matrices
	char ptr;
	double simij;
	for(int i=0;i<n1;i++){
		for(int j=0;j<n2;j++){
			scoreMatrix_p[i+1][j+1] = max(scoreMatrix_p[i][j+1]+weighted_gap_ext,
					   scoreMatrix_q[i][j+1]+weighted_gap_init,
					   scoreMatrix_r[i][j+1]+weighted_gap_init,
					   &ptr);
			pathMatrix_p[i+1][j+1]=ptr;
			scoreMatrix_q[i+1][j+1] = max(scoreMatrix_p[i+1][j]+weighted_gap_init,
						   scoreMatrix_q[i+1][j]+weighted_gap_ext,
						   scoreMatrix_r[i+1][j]+weighted_gap_init,
						   &ptr);
			pathMatrix_q[i+1][j+1]=ptr;
			simij = similarity(i,j);
			scoreMatrix_r[i+1][j+1] = max(scoreMatrix_p[i][j]+simij,
						   scoreMatrix_q[i][j]+simij,
						   scoreMatrix_r[i][j]+simij,
						   &ptr);
			pathMatrix_r[i+1][j+1]=ptr;
		}
	}


	#ifdef DEBUG
		printf("P:\n");
		print2Dvec(scoreMatrix_p);
		printMatrix(pathMatrix_p,n1+1,n2+1);
		printf("Q:\n");
		print2Dvec(scoreMatrix_q);
		printMatrix(pathMatrix_q,n1+1,n2+1);
		printf("R:\n");
		print2Dvec(scoreMatrix_r);
		printMatrix(pathMatrix_r,n1+1,n2+1);
	#endif
}

void ProfileAligner::backtrackNeedlemanWunsch(){
	char matrix;
	int n1 = p1.getLength();
	int n2 = p2.getLength();
	vector<string> s1(p1.getNumberOfSeqs()+p2.getNumberOfSeqs());
	vector<vector<double> > h1(p1.getNumberOfSeqs()+p2.getNumberOfSeqs());
	max_score = max(scoreMatrix_p[n1][n2],
			        scoreMatrix_q[n1][n2],
					scoreMatrix_r[n1][n2],
					& matrix);
	int x=n1;
	int y = n2;
	while(x>0 || y>0){
//		printf("%u,%u,%c\n",x,y,matrix);
		if(matrix=='R'){
			extendProfile(s1,h1,x-1,y-1,0);
			matrix = pathMatrix_r[x][y];
			x-=1;y-=1;
		}
		else if(matrix=='P'){
			extendProfile(s1,h1,x-1,y,2);
			matrix = pathMatrix_p[x][y];
			x-=1;
		}
		else if(matrix=='Q'){
			extendProfile(s1,h1,x,y-1,1);
			matrix = pathMatrix_q[x][y];
			y-=1;
		}
	}
	vector<string> alnids =  mergeVectors(p1.getIds(),p2.getIds());
	vector <vector<int> > pos = computeAlignmentPositions();
	Profile alnp(s1,alnids,h1,pos,(p1.getBp()+p2.getBp())/2.0);
	p = alnp;
	len_aln=p.getLength();
	p.reverseProfile();


}





void ProfileAligner::buildSmithWatermanScoreMatrix(){
	//-----------initialize matrices
	int n1 = p1.getLength();
	int n2 = p2.getLength();
	double sim;
	scoreMatrix_p[0][0]=0;
	scoreMatrix_q[0][0]=0;
	scoreMatrix_r[0][0]=0;
	for(int i=1;i<=n1;i++){
		scoreMatrix_p[i][0]= 0;
		scoreMatrix_q[i][0]= 0;
		scoreMatrix_r[i][0]= 0;
		pathMatrix_p[i][0]='P';
		pathMatrix_q[i][0]='N'; //N is not used
		pathMatrix_r[i][0]='N';
	}
	for(int i=1;i<=n2;i++){
			scoreMatrix_p[0][i]= 0;
			scoreMatrix_q[0][i]= 0;
			scoreMatrix_r[0][i]= 0;
			pathMatrix_q[0][i]='Q';
			pathMatrix_p[0][i]='N';
			pathMatrix_r[0][i]='N';
	}

	//------Fill the matrices
	char ptr;
	double maxv;
	double simij;
	for(int i=0;i<n1;i++){
		for(int j=0;j<n2;j++){
			maxv = max(scoreMatrix_p[i][j+1]+weighted_gap_ext,
					   scoreMatrix_q[i][j+1]+weighted_gap_init,
					   scoreMatrix_r[i][j+1]+weighted_gap_init,
					   &ptr);
			if(maxv<0){
				maxv = 0 ;
				ptr = 'N';
			}
			scoreMatrix_p[i+1][j+1]=maxv;
			pathMatrix_p[i+1][j+1]=ptr;
			maxv = max(scoreMatrix_p[i+1][j]+weighted_gap_init,
						   scoreMatrix_q[i+1][j]+weighted_gap_ext,
						   scoreMatrix_r[i+1][j]+weighted_gap_init,
						   &ptr);
			if(maxv<0){
				maxv = 0;
				ptr = 'N';
			}
			scoreMatrix_q[i+1][j+1]=maxv;
			pathMatrix_q[i+1][j+1]=ptr;
			
			simij = similarity(i,j);
			maxv = max(scoreMatrix_p[i][j]+simij,
						   scoreMatrix_q[i][j]+simij,
						   scoreMatrix_r[i][j]+simij,
						   &ptr);
			if(maxv<0){
				maxv = 0;
				ptr='N';
			}
			scoreMatrix_r[i+1][j+1] = maxv;
			pathMatrix_r[i+1][j+1]=ptr;
		}
	}

//	printf("%15c",' ');
//		for (int j = 0; j < scoreMatrix_p[0].size(); ++j)
//			printf("%15d",j);
//		printf("\n");
//		for (int i = 0; i < scoreMatrix_p.size(); ++i){
//			printf("%20d  ",i);
//			for (int j = 0; j < scoreMatrix_p[i].size(); ++j){
//					printf("%.2f,%.2f,%.2f  ",scoreMatrix_p[i][j],scoreMatrix_q[i][j],scoreMatrix_r[i][j]);
//			}
//	    printf("\n");
//		}

#ifdef DEBUG
	printf("P:\n");
	print2Dvec(scoreMatrix_p);
	printMatrix(pathMatrix_p,n1+1,n2+1);
	printf("Q:\n");
	print2Dvec(scoreMatrix_q);
	printMatrix(pathMatrix_q,n1+1,n2+1);
	printf("R:\n");
	print2Dvec(scoreMatrix_r);
	printMatrix(pathMatrix_r,n1+1,n2+1);
#endif
}


void ProfileAligner::backtrackSmithWaterman(){
	int n1 = p1.getLength();
	int n2 = p2.getLength();
	vector<string> s1(p1.getNumberOfSeqs()+p2.getNumberOfSeqs());
	vector<vector<double> > h1(p1.getNumberOfSeqs()+p2.getNumberOfSeqs());
	char matrix;
	int x,y,stop1,stop2;
	double score;
	max_score=0;

	for (int i=1; i<=n1;i++){
		for(int j=1;j<=n2;j++){
			if(scoreMatrix_p[i][j] > max_score){
				max_score = scoreMatrix_p[i][j];
				x=i;y=j;
				matrix='P';
			}
			if (scoreMatrix_q[i][j] > max_score){
				max_score = scoreMatrix_q[i][j];
				x=i;y=j;
				matrix='Q';
			}
			if (scoreMatrix_r[i][j] > max_score){
				max_score = scoreMatrix_r[i][j];
				x=i;y=j;
				matrix='R';
			}
		}
	}
	aln_stop[0]=x;aln_stop[1]=y;
	score = max_score;

	while(score>0){
		if(matrix=='R'){
			extendProfile(s1,h1,x-1,y-1,0);
			matrix = pathMatrix_r[x][y];
			x-=1;y-=1;
		}
		else if(matrix=='P'){
			extendProfile(s1,h1,x-1,y,2);
			matrix = pathMatrix_p[x][y];
			x-=1;
		}
		else if(matrix=='Q'){
			extendProfile(s1,h1,x,y-1,1);
			matrix = pathMatrix_q[x][y];
			y-=1;
		}
		if(matrix=='P'){
			score=scoreMatrix_p[x][y];
		}
		else if(matrix=='Q'){
			score = scoreMatrix_q[x][y];
		}
		else if(matrix=='R'){
			score = scoreMatrix_r[x][y];
		}
	}
	aln_start[0]=x;aln_start[1]=y;
	vector<string> alnids =  mergeVectors(p1.getIds(),p2.getIds());
	vector <vector<int> > pos = computeAlignmentPositions();
	Profile alnp(s1,alnids,h1,pos,(p1.getBp()+p2.getBp())/2.0);
	p = alnp;
	len_aln=p.getLength();
	p.reverseProfile();
}

float ProfileAligner::gap(int k){
	if (k<=0)
		return 0;
	else
		return (1-str_weight)*(gap_init+gap_ext*(k-1)) + str_weight*(gap_init_str+gap_ext_str*(k-1) + (shift_str*k)/2.0);
}

/*compute the value based on the given percentile
*to be used as MAX to convert distance of heights
to a measure of similarity in structural alignment.*/
double ProfileAligner::computeMaxHeight(){
	double max1=-INF,min1=INF,max2=-INF,min2=INF;
	for (int i = 0; i < p1.getNumberOfSeqs(); ++i) {
		for (int j = 0; j <	p1.getLength(); ++j) {
			if (p1.getHeights()[i][j]>max1) {
				max1=p1.getHeights()[i][j];
			}
			if (p1.getHeights()[i][j]<min1) {
				min1=p1.getHeights()[i][j];
			}
		}
	}
	for (int i = 0; i < p2.getNumberOfSeqs(); ++i) {
		for (int j = 0; j <	p2.getLength(); ++j) {
			if (p2.getHeights()[i][j]>max2) {
				max2=p2.getHeights()[i][j];
			}
			if (p2.getHeights()[i][j]<min2) {
				min2=p2.getHeights()[i][j];
			}
		}
	}

	return max(max2-min1,max1-min2);
}
/*add characters C to the final alignment from position x of p1
 * and position y of p2. If type=0 match, type=1: gap in p1, type=2: gap in p2*/
void ProfileAligner::extendProfile(vector<string> &s1,vector<vector<double> > &h1 ,int x, int y, int type){
	int i;
	int m1 = p1.getNumberOfSeqs();
	int m2 = p2.getNumberOfSeqs();

	/*-------------extend p and h----------*/
	if(type==0){
		for (i = 0; i < m1+m2; ++i) {
			if(i<m1){
				s1[i]+=p1.getSeqs()[i][x];
				h1[i].push_back(p1.getHeights()[i][x]);
			}
			else{
				s1[i]+=p2.getSeqs()[i-m1][y];
				h1[i].push_back(p2.getHeights()[i-m1][y]);
			}
		}
	}else if(type==1){
		for (i = 0; i < m1+m2; ++i) {
			if(i<m1){
				s1[i]+='-';
				h1[i].push_back(0);
			}
			else{
				s1[i] += p2.getSeqs()[i-m1][y];
				h1[i].push_back(p2.getHeights()[i-m1][y]);
			}
		}
	}else if(type==2){
		for (i = 0; i < m1+m2; ++i) {
			if(i<m1){
				s1[i]+=p1.getSeqs()[i][x];
				h1[i].push_back(p1.getHeights()[i][x]);
			}
			else{
				s1[i] += '-';
				h1[i].push_back(0);
			}
		}
	}
}

double ** ProfileAligner::buildSimilarityMatrix(){
	double ** sim=allocate2DdoubleArray(p1.getLength(),p2.getLength());
	double ** seqSim=allocate2DdoubleArray(p1.getLength(),p2.getLength());
	int i,j;
	int n=p1.getLength();
	int m=p2.getLength();
	double sps;
	clock_t t;
	t=clock();
	for ( i = 0; i < n; ++i) {
			for ( j = 0; j < m; ++j) {
				sps=0.0;
				for (int nuc0 = 0; nuc0 < 4; ++nuc0) {
						for (int nuc1 = 0; nuc1 < 4; ++nuc1) {
							//printf("%d,%d,%d,%d\n",i,j,nuc0,nuc1);
							sps += p1.getPSSM()[nuc0][i]*p2.getPSSM()[nuc1][j]*seq_similarity[nuc0][nuc1];
							//sps += p1.pssm[nuc0][i]*p2.pssm[nuc1][j]*seq_similarity[nuc0][nuc1];
						}
				}
//				printf("i:%d,j:%d,h[i]:%f,h[j]:%f,sim[i][j]:%f = %f x %f\n",i,j,p1.getSumOfHeights()[i],p2.getSumOfHeights()[j],sim[i][j],sps,mult_seq);
				sim[i][j] = str_weight*
						    (- fabs(p1.getSumOfHeights()[i] - p2.getSumOfHeights()[j])+shift_str) +
							(1.0-str_weight)*sps*mult_seq;
//
			}
	}
	t = clock()-t;
	cout<< ((float)t)/CLOCKS_PER_SEC << endl;
	fflush(stdout);
	return sim;
}

double ProfileAligner::similarity(int i,int j){
	double sps,sim;
	clock_t t;
	t=clock();
	sps=0.0;
	for (int nuc0 = 0; nuc0 < 4; ++nuc0) {
			for (int nuc1 = 0; nuc1 < 4; ++nuc1) {
				//printf("%d,%d,%d,%d\n",i,j,nuc0,nuc1);
				//sps += p1.getPSSM()[nuc0][i]*p2.getPSSM()[nuc1][j]*seq_similarity[nuc0][nuc1];
				//sps=1;
				sps += p1.pssm[nuc0][i]*p2.pssm[nuc1][j]*seq_similarity[nuc0][nuc1];
			}
	}
//				printf("i:%d,j:%d,h[i]:%f,h[j]:%f,sim[i][j]:%f = %f x %f\n",i,j,p1.getSumOfHeights()[i],p2.getSumOfHeights()[j],sim[i][j],sps,mult_seq);
	sim = str_weight*
			    (- fabs(p1.getSumOfHeights()[i] - p2.getSumOfHeights()[j])+shift_str) +
				(1.0-str_weight)*sps*mult_seq;
	//sim=1;
//

	t = clock()-t;
	//cout<< ((float)t)/CLOCKS_PER_SEC << endl;
	fflush(stdout);
	return sim;
}

double ProfileAligner::getSumOfPairsScore(int i,int j){
	double sps =0;
	vector<vector<double> > pssm1=p1.getPSSM();
	vector<vector<double> > pssm2=p2.getPSSM();
	for (int nuc0 = 0; nuc0 < 4; ++nuc0) {
		for (int nuc1 = 0; nuc1 < 4; ++nuc1) {
			sps += pssm1[nuc0][i]*pssm2[nuc1][j]*seq_similarity[nuc0][nuc1];
		}

	}
	return sps;
}

void ProfileAligner::computeQuality(){
	if(calc_prob){
		reportErr("partition function and probabilities not implemented yet!");
	}
	else{
		for (int i=0;i<len_aln;i++)
			quality+='|';
	}
}

void ProfileAligner::printAlignment(string format){
	sprintf(outtxt+strlen(outtxt),"#Score: %.2f\n\n",max_score );
	p.printProfile(format);
}

vector<vector<int> >  ProfileAligner::mergePositions(vector<vector<int> > v1,vector<vector<int> > v2){
	vector<vector<int> > v;
	v.insert( v.end(), v1.begin(), v1.end() );
	v.insert( v.end(), v2.begin(), v2.end() );
	return v;
}

vector<vector<int> > ProfileAligner::computeAlignmentPositions(){
	int i;
	int m1 = p1.getNumberOfSeqs();
	int m2 = p2.getNumberOfSeqs();
	vector<vector<int> > v(m1+m2,vector<int>(2,0));
	vector <vector<int> > pos1 = p1.getPositions();
	vector <vector<int> > pos2 = p2.getPositions();
	if(aln_type=='g'){
		v = mergePositions(pos1,pos2);
	}
	else{
		for (i = 0; i <m1+m2 ; ++i) {
			if(i<m1){
				v[i][0] = pos1[i][0]+ aln_start[0] ;
				v[i][1] = pos1[i][0]+ aln_stop[0] ;
			}
			else{
				v[i][0] = pos2[i-m1][0]+ aln_start[1] ;
				v[i][1] = pos2[i-m1][0]+ aln_stop[1] ;
			}
		}
	}
	return v;
}

double ProfileAligner::getScore(){
	return max_score;
}


