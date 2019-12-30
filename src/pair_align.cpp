#include "pair_align.h"


Alignment::Alignment(string seq1,string seq2,const vector<double> &h1 ,const vector<double> &h2, double b1, double b2, const Arguments &args) {
	seq1_raw = seq1;
	seq2_raw = seq2;
	len1_raw = seq1_raw.length();
	len2_raw = seq2_raw.length();
	id1= args.ids[0];
	id2 = args.ids[1];
	len_aln=0;
	gap_init = args.gap_init;
	gap_ext = args.gap_ext;
	gap_init_str = args.gap_init_str;
	gap_ext_str = args.gap_ext_str;
	str_weight = args.str_weight;
	aln_type = args.aln_type;
	seq1_aln = "";
	seq2_aln = "";
	vector< vector<double> > tmp0(len1_raw+1, vector<double>(len2_raw+1,0));
	scoreMatrix_p = tmp0;
	scoreMatrix_q = tmp0;
	scoreMatrix_r = tmp0;
	
	vector< vector<char> > tmp1(len1_raw+1, vector<char>(len2_raw+1,'N'));
	pathMatrix_p =tmp1;
	pathMatrix_q =tmp1;
	pathMatrix_r =tmp1;
	
	h1_raw = h1;
	h2_raw = h2;
	
	max_score=-INF;
	aln_start[0]=1;aln_start[1]=1;
	aln_stop[0]=len1_raw;aln_stop[1]=len2_raw;
	
	for (int i=0;i<=3;i++)
		for (int j=0;j<=3;j++)
			seq_similarity[i][j]=args.subMat[i][j];
	
	weighted_gap_init = (1-str_weight)*gap_init + str_weight*(gap_init_str);
	weighted_gap_ext = (1-str_weight)*gap_ext + str_weight*(gap_ext_str);
	outformat = args.outformat;
	fixed_gamma = args.fixed_gamma;
	//getBPprobs();
	bp1=b1;
	bp2=b2;
	scaleScores();
	//print base pair probabilities
	//printf("\nbase pair probabilities:\n");
	//printf("SEQ1:\n");
	//for (int i=0; i<len1_raw;i++)
	//{
		//for (int j=0;j<=len1_raw;j++)
			//printf	("%.5f\t",bppr1[i][j]);
		//printf("\n");
	//}
	//printf("SEQ2:\n");
		//for (int i=0; i<len2_raw;i++)
	//{
		//for (int j=0;j<=len2_raw;j++)
			//printf	("%.5f\t",bppr2[i][j]);
		//printf("\n");
	//}
}


Alignment::Alignment(string seq1,string seq2,const vector<double> &h1 ,const vector<double> &h2,
					double ** probmat1, double ** probmat2, const Arguments &args) {
	seq1_raw = seq1;
	seq2_raw = seq2;
	len1_raw = seq1_raw.length();
	len2_raw = seq2_raw.length();
	id1= args.ids[0];
	id2 = args.ids[1];
	len_aln=0;
	gap_init = args.gap_init;
	gap_ext = args.gap_ext;
	gap_init_str = args.gap_init_str;
	gap_ext_str = args.gap_ext_str;
	str_weight = args.str_weight;
	aln_type = args.aln_type;
	seq1_aln = "";
	seq2_aln = "";
	vector< vector<double> > tmp0(len1_raw+1, vector<double>(len2_raw+1,0));
	scoreMatrix_p = tmp0;
	scoreMatrix_q = tmp0;
	scoreMatrix_r = tmp0;
	
	vector< vector<char> > tmp1(len1_raw+1, vector<char>(len2_raw+1,'N'));
	pathMatrix_p =tmp1;
	pathMatrix_q =tmp1;
	pathMatrix_r =tmp1;
	
	bppr1 = allocate2DdoubleArray(len1_raw,len1_raw);
	bppr2 = allocate2DdoubleArray(len2_raw,len2_raw);
	bppr1 = probmat1;
	bppr2 = probmat2;
	
	h1_raw = h1;
	h2_raw = h2;
	
	max_score=-INF;
	aln_start[0]=1;aln_start[1]=1;
	aln_stop[0]=len1_raw;aln_stop[1]=len2_raw;
	
	for (int i=0;i<=3;i++)
		for (int j=0;j<=3;j++)
			seq_similarity[i][j]=args.subMat[i][j];
	
	weighted_gap_init = (1-str_weight)*gap_init + str_weight*(gap_init_str);
	weighted_gap_ext = (1-str_weight)*gap_ext + str_weight*(gap_ext_str);
	outformat = args.outformat;
	fixed_gamma = args.fixed_gamma;
	getBPprobs();
	scaleScores();
	//print base pair probabilities
	//printf("\nbase pair probabilities:\n");
	//printf("SEQ1:\n");
	//for (int i=0; i<len1_raw;i++)
	//{
		//for (int j=0;j<=len1_raw;j++)
			//printf	("%.5f\t",bppr1[i][j]);
		//printf("\n");
	//}
	//printf("SEQ2:\n");
		//for (int i=0; i<len2_raw;i++)
	//{
		//for (int j=0;j<=len2_raw;j++)
			//printf	("%.5f\t",bppr2[i][j]);
		//printf("\n");
	//}
	
}


Alignment::Alignment() {

}


Alignment::~Alignment() {

}

bool Alignment::operator < (const Alignment & aln) const{
	return (max_score < aln.max_score);
}

bool Alignment::operator > (const Alignment & aln) const{
	return (max_score > aln.max_score);
}

vector<int> Alignment::getAlignmentPositions(){
	vector<int> pos(4);
	pos[0] = aln_start[0];pos[1] = aln_stop[0];
	pos[2] = aln_start[1];pos[3] = aln_stop[1];
	return pos;
}

double Alignment::getScore(){
	return max_score;
}
string Alignment::getSeq1(){
	return seq1_aln;
}

string Alignment::getSeq2(){
	return seq2_aln;
}

int Alignment::getLength(){
	return len_aln;
}
char Alignment::getType(){
	return aln_type;
}

void Alignment::getBPprobs(){
	//bpppr[i][n] contains probability of i being base paired
	int i;
	bp1=0.0;
	for(i=0;i<len1_raw;++i)
		bp1+=bppr1[i][len1_raw];
	bp1/=len1_raw;
	
	bp2=0.0;
	for(i=0;i<len2_raw;++i)
		bp2+=bppr2[i][len2_raw];
	bp2/=len2_raw;
	return;
}

void Alignment::align(){
	if (aln_type=='g'){//global alignment
			buildNeedlemanWunschScoreMatrix();
			backtrackNeedlemanWunsch();
			computeQuality();
		}
		else if (aln_type=='l'){//local alignment
			buildSmithWatermanScoreMatrix();
			backtrackSmithWaterman();
			computeQuality();
		}
		else if (aln_type=='q'){//semi-global alignment
			buildSemiGlobalScoreMatrix();
			backtrackSemiGlobal();
			computeQuality();
		}
}
void Alignment::buildNeedlemanWunschScoreMatrix(){
	//-----------initialize matrices
	scoreMatrix_p[0][0]=0;
	scoreMatrix_q[0][0]=0;
	scoreMatrix_r[0][0]=0;
	for(int i=1;i<=len1_raw;i++){
		scoreMatrix_p[i][0]= gap(i);
		scoreMatrix_q[i][0]= -INF;
		scoreMatrix_r[i][0]= -INF;
		pathMatrix_p[i][0]='P';
		pathMatrix_q[i][0]='N'; //N is not used
		pathMatrix_r[i][0]='N';
	}
	for(int i=1;i<=len2_raw;i++){
			scoreMatrix_p[0][i]= -INF;
			scoreMatrix_q[0][i]= gap(i);
			scoreMatrix_r[0][i]= -INF;
			pathMatrix_q[0][i]='Q';
			pathMatrix_p[0][i]='N';
			pathMatrix_r[0][i]='N';
	}
	//------Fill the matrices
	char ptr;
	for(int i=0;i<len1_raw;i++){
		for(int j=0;j<len2_raw;j++){
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
			scoreMatrix_r[i+1][j+1] = max(scoreMatrix_p[i][j]+similarity(i,j,fixed_gamma),
						   scoreMatrix_q[i][j]+similarity(i,j,fixed_gamma),
						   scoreMatrix_r[i][j]+similarity(i,j,fixed_gamma),
						   &ptr);
			pathMatrix_r[i+1][j+1]=ptr;
		}
	}

	//printf("%15c",' ');
				//for (int j = 0; j < scoreMatrix_p[0].size(); ++j)
					//printf("%15d",j);
				//printf("\n");
				//for (int i = 0; i < scoreMatrix_p.size(); ++i){
					//printf("%20d  ",i);
					//for (int j = 0; j < scoreMatrix_p[i].size(); ++j){
							//printf("%.2f,%.2f,%.2f  ",scoreMatrix_p[i][j],scoreMatrix_q[i][j],scoreMatrix_r[i][j]);
					//}
			    //printf("\n");
				//}
				//printf("P:\n");
	#ifdef DEBUG
		printf("P:\n");
		print2Dvec(scoreMatrix_p);
		printMatrix(pathMatrix_p,len1_raw+1,len2_raw+1);
		printf("Q:\n");
		print2Dvec(scoreMatrix_q);
		printMatrix(pathMatrix_q,len1_raw+1,len2_raw+1);
		printf("R:\n");
		print2Dvec(scoreMatrix_r);
		printMatrix(pathMatrix_r,len1_raw+1,len2_raw+1);
	#endif
}

void Alignment::backtrackNeedlemanWunsch(){
	char matrix;
	max_score = max(scoreMatrix_p[len1_raw][len2_raw],
			        scoreMatrix_q[len1_raw][len2_raw],
					scoreMatrix_r[len1_raw][len2_raw],
					& matrix);
	int x=len1_raw;
	int y = len2_raw;
	while(x>0 || y>0){
//		printf("%u,%u,%c\n",x,y,matrix);
		//NOTE: sequence is 0-indexed but matrices have empty row and columns
		if(matrix=='R'){
			seq1_aln+=seq1_raw[x-1];
			seq2_aln+=seq2_raw[y-1];
			h1_aln.push_back(h1_raw[x-1]);
			h2_aln.push_back(h2_raw[y-1]);
			matrix = pathMatrix_r[x][y];
			x-=1;y-=1;
		}
		else if(matrix=='P'){
			seq1_aln+=seq1_raw[x-1];
			seq2_aln+='-';
			h1_aln.push_back(h1_raw[x-1]);
			h2_aln.push_back('-');
			matrix = pathMatrix_p[x][y];
			x-=1;
		}
		else if(matrix=='Q'){
			seq1_aln+='-';
			seq2_aln+=seq2_raw[y-1];
			h1_aln.push_back('-');
			h2_aln.push_back(h2_raw[y-1]);
			matrix = pathMatrix_q[x][y];
			y-=1;
		}
	}
	len_aln=h1_aln.size();
	reverse(seq1_aln.begin(), seq1_aln.end());
	reverse(seq2_aln.begin(), seq2_aln.end());
	reverse(h1_aln.begin(), h1_aln.end());
	reverse(h2_aln.begin(), h2_aln.end());
}

//while(x>0 || y>0){
////		printf("%u,%u,%c\n",x,y,matrix);
//		if(matrix=='R'){
//			seq1_aln+=seq1_raw[x-1];
//			seq2_aln+=seq2_raw[y-1];
//			h1_aln.push_back(h1_raw[x-1]);
//			h2_aln.push_back(h2_raw[y-1]);
//			matrix = pathMatrix_r[x][y];
//			x-=1;y-=1;
//		}
//		else if(matrix=='P'){
//			seq1_aln+=seq1_raw[x-1];
//			seq2_aln+='-';
//			h1_aln.push_back(h1_raw[x-1]);
//			h2_aln.push_back('-');
//			matrix = pathMatrix_p[x][y];
//			x-=1;
//		}
//		else if(matrix=='Q'){
//			seq1_aln+='-';
//			seq2_aln+=seq2_raw[y-1];
//			h1_aln.push_back('-');
//			h2_aln.push_back(h2_raw[y-1]);
//			matrix = pathMatrix_q[x][y];
//			y-=1;
//		}
//	}
//	len_aln=h1_aln.size();
//	reverse(seq1_aln.begin(), seq1_aln.end());
//	reverse(seq2_aln.begin(), seq2_aln.end());
//	reverse(h1_aln.begin(), h1_aln.end());
//	reverse(h2_aln.begin(), h2_aln.end());
//}

//#ifdef DEBUG
//	printf("P:\n");
//	printMatrix1(scoreMatrix_p,len1_raw+1,len2_raw+1);
//	printMatrix(pathMatrix_p,len1_raw+1,len2_raw+1);
//	printf("Q:\n");
//	printMatrix1(scoreMatrix_q,len1_raw+1,len2_raw+1);
//	printMatrix(pathMatrix_q,len1_raw+1,len2_raw+1);
//	printf("R:\n");
//	printMatrix1(scoreMatrix_r,len1_raw+1,len2_raw+1);
//	printMatrix(pathMatrix_r,len1_raw+1,len2_raw+1);
//#endif

void Alignment::buildSmithWatermanScoreMatrix(){
	//-----------initialize matrices
	scoreMatrix_p[0][0]=0;
	scoreMatrix_q[0][0]=0;
	scoreMatrix_r[0][0]=0;
	for(int i=1;i<=len1_raw;i++){
		scoreMatrix_p[i][0]= 0;
		scoreMatrix_q[i][0]= 0;
		scoreMatrix_r[i][0]= 0;
		pathMatrix_p[i][0]='P';
		pathMatrix_q[i][0]='N'; //N is not used
		pathMatrix_r[i][0]='N';
	}
	for(int i=1;i<=len2_raw;i++){
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
	for(int i=0;i<len1_raw;i++){
		for(int j=0;j<len2_raw;j++){
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

			maxv = max(scoreMatrix_p[i][j]+similarity(i,j,fixed_gamma),
						   scoreMatrix_q[i][j]+similarity(i,j,fixed_gamma),
						   scoreMatrix_r[i][j]+similarity(i,j,fixed_gamma),
						   &ptr);
			if(maxv<0){
				maxv = 0;
				ptr='N';
			}
			scoreMatrix_r[i+1][j+1] = maxv;
			pathMatrix_r[i+1][j+1]=ptr;
			
			//#ifdef DEBUG
			//printf("\n--------- i=%d , j=%d ---------------\n",i,j);
			//printf("P:\n");
			//printMatrix1(scoreMatrix_p,len1_raw+1,len2_raw+1);
			//printMatrix(pathMatrix_p,len1_raw+1,len2_raw+1);
			//printMatrix2(reset_p,len1_raw+1,len2_raw+1);
			//printf("Q:\n");
			//printMatrix1(scoreMatrix_q,len1_raw+1,len2_raw+1);
			//printMatrix(pathMatrix_q,len1_raw+1,len2_raw+1);
			//printMatrix2(reset_q,len1_raw+1,len2_raw+1);
			//printf("R:\n");
			//printMatrix1(scoreMatrix_r,len1_raw+1,len2_raw+1);
			//printMatrix(pathMatrix_r,len1_raw+1,len2_raw+1);
			//printMatrix2(reset_r,len1_raw+1,len2_raw+1);
			//#endif
		}
	}


	//printf("%15c",' ');
			//for (int j = 0; j < scoreMatrix_p[0].size(); ++j)
				//printf("%15d",j);
			//printf("\n");
			//for (int i = 0; i < scoreMatrix_p.size(); ++i){
				//printf("%20d  ",i);
				//for (int j = 0; j < scoreMatrix_p[i].size(); ++j){
						//printf("%.2f,%.2f,%.2f  ",scoreMatrix_p[i][j],scoreMatrix_q[i][j],scoreMatrix_r[i][j]);
				//}
		    //printf("\n");
			//}
			//printf("P:\n");

			//printMatrix(pathMatrix_p,len1_raw+1,len2_raw+1);
			//printf("Q:\n");

			//printMatrix(pathMatrix_q,len1_raw+1,len2_raw+1);
			//printf("R:\n");

			//printMatrix(pathMatrix_r,len1_raw+1,len2_raw+1);
}


void Alignment::backtrackSmithWaterman(){

	char matrix;
	int x,y,stop1,stop2;
	double score;
	max_score=0;

	for (int i=1; i<=len1_raw;i++){
		for(int j=1;j<=len2_raw;j++){
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
//	if(matrix>0){
//		seq1_aln+=seq1_raw[x-1];
//		seq2_aln+=seq2_raw[y-1];
//		h1_aln.push_back(h1_raw[x-1]);
//		h2_aln.push_back(h2_raw[y-1]);
//		x--;y--;
//	}
	while(score>0){
		if(matrix=='R'){
			seq1_aln+=seq1_raw[x-1];
			seq2_aln+=seq2_raw[y-1];
			h1_aln.push_back(h1_raw[x-1]);
			h2_aln.push_back(h2_raw[y-1]);
			matrix = pathMatrix_r[x][y];
			x-=1;y-=1;
		}
		else if(matrix=='P'){
			seq1_aln+=seq1_raw[x-1];
			seq2_aln+='-';
			h1_aln.push_back(h1_raw[x-1]);
			h2_aln.push_back('-');
			matrix = pathMatrix_p[x][y];
			x-=1;
		}
		else if(matrix=='Q'){
			seq1_aln+='-';
			seq2_aln+=seq2_raw[y-1];
			h1_aln.push_back('-');
			h2_aln.push_back(h2_raw[y-1]);
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
	aln_start[0]=x+1;aln_start[1]=y+1;
	len_aln=h1_aln.size();
	//printf("****%s\n%s\n;",seq1_aln.c_str(),seq2_aln.c_str());
	reverse(seq1_aln.begin(), seq1_aln.end());
	reverse(seq2_aln.begin(), seq2_aln.end());
	reverse(h1_aln.begin(), h1_aln.end());
	reverse(h2_aln.begin(), h2_aln.end());
	//printf("****%s\n%s\n;",seq1_aln.c_str(),seq2_aln.c_str());
}

void Alignment::buildSemiGlobalScoreMatrix(){
	//-----------initialize matrices
	scoreMatrix_p[0][0]=0;
	scoreMatrix_q[0][0]=0;
	scoreMatrix_r[0][0]=0;
	for(int i=1;i<=len1_raw;i++){
		scoreMatrix_p[i][0]= -INF;
		scoreMatrix_q[i][0]= gap(i);
		scoreMatrix_r[i][0]= -INF;
		pathMatrix_p[i][0]='Q';
		pathMatrix_q[i][0]='N'; //N is not used
		pathMatrix_r[i][0]='N';
	}
	for(int i=1;i<=len2_raw;i++){
			scoreMatrix_p[0][i]= 0;
			scoreMatrix_q[0][i]= -INF;
			scoreMatrix_r[0][i]= -INF;
			pathMatrix_q[0][i]='P';
			pathMatrix_p[0][i]='N';
			pathMatrix_r[0][i]='N';
	}
	//------Fill the matrices
	char ptr;
	for(int i=0;i<len1_raw;i++){
		for(int j=0;j<len2_raw;j++){
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
			scoreMatrix_r[i+1][j+1] = max(scoreMatrix_p[i][j]+similarity(i,j,fixed_gamma),
						   scoreMatrix_q[i][j]+similarity(i,j,fixed_gamma),
						   scoreMatrix_r[i][j]+similarity(i,j,fixed_gamma),
						   &ptr);
			pathMatrix_r[i+1][j+1]=ptr;
		}
	}


	#ifdef DEBUG
		printf("P:\n");
		printMatrix1(scoreMatrix_p,len1_raw+1,len2_raw+1);
		printMatrix(pathMatrix_p,len1_raw+1,len2_raw+1);
		printf("Q:\n");
		printMatrix1(scoreMatrix_q,len1_raw+1,len2_raw+1);
		printMatrix(pathMatrix_q,len1_raw+1,len2_raw+1);
		printf("R:\n");
		printMatrix1(scoreMatrix_r,len1_raw+1,len2_raw+1);
		printMatrix(pathMatrix_r,len1_raw+1,len2_raw+1);
	#endif
}

void Alignment::backtrackSemiGlobal(){

	char matrix;
	int x,y,stop1,stop2;
	double score;
	max_score=-INF;
	
	for (int i=1; i<=len2_raw;i++){
			if(scoreMatrix_p[len1_raw][i] > max_score){
				max_score = scoreMatrix_p[len1_raw][i];
				x=len1_raw;y=i;
				matrix='P';
			}
			if (scoreMatrix_q[len1_raw][i] > max_score){
				max_score = scoreMatrix_q[len1_raw][i];
				x=len1_raw;y=i;
				matrix='Q';
			}
			if (scoreMatrix_r[len1_raw][i] > max_score){
				max_score = scoreMatrix_r[len1_raw][i];
				x=len1_raw;y=i;
				matrix='R';
			}
	}
	aln_stop[0]=x;aln_stop[1]=y;
	score = max_score;

	while(x>0){
		//printf("%u,%u,%c\n",x,y,matrix);
		if(matrix=='R'){
			seq1_aln+=seq1_raw[x-1];
			seq2_aln+=seq2_raw[y-1];
			h1_aln.push_back(h1_raw[x-1]);
			h2_aln.push_back(h2_raw[y-1]);
			matrix = pathMatrix_r[x][y];
			x-=1;y-=1;
		}
		else if(matrix=='P'){
			seq1_aln+=seq1_raw[x-1];
			seq2_aln+='-';
			h1_aln.push_back(h1_raw[x-1]);
			h2_aln.push_back('-');
			matrix = pathMatrix_p[x][y];
			x-=1;
		}
		else if(matrix=='Q'){
			seq1_aln+='-';
			seq2_aln+=seq2_raw[y-1];
			h1_aln.push_back('-');
			h2_aln.push_back(h2_raw[y-1]);
			matrix = pathMatrix_q[x][y];
			y-=1;
		}
	}
	aln_start[0]=x+1;aln_start[1]=y+1;
	len_aln=h1_aln.size();
	reverse(seq1_aln.begin(), seq1_aln.end());
	reverse(seq2_aln.begin(), seq2_aln.end());
	reverse(h1_aln.begin(), h1_aln.end());
	reverse(h2_aln.begin(), h2_aln.end());
}

float Alignment::gap(int k){
	if (k<=0)
		return 0;
	else
		return (1-str_weight)*(gap_init+gap_ext*(k-1)) + str_weight*(gap_init_str+gap_ext_str*(k-1) + (shift_str*k)/2.0);//maximum score of structure similarity is SHIFT
}


double Alignment::similarity(int i , int j, bool fixed_gamma){
	double Gamma;
	float str_sim;
	float scaled_str_sim;
	
	int x = char2num(seq1_raw[i]);
	int y = char2num(seq2_raw[j]);
	
	if(fixed_gamma)
		Gamma = str_weight;
	else
		Gamma = ((1-bppr1[i][len1_raw]) + (1-bppr2[j][len2_raw]))/2.0;
		
	scaled_str_sim = -fabs(h1_raw[i]-h2_raw[j]) + shift_str;
	double dist = (Gamma*scaled_str_sim +(1-Gamma)*seq_similarity[x][y]*mult_seq);
	//printf("%d,%d : %.2f,%.2f,%.2f,%.2f\n",i,j,str_sim,scaled_str_sim,dist,Gamma);
	return dist;
}


void Alignment::computeQuality(){
	for (int i=0;i<len_aln;i++)
			quality+='|';
	//if(calc_prob){
		//reportErr("partition function and probabilities not implemented yet!");
	//}
	//else{
		//for (int i=0;i<len_aln;i++)
			//quality+='|';
	//}
}


void Alignment::scaleScores(){
	vector<double> nucpr1, nucpr2;
	double mu_seq=0.0,std_seq=0.0,mu_str=0.0,std_str=0.0;
	mult_seq=1.0;shift_str=1.0;
	getNucleotideProbabilities(seq1_raw,nucpr1);
	getNucleotideProbabilities(seq2_raw,nucpr2);
	computeMeanStdSeqScores(nucpr1,nucpr2,mu_seq,std_seq);
	computeMeanStdStrScores(nucpr1,nucpr2,mu_str,std_str);
	//printf("museq=%f;stdseq=%f\n",mu_seq,std_seq);
	//printf("mustr=%f;stdstr=%f\n",mu_str,std_str);
	if (!equal(std_seq,0.0) and !equal(std_str,0.0)){
		mult_seq = std_str/std_seq;
	}
	shift_str = mult_seq*mu_seq - mu_str;
	//printf("mul=%f;shift=%f\n",mult_seq,shift_str);
	return;
}


void Alignment::computeMeanStdSeqScores(vector<double> p,vector<double> q, double & mu, double & std){
	int i,j;
	for(i=0;i<p.size();++i){
		for(j=0;j<q.size();++j){
			mu+= p[i]*q[j]*seq_similarity[i][j];
		}
	}
	for(i=0;i<p.size();++i){
		for(j=0;j<q.size();++j){
			std+= p[i]*q[j]*(seq_similarity[i][j]-mu)*(seq_similarity[i][j]-mu);
		}
	}
	std = sqrt(std);
	return;	
}

void Alignment::computeMeanStdStrScores(vector<double> p,vector<double> q, double & mu, double & std){
	//Assuming heights of a fixed structure(not the ensemble) structure similarity can be -2,-1,0.
	int i,j;
	//printf("i:%d,j:%d,bp1:%f,bp2:%f\n",i,j,bp1,bp2);

	double pl1=bp1/2.0,pl2=bp2/2.0;
	mu = -1.0*(pl1*(1.0-bp2) + pl2*(1.0-bp1) + pl1*(1.0-bp2) + pl2*(1.0-bp1))+ //#"(" with "." or ")" with "."
	      -2.0*(pl1*pl2 + pl1*pl2); //	#"(" with ")"
	std = sqrt((pl1*pl2 + pl1*pl2 + (1.0-bp1)*(1-bp2))*(-mu)*(-mu) + //if the score is 0
			(pl1*(1.0-bp2) + pl2*(1.0-bp1) + pl1*(1.0-bp2) + pl2*(1.0-bp1))*(-1-mu)*(-1-mu)+ //if the score is -1
			(pl1*pl2 + pl1*pl2)*(-2-mu)*(-2-mu));//if the score is -2
	return;	
}


void Alignment::printAlignment(){
	string tmpstr1,tmpstr2;
	int gapnum1,gapnum2;
	char tmp1[2000];
	char tmp2[2000];
	string hh1="",hh2="";
	//printf("#%s: %s\n",id1.c_str(),seq1_raw.c_str());
	//printf("#%s: %s\n",id2.c_str(),seq2_raw.c_str());
	printf("\n#Scaling factors:\nmultiplication:%f\nshift:%f \n",mult_seq,shift_str);
	printf("\n#Score: %.2f\n\n",max_score);
	//print alignment
	printf("\n#Alignment\n");
	if(outformat=="clustal"){
		int c = len_aln / LINLEN;
		int d = len_aln % LINLEN;
		int offset1 = aln_start[0];
		int offset2 = aln_start[1];
		int start=0;
		for (int i=1;i<=c;i++){
			tmpstr1 = seq1_aln.substr(start,LINLEN);
			gapnum1 = count(tmpstr1.begin(), tmpstr1.end(), '-');
			tmpstr2 = seq2_aln.substr(start,LINLEN);
			gapnum2 = count(tmpstr2.begin(), tmpstr2.end(), '-');
			
			printf("%-15s %-5d %s  %-5d\n",id1.c_str(),offset1,seq1_aln.substr(start,LINLEN).c_str(),offset1+LINLEN-1-gapnum1);
			printf("%-20s  %s \n" ,"",quality.substr(start,LINLEN).c_str());
			printf("%-15s %-5d %s  %-5d\n\n",id2.c_str(),offset2,seq2_aln.substr(start,LINLEN).c_str(),offset2+LINLEN-1-gapnum2);
			start+=LINLEN;
			offset1+=LINLEN-gapnum1;
			offset2+=LINLEN-gapnum2;
		}
		if(d>0){
			tmpstr1 = seq1_aln.substr(start,d);
			gapnum1 = count(tmpstr1.begin(), tmpstr1.end(), '-');
			tmpstr2 = seq2_aln.substr(start,d);
			gapnum2 = count(tmpstr2.begin(), tmpstr2.end(), '-');
			printf("%-15s %-5d %s  %-5d\n",id1.c_str(),offset1,seq1_aln.substr(start,d).c_str(),offset1+d-1-gapnum1);
			printf("%-20s  %s \n" ,"",quality.substr(start,d).c_str());
			printf("%-15s %-5d %s  %-5d\n\n",id2.c_str(),offset2,seq2_aln.substr(start,d).c_str(),offset2+d-1-gapnum2);
		printf("\n");
		}
	}
	else if(outformat=="fasta"){
		int c = len_aln / LINLEN;
		int d = len_aln % LINLEN;
		int start=0;
		printf(">%s\n",id1.c_str() );
		for (int i=1;i<=c;i++){
			printf("%s\n",seq1_aln.substr(start,LINLEN).c_str());
			start+=LINLEN;
		}
		if(d>0){
			printf("%s\n",seq1_aln.substr(start,LINLEN).c_str());
		}
		start=0;
		printf(">%s\n",id2.c_str() );
		for (int i=1;i<=c;i++){
			printf("%s\n",seq2_aln.substr(start,LINLEN).c_str());
			start+=LINLEN;
		}
		if(d>0){
			printf("%s\n",seq2_aln.substr(start,LINLEN).c_str());
		}
	}
}

string Alignment::sprintAlignment(){
	string tmpstr1,tmpstr2;
	int gapnum1,gapnum2;
	char buff[10000]="";
	string com = "#";
	const int LINLEN=60;
	char tmp1[2000]="";
	char tmp2[2000]="";
	string hh1="",hh2="";

	sprintf(buff+strlen(buff),"\n#Scaling factors:\nmultiplication:%f\nshift:%f \n",mult_seq,shift_str);
	sprintf(buff+strlen(buff),"\n#Score: %.2f\n\n",max_score);
	
	//print alignment
	sprintf(buff+strlen(buff),"\n#Alignment\n");
	if(outformat=="clustal"){
		sprintf(buff+strlen(buff),"CLUSTALW\n\n");
		int c = len_aln / LINLEN;
		int d = len_aln % LINLEN;
		int offset1 = aln_start[0];
		int offset2 = aln_start[1];
		int start=0;
		for (int i=1;i<=c;i++){
			tmpstr1 = seq1_aln.substr(start,LINLEN);
			gapnum1 = count(tmpstr1.begin(), tmpstr1.end(), '-');
			tmpstr2 = seq2_aln.substr(start,LINLEN);
			gapnum2 = count(tmpstr2.begin(), tmpstr2.end(), '-');
			
			sprintf(buff+strlen(buff),"%-15s %s  %-5d\n",id1.c_str(),seq1_aln.substr(start,LINLEN).c_str(),offset1+LINLEN-1-gapnum1);
			sprintf(buff+strlen(buff),"%-15s %s  %-5d\n\n",id2.c_str(),seq2_aln.substr(start,LINLEN).c_str(),offset2+LINLEN-1-gapnum2);
			
			//sprintf(buff+strlen(buff),"%-15s %-5d %s  %-5d\n",id1.c_str(),offset1,seq1_aln.substr(start,LINLEN).c_str(),offset1+LINLEN-1-gapnum1);
			//sprintf(buff+strlen(buff),"%-20s  %s \n" ,"",quality.substr(start,LINLEN).c_str());
			//sprintf(buff+strlen(buff),"%-15s %-5d %s  %-5d\n\n",id2.c_str(),offset2,seq2_aln.substr(start,LINLEN).c_str(),offset2+LINLEN-1-gapnum2);
			
			start+=LINLEN;
			offset1+=LINLEN-gapnum1;
			offset2+=LINLEN-gapnum2;
		}
		if(d>0){
			tmpstr1 = seq1_aln.substr(start,d);
			gapnum1 = count(tmpstr1.begin(), tmpstr1.end(), '-');
			tmpstr2 = seq2_aln.substr(start,d);
			gapnum2 = count(tmpstr2.begin(), tmpstr2.end(), '-');
			sprintf(buff+strlen(buff),"%-15s %s  %-5d\n",id1.c_str(),seq1_aln.substr(start,d).c_str(),offset1+d-1-gapnum1);
			sprintf(buff+strlen(buff),"%-15s %s  %-5d\n\n",id2.c_str(),seq2_aln.substr(start,d).c_str(),offset2+d-1-gapnum2);
			//sprintf(buff+strlen(buff),"%-15s %-5d %s  %-5d\n",id1.c_str(),offset1,seq1_aln.substr(start,d).c_str(),offset1+d-1-gapnum1);
			//sprintf(buff+strlen(buff),"%-20s  %s \n" ,"",quality.substr(start,d).c_str());
			//sprintf(buff+strlen(buff),"%-15s %-5d %s  %-5d\n\n",id2.c_str(),offset2,seq2_aln.substr(start,d).c_str(),offset2+d-1-gapnum2);
			
		sprintf(buff+strlen(buff),"\n");
		}
	}
	else if(outformat=="fasta"){
		int c = len_aln / LINLEN;
		int d = len_aln % LINLEN;
		int start=0;
		sprintf(buff+strlen(buff),">%s\n",id1.c_str() );
		for (int i=1;i<=c;i++){
			sprintf(buff+strlen(buff),"%s\n",seq1_aln.substr(start,LINLEN).c_str());
			start+=LINLEN;
		}
		if(d>0){
			sprintf(buff+strlen(buff),"%s\n",seq1_aln.substr(start,LINLEN).c_str());
		}
		start=0;
		sprintf(buff+strlen(buff),">%s\n",id2.c_str() );
		for (int i=1;i<=c;i++){
			sprintf(buff+strlen(buff),"%s\n",seq2_aln.substr(start,LINLEN).c_str());
			start+=LINLEN;
		}
		if(d>0){
			sprintf(buff+strlen(buff),"%s\n",seq2_aln.substr(start,LINLEN).c_str());
		}
	}
	return buff;
}

