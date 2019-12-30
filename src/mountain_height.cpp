#include "mountain_height.h"


const int DANGLE=2;
double T= 37.;
double kT = (T + 273.15)*1.98717/1000.; /* kT in kcal/mol */
int engFlag=2004;//Turner2004; Turner1999 and Andronescu072007 also possible

void BasePairProbabilities(char * seq , int n, double * mfe, char* mfeStr, double *ensEng,char * argv0, double ** bppr)
{
	readEnergyFile(argv0);
	FLT_OR_DBL *bppm;
	dangles = DANGLE;
	temperature = T;
	*mfe = fold(seq,mfeStr);
	//printf("seq:%s,temp:%f,d:%d,mfe:%f\n",seq,T,dangleFlag,*mfe);
	pf_scale = exp((-*mfe/kT)/n);
	//update_pf_params(n);
	*ensEng = pf_fold(seq, NULL);
	bppm = export_bppm();

	//double **bppr = new double *[n];
	//for(int i = 0; i < n; i++)
	//	bppr[i] = new double[n+1];
	for(int i=0;i<n;i++)
		for(int j=0;j<n+1;j++)
			bppr[i][j]=0;
	double pr = 0;
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			bppr[i][j]=bppm[iindx[i+1] - (j+1)  ];
			bppr[j][i]=bppm[iindx[i+1] - (j+1)  ];
		}
	}
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			pr +=bppr[i][j];
		}
		bppr[i][n]= 1-pr;
		pr=0;
	}
	//print base pair probabilities
	//printf("\nbase pair probabilities:\n");
	//for (int i=0; i<n;i++)
	//{
		//for (int j=0;j<n+1;j++)
			//printf	("%.5f\t",bppr[i][j]);
		//printf("\n");
	//}
	free_arrays();
    free_pf_arrays();
    free_co_arrays();
    free_co_pf_arrays();
}


//Given h=[h_0,...,h_n-1] return h' where h'[0]=0 , h'[k]=h[k]-h[k-1] for 0<k<n
vector<double> heightDiff(vector<double> h){
	vector<double> h_diff;
	h_diff.push_back(h[0]);
	for(int i=1;i<h.size();++i){
		h_diff.push_back(h[i]-h[i-1]);
	}
	return h_diff;
}

vector<double> expectedModifiedHeight(double ** bppr, int n)
{
	vector<double> h;
	for(int i=0;i<n;i++)
		h.push_back(0);
	//--------base case--------
	for (int i=1;i<n;i++)
		h[0] += bppr[0][i]*(1/float(i));

	for(int k=0;k<n-1;k++)
	{
		//printf("K=%d\n",k);
		h[k+1] += h[k];
		for (int y=k+2;y<n;y++){
			h[k+1] += bppr[k+1][y]*(1/float((y-k-1)));
			//printf("+: %d,%d,%lf\n",k+1,y,bppr[k+1][y]);
		}
		for (int x=0; x<=k ; x++){
			h[k+1] -= bppr[x][k+1] * (1/float((k+1-x)));
			//printf("-: %d,%d,%lf\n",x,k+1,bppr[x][k+1]);
		}
	}
//	printf("\n****************expected modified height at each position***********************\n");
	//for(int i=0;i<n;i++)
		//h[i] = h[i];
	return h;
}

vector<double> expectedHeight(double **bppr , int n)
{
	vector<double> h;
	int i;
	for(i=0;i<n;i++)
		h.push_back(0);
	//--------base case--------
	for (i=1;i<n;i++)
		h[0] += bppr[0][i];

	for(i=0;i<n-1;i++)
	{
		//printf("i=%d\n",i);
		h[i+1] += h[i];
		for (int y=i+2;y<n;y++){
			h[i+1] += bppr[i+1][y];
			//printf("+: %d,%d,%lf\n",i+1,y,bppr[i+1][y]);
		}
		for (int x=0; x<=i ; x++){
			h[i+1] -= bppr[x][i+1];
			//printf("-: %d,%d,%lf\n",x,i+1,bppr[x][i+1]);
		}
	}
	//printf("\n****************expected height at each position***********************\n");
	//for(int i=0;i<n;i++)
			//printf("%d\t%lf\n",i,h[i]);
	return h;
}

double basePairProbability(double **bppr , int n)
{
	double bpprob=0.0;
	for(int i=0;i<n;i++){
		bpprob += 1.0 - bppr[i][n];
	}	
	return bpprob/n;
}

void readEnergyFile(char * argv0){
	//-------load the correct energy parameter file
	char * execPath = getExecPath(argv0);
	//printf("energy file path: %s\n",execPath);
	char turner99[100] = "/energy_params/rna_turner1999.par";
	char turner04[100] = "/energy_params/rna_turner2004.par";
	char andronescu07[100] = "/energy_params/rna_andronescu2007.par";

	if(engFlag==1999)
		strcat(execPath,turner99);
	else if (engFlag==2004)
		strcat(execPath , turner04);
	else if (engFlag==2007)
		strcat(execPath , andronescu07);
	read_parameter_file(execPath);
}

void printHeights(const vector<double> & h){	
	printf("[");
	int i;
	for(i=0;i<h.size()-1;i++)
		printf("%.2f,",h[i]);
	printf("%.2f]\n",h[i]);
}

string sprintHeights(const vector<double>& h){	
	char buff[50000]="";
	int i;
	sprintf(buff+strlen(buff),"[");
	for(i=0;i<h.size()-1;i++)
		sprintf(buff+strlen(buff),"%.2f,",h[i]);
	sprintf(buff+strlen(buff),"%.2f]\n",h[i]);

	return buff;
}

void runAlifold(const char ** seqs, char * str, double &mfe ){
	mfe = alifold( seqs, str );
	free_alifold_arrays();
}
