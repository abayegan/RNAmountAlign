#include "aux.h"

double round_percision(double x,int percision){
	return round(pow(10,percision)*x)/pow(10,percision);
}

double readBasePairProb(string fname,double pa,double pc,double pg,double pu){
	ifstream infile;
	infile.open(fname.c_str());
	double a,c,g,u,op;
	if (!infile) reportErr("file %s could not be found!");
	while (infile >> a >> c >> g >> u >> op ) {
		if(equal(a,pa) && equal(c,pc) && equal(g,pg) && equal(u,pu))
			break;
	}
	return op*2.0;
}

bool equal(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

void writeVec(char * fname, vector<double> vec){
	ofstream f(fname);
    for(vector<double>::iterator it = vec.begin(); it != vec.end(); ++it)
		f << *it << '\n';
}

float getGCcontent(string seq){
	float gccont=0.0;
	for (string::iterator it=seq.begin();it!=seq.end();it++){
		switch(*it){
			case 'C':
			case 'G':
			case 'c':
			case 'g':
				gccont+=1;
				break;
		}
	}
	//printf("seq:%s ,gc:%f\n",seq.c_str(),gccont/seq.size());
	return (gccont/seq.size());
}

vector<string> mergeVectors(vector<string> v1,vector<string> v2){
	vector<string>v;
	v.insert( v.end(), v1.begin(), v1.end() );
	v.insert( v.end(), v2.begin(), v2.end() );
	return v;
}

double ** allocate2DdoubleArray(int row, int col){
	double ** array = new double * [row];
	for(int i=0; i<row;i++){
		array[i] = new double[col];
	}
	return array;
}

char ** allocate2DcharArray(int row, int col){
	char ** array = new char * [row];
	for(int i=0; i<row;i++){
		array[i] = new char[col];
	}
	return array;
}

void free2DdoubleArray(double ** array,int row){
	for(int i = 0; i < row; ++i) {
	    delete [] array[i];
	}
	delete [] array;
}

void free2DcharArray(char ** array,int row){
	for(int i = 0; i < row; ++i) {
	    delete [] array[i];
	}
	delete [] array;
}



//NOTE: order of the parameters are important in order to
//return the correct pointer.
double max( double p, double q, double r, char * ptr )
{
	const int INF = numeric_limits<int>::max();
	double  max=-INF;
	if(p>max){
		max = p ;
		*ptr = 'P' ;
	}
	if(q>max){
		max = q ;
		*ptr = 'Q' ;
	}
	if (r>=max){
		max = r ;
		*ptr = 'R' ;
	}
//	printf("%f\t%f\t%f\t%f\t%c\n",p,q,r,max,*ptr);
	return max;
}

void print2Dvec(vector<vector<double> > vec){
	printf("%5c",' ');
	for (int j = 0; j < vec[0].size(); ++j)
		printf("%5d",j);
	printf("\n");
	for (int i = 0; i < vec.size(); ++i){
		printf("%5d",i);
		for (int j = 0; j < vec[i].size(); ++j){
				printf("%5.2f",vec[i][j]);
		}
    printf("\n");
	}
}


void printMatrix(vector<vector<char> > mat,int n,int m){
	printf("%5c",' ');
	for (int j = 0; j < m; ++j)
		printf("%5d",j);
	printf("\n");
	for(int i=0;i<n;i++){
		printf("%5d",i);
		for(int j=0;j<m;j++){
			printf("%5c",mat[i][j]);
		}
		printf("\n");
	}
}

void printMatrix1(double ** mat,int n,int m){
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++)
			printf("%.1f ",mat[i][j]);
		cout<<endl;
	}
}

double getPercentileFromSortedVector(vector<double> vec,int percentile){
	int idx = ceil((percentile/100.0) * vec.size());
	return vec.at(max(idx-1,0));
}

void reportErr(const char * msg)
{
    cerr << "ERROR: " << msg << "\n";
    exit(EXIT_FAILURE);
}

void reportWarning(const char * msg)
{
    cerr << "WARNING: " << msg << "\n";
}

int char2num(char c){
	int num;
	if (c=='A' || c=='a' )
		num=0;
	else if(c=='C' || c=='c')
		num=1;
	else if(c=='G' || c=='g')
		num=2;
	else if(c=='U' || c=='u' || c=='T' || c=='t')
		num=3;
	else
		reportErr("input sequences should only contain A,C,G,U !");
	return num;
}

void getNucleotideProbabilities(string seq,vector<double> & nucpr){
	float suma,sumc,sumg,sumu;
	suma=sumc=sumg=sumu=0.0;
	for(string::iterator it=seq.begin(); it!=seq.end(); ++it){
		switch (*it){
			case 'A':
			case 'a':
				suma+=1;
				break;
			case 'C':
			case 'c':
				sumc+=1;
				break;
			case 'G':
			case 'g':
				sumg+=1;
				break;
			case 'U':
			case 'T':
			case 't':
			case 'u':
				sumu+=1;
				break;
			}
	}
	nucpr.push_back(suma/seq.size());nucpr.push_back(sumc/seq.size());
	nucpr.push_back(sumg/seq.size());nucpr.push_back(sumu/seq.size());
	//printf("Nucleotide probabilities:%f\t%f\t%f\t%f\n",nucpr[0],nucpr[1],nucpr[2],nucpr[3]);
}

void generateRandomRNA(char * seq, int length,float pa,float pc,float pg,float pu){
	string NUC = "ACGU";
	float p[4] = {pa,pc,pg,pu};
	for(int i=0;i<length;i++){
		seq[i]= NUC[rouletteWheel(p,4)];
	}
	seq[length]='\0';
	//printf("random RNA:%s \n",seq);
}

int rouletteWheel(float *p, int n){
	float rnd;
	float sum=0.0;
	int i;
	//rnd = sre_random();
	rnd = rand() / double(RAND_MAX);
	//printf("%f\n",rnd);
	for (i = 0; i < n; i++)
	{
	  sum += p[i];
	  if (rnd < sum) return i;
	}
	return (int) (sre_random() * n);
}

//taken from squid library by S. Eddy
static int sre_randseed = 42;	/* default seed for sre_random()   */
double sre_random(){
  static long  rnd1;		/* random number from LCG1 */
  static long  rnd2;            /* random number from LCG2 */
  static long  rnd;             /* random number we return */
  static long  tbl[64];		/* table for Bays/Durham shuffle */
  long x,y;
  int i;

  /* Magic numbers a1,m1, a2,m2 from L'Ecuyer, for 2 LCGs.
   * q,r derive from them (q=m/a, r=m%a) and are needed for Schrage's algorithm.
   */
  long a1 = 40014;
  long m1 = 2147483563;
  long q1 = 53668;
  long r1 = 12211;

  long a2 = 40692;
  long m2 = 2147483399;
  long q2 = 52774;
  long r2 = 3791;

  if (sre_randseed > 0)
    {
      rnd1 = sre_randseed;
      rnd2 = sre_randseed;
				/* Fill the table for Bays/Durham */
      for (i = 0; i < 64; i++) {
	x    = a1*(rnd1%q1);   /* LCG1 in action... */
	y    = r1*(rnd1/q1);
	rnd1 = x-y;
	if (rnd1 < 0) rnd1 += m1;

	x    = a2*(rnd2%q2);   /* LCG2 in action... */
	y    = r2*(rnd2/q2);
	rnd2 = x-y;
	if (rnd2 < 0) rnd2 += m2;

	tbl[i] = rnd1-rnd2;
	if (tbl[i] < 0) tbl[i] += m1;
      }
      sre_randseed = 0;		/* drop the flag. */
    }/* end of initialization*/


  x    = a1*(rnd1%q1);   /* LCG1 in action... */
  y    = r1*(rnd1/q1);
  rnd1 = x-y;
  if (rnd1 < 0) rnd1 += m1;

  x    = a2*(rnd2%q2);   /* LCG2 in action... */
  y    = r2*(rnd2/q2);
  rnd2 = x-y;
  if (rnd2 < 0) rnd2 += m2;

   			/* Choose our random number from the table... */
  i   = (int) (((double) rnd / (double) m1) * 64.);
  rnd = tbl[i];
			/* and replace with a new number by L'Ecuyer. */
  tbl[i] = rnd1-rnd2;
  if (tbl[i] < 0) tbl[i] += m1;

  return ((double) rnd / (double) m1);
}


char* getExecPath(char* argv0){
	char* exec_path=(char*) malloc (sizeof(char)*PATH_MAX);
	exec_path[0]='.';
	exec_path[1]='\0';
	int path_length=1;
	FILE* file;
	#if defined(_WIN32)
	/* GetModuleFileName */
	/* _pgmptr */
	#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix)
	char buff[PATH_MAX];
	int bufsize = PATH_MAX-1;
	if(file = fopen("/proc/self/exe", "r")){
		fclose(file);
		ssize_t len = readlink("/proc/self/exe", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}
	}
	else if(file = fopen("/proc/curproc/file", "r")){
		fclose(file);
		ssize_t len = readlink("/proc/curproc/file", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}
	}
	else if(file = fopen("/proc/self/path/a.out", "r")){
		fclose(file);
		ssize_t len = readlink("/proc/self/path/a.out", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}
	}
	else{
		strcpy(exec_path, argv0);
	}
	char* slash_pointer= strrchr(exec_path,'/');
	int slash_position = (int)(slash_pointer - exec_path);
	if(slash_position != path_length){
		exec_path[slash_position]='\0';
	}
	else{
		exec_path[0]='.';
		exec_path[1]='\0';
	}

	#endif


	return exec_path;
}

