/*
 * Profile.cpp
 *
 *  Created on: Dec 18, 2016
 *      Author: Amir Bayegan
 */

#include "profile.h"

Profile::Profile(const vector<string> &pr1,
		const vector<string> &i1,
		const vector<vector<double> > &hh1, vector< vector <int> > coordinates, double b) {
	p = pr1;
	n = checkLength();
	m = p.size();
	ids= i1;
	h = hh1;
	computePSSM();
	h_sum = computeSumOfHeights();
	positions = coordinates;
	bp=b;

	//checkPositions();
}

Profile::Profile() {

}
Profile::~Profile() {

}
Profile::Profile(const Profile & that){
	p = that.p;
	n=that.n;
	m = that.m;
	ids = that.ids;
	h=that.h;
	pssm=that.pssm;
	h_sum = that.h_sum;
	positions = that.positions;
	bp = that.bp;
}

//Profile::Profile& operator=(const Profile &){
//
//}

void Profile::setLength(int nn){
	n=nn;
}

int Profile::getLength(){
	return n;
}

double Profile::getBp(){
	return bp;
}

int Profile::getNumberOfSeqs(){
	return m;
}

vector<string> Profile::getSeqs(){
	return p;
}

vector<vector<double> > Profile::getHeights(){
	return h;
}

vector<string> Profile::getIds(){
	return ids;
}

vector< vector<double> > Profile::getPSSM(){
	return pssm;
}

void Profile::checkPositions(){
	int start,stop;
	for (int i = 0; i <m ; ++i) {
		start = positions[i][0];
		stop = positions[i][1];
		if((stop-start)!=n){
			reportErr("Error in Profile::checkPositions()!\n");
		}
	}
}
/*Sequences in the profile must have the same length. Return length of the alignment.  */
int Profile::checkLength(){
	int n=p[0].length();
	for (int i = 1; i < p.size(); ++i) {
		if (p[i].length()!=n) {
			reportErr("Error in Profile::checkLength()! Sequences in the profile must have the same length.");
		}
	}
	return n;
}
vector<vector<int> > Profile::getPositions(){
	return positions;
}
void Profile::computePSSM(){
	int n = p[0].length();
	int m = p.size();
	float suma,sumc,sumg,sumu,sumgap;
	
	vector<double> a,c,g,u,gap;
	//printProfile("fasta");
	
	for (int j = 0; j < n; ++j) {
		suma=0,sumc=0,sumg=0,sumu=0,sumgap=0;
		for (int i = 0; i < m; ++i) {
			if (p[i][j]=='A') {
				suma+=1;
			} else if (p[i][j]=='C'){
				sumc+=1;
			} else if (p[i][j]=='G'){
				sumg+=1;
			} else if (p[i][j]=='U'){
				sumu+=1;
			}else if (p[i][j]=='-'){
				sumgap+=1;
			}
		}
		a.push_back(suma/m);
		c.push_back(sumc/m);
		g.push_back(sumg/m);
		u.push_back(sumu/m);
		gap.push_back(sumgap/m);
	}
	pssm.push_back(a);
	pssm.push_back(c);
	pssm.push_back(g);
	pssm.push_back(u);
	pssm.push_back(gap);
	//print2Dvec(pssm);
}

vector<double> Profile::getSumOfHeights(){
	return h_sum;
}
vector<double> Profile::computeSumOfHeights(){
	int n1 = h[0].size();
	int m1 = h.size();
	float sum;
	vector<double> h_sum(n,0);
	
	if (n!=n1 || m!=m1){
		reportErr("error in the heights of profiles!\n");
	}
	
	for (int j = 0; j < n; ++j) {
		sum=0;
		for (int i = 0; i < m; ++i) {
			sum+=h[i][j];
		}
		h_sum[j] =sum/m;
	}
	return h_sum;
}

void Profile::reverseProfile(){
	/*positions remains untouched */
	int i,j;
	for ( i = 0; i < m; ++i) {
		reverse(p[i].begin(), p[i].end());
		reverse(h[i].begin(), h[i].end());
	}
	for (i = 0; i < 5; ++i) {
		reverse(pssm[i].begin(),pssm[i].end());
	}
	reverse(h_sum.begin(),h_sum.end());
}


void Profile::printProfile(string format){
	char tmp1[2000];
	int gapnum;
	string tmpstr;
	int i,j;
	vector <int> offset(m);
	for(i=0;i<m;i++){
		offset[i]=positions[i][0]+1;
	}
	//sprintf(buff,"");
	if(verbose){
		sprintf(outtxt+strlen(outtxt),"#Expected heights:\n");
		for(i=0;i<m;i++){
			string hh1="";
			for(j=0;j<n;j++){
				sprintf(tmp1,"%.2f,",h[i][j]);
				hh1+=tmp1;
			}
			sprintf(outtxt+strlen(outtxt),"#%s:[%s]\n",ids[i].c_str(),hh1.c_str());
		}
	}
	sprintf(outtxt+strlen(outtxt),"\n#Alignment\n");
	sprintf(outtxt+strlen(outtxt),"CLUSTALW\n\n");
	if (format=="clustal")
	{
		int c = n / LINLEN;
		int d = n % LINLEN;
		int start=0;
		for(j=0;j<c;j++){
			for(i=0;i<m;i++){
				tmpstr=p[i].substr(start,LINLEN);
				gapnum=count(tmpstr.begin(),tmpstr.end(),'-');
				//printf("%-20s %-5d %s  %-5d\n",ids[i].c_str(),offset[i],p[i].substr(start,LINLEN).c_str(),offset[i]+LINLEN-1-gapnum);
				sprintf(outtxt+strlen(outtxt),"%-15s %s  %-5d\n",ids[i].c_str(),p[i].substr(start,LINLEN).c_str(),offset[i]+LINLEN-1-gapnum);
				offset[i]+=LINLEN - gapnum;
			}
			start+=LINLEN;
			sprintf(outtxt+strlen(outtxt),"\n\n");
		}
		if(d>0){
			for(i=0;i<m;i++){
				tmpstr=p[i].substr(start,LINLEN);
				gapnum=count(tmpstr.begin(),tmpstr.end(),'-');
				//printf("%-15s %-5d %s  %-5d\n",ids[i].c_str(),offset[i],p[i].substr(start,d).c_str(),offset[i]+d-1-gapnum);
				sprintf(outtxt+strlen(outtxt),"%-15s %s  %-5d\n",ids[i].c_str(),p[i].substr(start,d).c_str(),offset[i]+d-1-gapnum);
			}
		}
		sprintf(outtxt+strlen(outtxt),"\n");
	}
	else if(format=="fasta"){
		int c = n / LINLEN;
		int d = n % LINLEN;
		int j;
		for(i=0;i<m;i++){
			sprintf(outtxt+strlen(outtxt),">%s\n",ids[i].c_str() );
			j=c;
			int offset1=0;
			int start=0;
			while(j>0){
				sprintf(outtxt+strlen(outtxt),"%s\n",p[i].substr(start,LINLEN).c_str());
				start+=LINLEN;
				j--;
			}
			if(d>0){
			sprintf(outtxt+strlen(outtxt),"%s\n",p[i].substr(start,d).c_str());
			}
		}
	}
}
