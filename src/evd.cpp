#include "evd.h"

EVD::EVD(vector<double> dist, int C) {
	lambda=mu=sum_of_elements=num_iterations=0;
	data=dist;
	size = data.size();
	for(vector<double>::iterator it = data.begin(); it != data.end(); ++it)
		sum_of_elements += *it;
	CENSORED = C;
}

EVD::~EVD() {

}

void EVD::printHistogram(double min, double max, double bin_size){
	int total = 0,relative;
	for(int i=0; i<histogram.size();++i)
		total +=histogram[i].size();
	printf("#HISTOGRAM\nTOTAL:%d\n-",total);
	for(int i=0; i<histogram.size();++i){
		printf("%d\t%zu\t\t[%.2f,%.2f)\t|",i,histogram[i].size(),min+(i*bin_size),min+(i+1)*bin_size);
		relative = ceil(histogram[i].size()*100.0/total);
		for (int j=0;j<relative;++j){
			printf("*");
		}
		printf("\n");
	}
	return;
}	

void EVD::makeHistogram(){
	int i;
	const int NUM_BINS = 20;
	double min=9999999,max=-9999999;
	double bin_size=0;
	int bin_index;
	histogram.resize(NUM_BINS);
	for(vector<double>::iterator it=data.begin();it!=data.end();++it){
			if(max<*it){
				max=*it;}
			if(min>*it){
				min=*it;}
	}
	bin_size = (max - min)/NUM_BINS;
	for(vector<double>::iterator it=data.begin();it!=data.end();++it){
		if(fabs(*it-max) <0.0001)
			bin_index = NUM_BINS-1;
		else
			bin_index = floor(fabs(*it-min)/bin_size);
		histogram[bin_index].push_back(*it);
	}
	printHistogram(min,max,bin_size);
	return;
}

double EVD::pvalue(double x, double mu, double lambda){
	return exp(-exp(-lambda*(x-mu)));
}

int EVD::totalFit(vector<double> data){
	int i;
	int max_it_num=100;
	double thresh = 0.00001;
	double lam,x,expsum,xexpsum,x2expsum,e,f,fprime;
	int n =data.size();
	//Initialize lambda
	lam = 0.2;

	//update lambda until convergence or many iterations
	for(num_iterations=0;num_iterations<max_it_num;num_iterations++){
		expsum=xexpsum=x2expsum=0.0;
		//compute equation 10 and 11 in Eddy
		for (i = 0; i < n; ++i) {
			x=data[i];
			e = exp(-1.0*lam*x);
			expsum += e ;
			xexpsum += x*e;
			x2expsum += x*x*e;
		}
		f = (1.0/lam) - (sum_of_elements/(float)n) + (xexpsum/expsum);
		//printf("%f\n",f);
		if (fabs(f)<thresh)
			{
				//printf("lam = %f\n",lam);
				break;
			}
		fprime = (pow(xexpsum,2)/pow(expsum,2)) - (x2expsum/expsum) - (1.0/(lam*lam));
		lam = lam - f / fprime;
		if (lam <= 0.) lam = 0.001;
		//printf("%d,%f,%f\n",num_iterations,f,(1.0/size)*expsum);
	}
	if(num_iterations==max_it_num){
		printf("WARNING: max likelihood did not converge after %d iterations\n",max_it_num);
		return 0;
	}
	//plug lambda to equation 9 of Eddy and compute mu
	mu = -(1.0/lam) * log(expsum/(float)n);
	lambda = lam;
	return 1;
	
}


int EVD::censoredFit(vector<double> data, double c, int z){
	int i;
	int max_it_num=100;
	double thresh = 0.00001;
	double lam,x,expsum,xexpsum,x2expsum,e,f,fprime;
	double psx;
	double explamc ;
	
	//Initialize lambda
	lam = 0.2;

	
	//update lambda until convergence or many iterations
	for(num_iterations=0;num_iterations<max_it_num;num_iterations++){
		expsum=xexpsum=x2expsum=0.0;
		explamc = exp(-1.0*lambda*c);

		//compute equation 16 and 17 in Eddy
		for (i = 0; i < data.size(); ++i) {
			x=data[i];
			e = exp(-1.0*lam*x);
			expsum += e ;
			xexpsum += x*e;
			x2expsum += x*x*e;
		}
		f = (1.0/lam) - (sum_of_elements/size) + ((xexpsum+z*c*explamc)/(z*explamc+expsum));
		//printf("%f\n",f);
		if (fabs(f)<thresh)
			{
				//printf("lam = %f\n",lam);
				break;
			}
		fprime = pow(xexpsum + z*c*explamc,2)/pow(expsum+z*explamc,2) - (x2expsum+z*c*c*explamc)/(expsum+z*explamc) - 1.0/(lam*lam);
		lam = lam - f / fprime;
		if (lam <= 0.) lam = 0.001;
		//printf("%d,%f,%f\n",num_iterations,lam,f / fprime);
	}
	if(num_iterations==max_it_num){
		printf("max likelihood did not converge after %d number of iterations\n",max_it_num);
		return 0;
	}
	//plug lambda to equation 9 of Eddy and compute mu
	mu = -(1.0/lam) * log((expsum+z*explamc)/data.size());
	lambda = lam;
	//printf("%d\n",data.size());

	return 1;
}


//Based on maximum likelihood method from Mott(1992) and Lawless(1982).
//Adapted from
//     Maximum likelihood fitting of extreme value distribution, S. Eddy(1997)
int EVD::fit(){
	
	if(CENSORED==1){
		const int max_it_num=100;
		double pval;
		double lowbound=0.,highbound=1.;
		int i,j,num_iterations,peak;
		vector<double> censored_data, new_censored_data;
		int prev_count=0;
		double prev_mu=-99999,prev_lam=-99999;
		const float TH= 0.0001;
		
		//printf("new highbound:%f,new lowbound:%f\n",highbound,lowbound);
		censored_data = data;
		totalFit(censored_data);
		//printf("mu=%f,lambda=%f,\n",mu,lambda);
		//printf("*************\n");
		
		for(num_iterations=1;num_iterations<=max_it_num;num_iterations++){
			new_censored_data.clear();
			//construct the censored vector
			lowbound = (double) (1./(censored_data.size()+1.));
			highbound = (double) (censored_data.size()/(censored_data.size()+1.));
			//printf("******CENSORED******\niteration=%d\n",num_iterations);
			//printf("mu=%f,lambda=%f,\n",mu,lambda);
			//printf("new highbound:%f,new lowbound:%f\n",highbound,lowbound);
			//printf("data_size:%d\n",censored_data.size());
			for(i=0;i<censored_data.size();++i){
				pval = pvalue(censored_data[i],mu,lambda);
				if(pval <= highbound && pval >= lowbound)
					new_censored_data.push_back(censored_data[i]);
			}
			//printf("new_data_size:%d\n",new_censored_data.size());
				
			
			if((new_censored_data.size()-censored_data.size())==0 || (fabs(lambda-prev_lam)<TH && fabs(mu-prev_mu)<TH) ) break;
			
			//save the parameters before fitting the new data
			prev_lam = lambda;
			prev_mu=mu;
			prev_count = censored_data.size();
			
			totalFit(new_censored_data);
			
			}

	}
	else{
		totalFit(data);
	}
	return 1;
}

