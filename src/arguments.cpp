#include "arguments.h"


Arguments::Arguments(){
	gap_init = -3;
	gap_ext = -1;
	gap_init_str = -1000;
	gap_ext_str = -1000;
	str_weight = 0.5;
	percentile = 0;
	print_prob = false;
	normal = false;
	eps_name="";
	fasta_name="";
	query_name="";
	target_name="";
	output_name="";
	outformat ="clustal";
	verbose=false;
	window_size = 300;
	step_size = 200;
	rnd_num = 500;
	gc_seg = 0.05;
	fixed_gamma=true;
	modified_height = false;
	stats=false;
	evd=false;
	matrix="matrices/RIBOSUM85-60.mat";
	alifold=false;
	qsearch=false;
}
void Arguments::parseArgs(const int argc, char * argv[]){
		if(!(strstr(argv[0],"RNAmountAlignScan"))){
			aln_type='g';
			qsearch=false;
		}
		else{
			aln_type='q';
			qsearch=true;
		}
		argv0=argv[0];
		char * execPath = getExecPath(argv0);
		strcat(execPath,"/");
		strcat(execPath,matrix);
		matrix=execPath;
		
		if(argc<3){
			printf("error in the input argument!\n");
			printUsage(argv[0]);
		}
		for(int i=1;i<argc;i++){
			if(strcmp(argv[i],"-s")==0){
				if(i+2<argc && argv[i+1][0]!='-' && argv[i+2][0]!='-'){
					seqs.push_back(argv[i+1]);
					seqs.push_back(argv[i+2]);
					ids.push_back("seq1");
					ids.push_back("seq2");
				}
				else{
					printf("error in -s argument!\n");
					printUsage(argv[0]);
				}
			}
			if(strcmp(argv[i],"-f")==0){
				if(i+1<argc && argv[i+1][0]!='-'){
					fasta_name = argv[i+1];
				}
				else{
					printf("error in -f argument!\n");
					printUsage(argv[0]);
				}
			}
			if(strcmp(argv[i],"-qf")==0){
				if(i+1<argc && argv[i+1][0]!='-'){
					query_name = argv[i+1];
				}
				else{
					printf("error in -qf argument!\n");
					printUsage(argv[0]);
				}
			}
			if(strcmp(argv[i],"-tf")==0){
				if(i+1<argc && argv[i+1][0]!='-'){
					target_name = argv[i+1];
				}
				else{
					printf("error in -tf argument!\n");
					printUsage(argv[0]);
				}
			}
			
			else if(strcmp(argv[i],"-o")==0){
				if(i+1<argc && argv[i+1][0]!='-' ){
					output_name = argv[i+1];
				}
				else{
					printf("error in -o argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-gi")==0){
				if(i+1<argc){
					gap_init = strtof(argv[i+1],NULL);
					gap_init_str = gap_init;
				}
				else{
					printf("error in -gi argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-ge")==0){
				if(i+1<argc){
					gap_ext = strtof(argv[i+1],NULL);
					gap_ext_str = gap_ext;
				}
				else{
					printf("error in -ge argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-c")==0){
				if(i+1<argc){
					gap_init_str = strtof(argv[i+1],NULL);
				}
				else{
					printf("error in -c argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-d")==0){
				if(i+1<argc){
					gap_ext_str = strtof(argv[i+1],NULL);
				}
				else{
					printf("error in -d argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-gamma")==0){
				if(i+1<argc && argv[i+1][0]!='-' ){
					str_weight = strtof(argv[i+1],NULL);
				}
				else{
					printf("error in -gamma argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-dyn")==0){
					fixed_gamma = false;
			}
			else if(strcmp(argv[i],"-mod")==0){
					modified_height = true;
			}
			else if(strcmp(argv[i],"-t")==0){
				if(i+1<argc && argv[i+1][0]!='-' ){
					percentile = atoi(argv[i+1]);
				}
				else{
					printf("error in -t argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-p")==0){
				print_prob = true;
			}
			else if(strcmp(argv[i],"-n")==0){
				normal = true;
			}
			else if(strcmp(argv[i],"-m")==0){
				if(i+1<argc && argv[i+1][0]!='-'){
					matrix= argv[i+1];
				}
				else{
					printf("error in -m argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-format")==0){
				if(i+1<argc && argv[i+1][0]!='-'){
					outformat= argv[i+1];
				}
				else{
					printf("error in -format argument!\n");
					printUsage(argv[0]);
				}
				if(outformat!="clustal" && outformat!="fasta"){
					printf("error in -format argument!\n");
					printUsage(argv[0]);
				}
			}
			else if(strcmp(argv[i],"-semi")==0){
				aln_type = 'q';
			}
			else if(strcmp(argv[i],"-global")==0){
				aln_type = 'g';
			}
			else if(strcmp(argv[i],"-local")==0){
				aln_type = 'l';
			}
			else if(strcmp(argv[i],"-alifold")==0){
				alifold = true;
			}
			else if(strcmp(argv[i],"-window")==0){
				window_size = atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-step")==0){
				step_size = atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-num")==0){
				rnd_num = atoi(argv[i+1]);
			}
			else if(strcmp(argv[i],"-gc")==0){
				gc_seg = atoi(argv[i+1])/100.0;
			}
			//else if(strcmp(argv[i],"-eps")==0){
				//if(i+1<argc && argv[i+1][0]!='-'){
					//eps_name= argv[i+1];
				//}
				//else{
					//printf("error in -eps argument!\n");
					//printUsage(argv[0]);
				//}
			//}
			else if(strcmp(argv[i],"-stat")==0){
				stats=true;
			}
			else if(strcmp(argv[i],"-evd")==0){
				evd=true;
			}
			else if(strcmp(argv[i],"-v")==0){
				verbose = true;
			}
			else if(strcmp(argv[i],"-h")==0){
				printUsage(argv[0]);
			}
		}
	if(!(strstr(argv[0],"mountAlign") || strstr(argv[0],"Scan"))){
		reportErr("Please run executables RNAmountAlign or RNAmountAlignScan.\n");
		printUsage(argv[0]);
	}
	if((strstr(argv[0],"mountAlign")) && fasta_name!=""){
		readFasta(fasta_name,seqs,ids);
	}
	else if(strstr(argv[0],"RNAmountAlignScan")){
		if(aln_type=='g'){
			reportWarning("you are performing query search using global alignment.\n");
		}
		else if(aln_type=='l'){
			reportWarning("you are performing query search using local alignment.\n");
		}
		if(target_name!="" && query_name!=""){
			readFasta(query_name,seqs,ids);
			readFasta(target_name,seqs,ids);
		}
		else{
			reportErr("For query search define input fasta files using -qf and -tf flags.\n\n");
			printUsage(argv[0]);
		}
		if(seqs.size()>2) 
			reportErr("Each file must contain a single sequence\n");
	}
	loadRibosumMatrix(matrix);
	if(gap_init_str==-1000){ //the same structure gap penalties as sequence
		gap_init_str = gap_init;
		gap_ext_str = gap_ext;
	}

	if(seqs.size()<2){
		reportErr("Input sequences not defined\n");
	}
	if(aln_type=='q' && seqs[0].length()>seqs[1].length()){
		reportErr("In semi-global query(1st seq) must be shorter than target(2nd seq)\n");
	}
	
}

void Arguments::printUsage(const char * argv0){
	if(!(strstr(argv0,"RNAmountAlignScan"))){
		printf("\nUsage: %s [options]\n",argv0);
		printf("REQUIRED:\n"
		"\tRNAmountAlign:\n"
		"\t\t-f <string>	the input fasta file containing two sequences \n"
		"\t\tOR\n"
		"\t\t-s <string> <string>	provide sequence1 and sequence2\n\n"
		
		"\nOPTIONS:\n"
		
		"\nALIGNMENT:\n"
		"\t-gi <float>	Gap initiation penalty for sequence alignment(Default:-3).\n"
		"\t\t\tA negative value should be provided\n\n"
		
		"\t-ge <float>	Gap extension penalty for sequence alignment(Default:-1).\n"
		"\t\t\tA negative value should be provided\n\n"
		
		"\t-gamma <float>	Weight of the structural homology. Must be in [0,1](Default:0.5)\n"
		"\t\t\tsimilarity = gamma*str_sim + (1-gamma)*seq_sim\n\n"
		
		"\t-m <string>	Similarity matrix file in RIBOSUM format.(Default:RIBOSUM85-60.mat)\n"
		"\t\t\tAll RIBOSUM files are included in ./matrices directory\n\n"
		
		"\t-semi 		Perform semi-global alignment.\n" 
		"\t\t\tBoth gap ends of the first sequence will be free of penalty.\n\n"
		
		"\t-local 		Perform local alignment\n\n"
		
		"\t-global 		Perform global alignment.(Default alignment type)\n\n"
		
		"\t-alifold		Output the consensus structure of the alignment from RNAalifold\n\n"
		
		"\nSTATISTICS:\n"
		"\t\tKA=Karlin-Altschul; EVD=extreme value distribution, ND=normal distribution\n\n"
		"\t-stat		Report statistics based on the alignment type(Default: off). \n"
		"\t\t\tFor local alignments: E-value from Karlin-Altschul or EVD fitting\n"
		"\t\t\tFor global alignments: p-value from ND fitting\n"
		"\t\t\tFor semi-global alignments: p-value from ND fitting\n\n"
		
		"\t-evd		EVD fitting will be used for local alignments instead of KA\n\n"
		
		"\t-num <int>	Number of random sequences generated for fitting.(Default:500)\n\n"
		
		"\t-gc <int>	Size of GC bins (an integer between [0-100]). (Default: 5)\n"
		"\t\t\tThis is used only with with parameter fitting and not KA.\n\n"
		
		"\nOUTPUT:\n"
		"\t-o <string>	Write the output to a file.\n"
		"\t\t\tIf not used the output will be printed to stdout.\n\n"
		
		"\t-format <clustal|fasta>	Format of the alignment output. (Default: clustal)\n\n"
		
		"\t-v		Verbose output. Prints MFE structures, the ensemble expected and incremental heights.\n\n"
		"\t-h		Print help\n\n");
	}
	else{
		printf("\nUsage: %s [options]\n",argv0);
		printf("REQUIRED:\n"
		"\tRNAmountAlignScan:\n"
		"\t\t-qf <string>	fasta file containing query sequence\n"
		"\t\tAND\n"
		"\t\t-tf <string>	fasta file containing target sequence\n\n"
		
		"\nOPTIONS:\n"
		
		"\nSLIDING\n"
		"\t -window <int>	Size of the sliding window(Default:300).\n\n"
		
		"\t -step <int>	Step size for incrementing the window start(Default:200).\n\n"
		
		
		"\nALIGNMENT:\n"
		"\t-gi <float>	Gap initiation penalty for sequence alignment(Default:-3).\n"
		"\t\t\tA negative value should be provided\n\n"
		
		"\t-ge <float>	Gap extension penalty for sequence alignment(Default:-1).\n"
		"\t\t\tA negative value should be provided\n\n"
		
		"\t-gamma <float>	Weight of the structural homology. Must be in [0,1](Default:0.5)\n"
		"\t\t\tsimilarity = gamma*str_sim + (1-gamma)*seq_sim\n\n"
		
		"\t-m <string>	Similarity matrix file in RIBOSUM format.(Default:RIBOSUM85-60.mat)\n"
		"\t\t\tAll RIBOSUM files are included in ./matrices directory\n\n"
		
		"\t-semi 		Perform semi-global alignment.\n" 
		"\t\t\tBoth gap ends of the first sequence will be free of penalty.\n\n"
		
		"\t-local 		Perform local alignment\n\n"
		
		"\t-global 		Perform global alignment.(Default alignment type)\n\n"
		
		
		"\nSTATISTICS:\n"
		"\t\tKA=Karlin-Altschul; EVD=extreme value distribution, ND=normal distribution\n\n"
		"\t-stat		Report statistics based on the alignment type(Default: off). \n"
		"\t\t\tFor local alignments: E-value from Karlin-Altschul(default) or EVD fitting\n"
		"\t\t\tFor global alignments: p-value from ND fitting\n"
		"\t\t\tFor semi-global alignments: p-value from ND fitting\n\n"
		
		"\t-evd		EVD fitting will be used for local alignments instead of KA\n\n"
		
		"\t-num <int>	Number of random sequences generated for fitting.(Default:500)\n\n"
		
		"\t-gc <int>	Size of GC bins (an integer between [0-100]). (Default: 5)\n"
		"\t\t\tThis is used only with with parameter fitting and not KA.\n\n"
		
		"\nOUTPUT:\n"
		"\t-o <string>	Write the output to a file.\n"
		"\t\t\tIf not used the output will be printed to stdout.\n\n"
		
		"\t-format <clustal|fasta>	Format of the alignment output. (Default: clustal)\n\n"
		
		"\t-v		Verbose output. Prints MFE structures, the ensemble expected and incremental heights.\n\n"
		"\t-h		Print help\n\n");
	}
	
	exit(1);
}

void Arguments::printArgs(){
	string matt(matrix);
	printf("#PARAMETERS:\ngap init:%.1f;gap ext:%.1f;structure weight:%.2f\nseq sim matrix=%s\n",
			gap_init,gap_ext,str_weight,matt.substr(matt.find_last_of("/\\") + 1).c_str());
	string type="";
	if(aln_type=='l')
		type="local";
	else if(aln_type=='g')
		type="global";
	else if(aln_type=='q')
		type="semi-global";
	printf("Type = %s\n\n",type.c_str());
}

string Arguments::sprintArgs(){
	char txt[1000]={0};
	string matt(matrix);
	sprintf(txt+strlen(txt),"#PARAMETERS:\ngap init:%.1f; gap ext:%.1f; structure weight:%.2f\nseq sim matrix: %s\n",
			gap_init,gap_ext,str_weight,matt.substr(matt.find_last_of("/\\") + 1).c_str());
	string type="";
	if(aln_type=='l')
		type="local";
	else if(aln_type=='g')
		type="global";
	else if(aln_type=='q')
		type="semi-global";
	sprintf(txt+strlen(txt),"Type = %s\n\n",type.c_str());
	return txt;
}

void Arguments::readFasta(const string fname, vector<string> & seqs, vector<string> & ids){
	    ifstream input(fname.c_str());
	    if(!input.good()){
	    	char err[500];
	    	sprintf(err,"Could not open file %s\n", fname.c_str());
	    	reportErr(err);
	    }
	    vector<string> strs;
	    string line, name, content,str;
	    while(getline( input, line ).good() ){
	        if( line.empty() || line[0] == '>' ){
	            if( !name.empty() ){
	                seqs.push_back(content);
	                if (!str.empty())
						strs.push_back(str);
					else
						strs.push_back("MFE");
	                name.clear();
	            }
	            if( !line.empty() ){
	                name = line.substr(1);
	                ids.push_back(name);
	            }
	            content.clear();
	            str.clear();
	        } else if( !name.empty() ){
	            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
	                name.clear();
	                content.clear();
	            }
	            else if (line[0]=='(' || line[0]==')' || line[0]=='.' )
					str += line;
				else
	                content += line;
	            }
	        }
	    if( !name.empty() ){ // Print out what we read from the last entry
	       seqs.push_back(content);
	       if (!str.empty())
				strs.push_back(str);
			else
				strs.push_back("MFE");
	    }
	    if (seqs.size() >= MAXSEQNUM){
	    			char err[500];
	    			sprintf(err,"too much sequences in the fasta file. Maximum is %d\n", MAXSEQNUM);
	    			reportErr(err);
	    		}
		return;
}
void Arguments::loadRibosumMatrix(const char * fname){
	string line;
	int i,linenum=1;
	vector<string> words;
	vector<string> nucs;
	ifstream file (fname);
	if (file.is_open())
	{
		while ( getline (file,line) )
		{
			istringstream iss(line);
			if(linenum == 3){//A C G U
				copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(nucs));
			}
			else if(3<linenum && linenum<8){
				copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(words));
				for(i=1 ; i<words.size(); ++i){
					subMat[char2num(words[0][0])][char2num(nucs[i-1][0])] = atof(words[i].c_str());
				}
				words.clear();	
			}	
			linenum++;
		}
	file.close();
	}
	else reportErr("Unable to open RIBOSUM matrix!");
	for(i=0;i<=3;i++)
		for(int j=i+1;j<=3;j++)
			subMat[i][j] = subMat[j][i];
	//for(int i=0;i<=3;i++){
					//for(int j=0;j<=3;j++)
						//printf("%.2f\t",subMat[i][j]);
						//cout<<endl;
	//}	
}
