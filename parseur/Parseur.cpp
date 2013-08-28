#include"Parseur.hpp"

Parseur::Parseur(unsigned int argc, char* argv[]):
	locked(false)
{
	if(argv[argc-1][0] == '-'){
		std::cerr<<"warning : Parseur : if the last argument has no var and starts with a -, the prog will crash"<<std::endl;
	}
	init(argc-1, argv);
}

Parseur::Parseur(unsigned int argc, char* argv[], unsigned int N_default):
	locked(false)
{
	if( 2*N_default < argc ) { 
		if(argv[argc-1][0] == '-'){
			std::cerr<<"warning : Parseur : if the last argument has no var and starts with a -, the prog will crash"<<std::endl;
		}
		init(2*N_default,argv);

		for(unsigned int i(2*N_default+1);i<argc;i++){
			val.push_back(argv[i]);
		}
		//for(unsigned int i(0); i<val.size();i++){
			//std::cout<<val[i]<<std::endl;
		//}
	} else {
		locked=true;
			std::cerr<<"warning : Parseur : you should give "<< N_default<<" arguments"<<std::endl;
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<var.size();i++){
		if(!used[i]){ std::cerr<<"Parseur : variable "<<var[i]<<" was given as input but not used"<<std::endl;}
	}
}

void Parseur::init(unsigned int N, char* argv[]){
	unsigned int j(0);
	for(unsigned int i(0);i<N;i++){
		if(argv[i+1][0]=='-'){
			std::string tmp(argv[i+1]);
			var.push_back(tmp.substr(1));
			val.push_back(argv[i+2]);
			i++;
		} else {
			std::stringstream ss;
			ss<<j;
			var.push_back(ss.str());
			val.push_back(argv[i+1]);
			j++;
		}
		used.push_back(false);
	}
}
