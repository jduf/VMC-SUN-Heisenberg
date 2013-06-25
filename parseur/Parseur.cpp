#include"Parseur.hpp"

Parseur::Parseur(unsigned int argc, char* argv[]):
	locked(false)
{
	std::cerr<<"Parseur : if the last argument has no var and starts with a -, the prog will crash"<<std::endl;
	unsigned int j(0);
	for(unsigned int i(1);i<argc;i++){
		if(argv[i][0]=='-'){
			std::string tmp(argv[i]);
			var.push_back(tmp.substr(1));
			val.push_back(argv[i+1]);
			i++;
		} else {
			std::stringstream ss;
			ss<<j;
			var.push_back(ss.str());
			val.push_back(argv[i]);
			j++;
		}
		unused.push_back(true);
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<var.size();i++){
		if(unused[i]){ std::cerr<<"Parseur : variable "<<var[i]<<" was given as input but not used"<<std::endl;}
	}
}

